#include "ancestry.h"

// ============================================================================
// User interrupts
// ============================================================================
// Long loops call R_CheckUserInterrupt() every INTERRUPT_EVERY rows so Ctrl-C /
// Esc works. It long-jumps back to R, which is safe wherever only PROTECTed
// objects and R_alloc() memory are live (R reclaims both). Code holding a
// FILE* uses interrupt_pending() instead, so it can close the file first.
#define INTERRUPT_EVERY 1024
#define INTERRUPT_EVERY_CELLS (1 << 20)

static void check_interrupt_fn(void *dummy) {
    (void) dummy;
    R_CheckUserInterrupt();
}

static int interrupt_pending(void) {
    return !R_ToplevelExec(check_interrupt_fn, NULL);
}

// ============================================================================
// Helper: Count ancestry-genotype combinations for a variant
// ============================================================================
// Returns counts for a single variant (row):
// N1 = pt==3 & gt==2, N2 = pt==3 & gt==1
// N4 = pt==2 & gt==2, N5 = pt==2 & gt==1 (singleton case)
// N7 = pt==1 & gt==2, N8 = pt==1 & gt==1
typedef struct {
    int N1, N2, N4, N5, N7, N8;
} VariantCounts;

static VariantCounts count_variant_combinations(int *gt_ptr, int *an_ptr, int nrow, int ncol, int row) {
    VariantCounts cnt = {0, 0, 0, 0, 0, 0};
    
    for (int j = 0; j < ncol; j++) {
        int gt = gt_ptr[row + nrow * j];
        int an = an_ptr[row + nrow * j];
        
        if (an == 3 && gt == 2) cnt.N1++;
        else if (an == 3 && gt == 1) cnt.N2++;
        else if (an == 2 && gt == 2) cnt.N4++;
        else if (an == 2 && gt == 1) cnt.N5++;
        else if (an == 1 && gt == 2) cnt.N7++;
        else if (an == 1 && gt == 1) cnt.N8++;
    }
    
    return cnt;
}

// ============================================================================
// count_ancestry_codes: Count occurrences of a code per row
// ============================================================================
static SEXP count_ancestry_codes_c(SEXP mat, SEXP code) {
    SEXP mat_int = PROTECT(coerceVector(mat, INTSXP));
    SEXP dim = PROTECT(getAttrib(mat_int, R_DimSymbol));
    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];
    int target_code = INTEGER(code)[0];
    
    SEXP result = PROTECT(allocVector(INTSXP, nrow));
    int *res_ptr = INTEGER(result);
    int *mat_ptr = INTEGER(mat_int);
    
    for (int i = 0; i < nrow; i++) {
        if (i % INTERRUPT_EVERY == 0) R_CheckUserInterrupt();
        res_ptr[i] = 0;
        for (int j = 0; j < ncol; j++) {
            if (mat_ptr[i + nrow * j] == target_code) {
                res_ptr[i]++;
            }
        }
    }
    
    UNPROTECT(3);
    return result;
}

// ============================================================================
// split_by_ancestry: Split genotype matrix by ancestry with proper p1/p2
// ============================================================================
// gla: 2 x 2 numeric matrix (rows: arm 0=p/1=q, cols: 0=AFR/1=EUR) of
//   per-arm global local ancestry proportions, or a zero-length vector to
//   disable GLA shrinkage entirely (reproduces the pre-shrinkage behavior
//   bit-for-bit). arm_id: integer vector (length nrow) giving each variant's
//   arm index into `gla`'s rows; ignored when gla is empty.
static SEXP split_by_ancestry_c(SEXP gt_genotype, SEXP ancestry, SEXP gla, SEXP arm_id) {
    SEXP gt_int = PROTECT(coerceVector(gt_genotype, INTSXP));
    SEXP an_int = PROTECT(coerceVector(ancestry, INTSXP));
    SEXP dim = PROTECT(getAttrib(gt_int, R_DimSymbol));

    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];

    int have_gla = length(gla) > 0;
    SEXP gla_real = PROTECT(coerceVector(gla, REALSXP));
    SEXP arm_int  = PROTECT(coerceVector(arm_id, INTSXP));
    double *gla_ptr = REAL(gla_real);
    int *arm_ptr    = INTEGER(arm_int);
    int n_arms      = have_gla ? INTEGER(getAttrib(gla_real, R_DimSymbol))[0] : 0;

    // Allocate output matrices
    SEXP african = PROTECT(allocMatrix(REALSXP, nrow, ncol));
    SEXP european = PROTECT(allocMatrix(REALSXP, nrow, ncol));

    int *gt_ptr = INTEGER(gt_int);
    int *an_ptr = INTEGER(an_int);
    double *afr_ptr = REAL(african);
    double *eur_ptr = REAL(european);

    // Process each variant (row)
    for (int i = 0; i < nrow; i++) {
        if (i % INTERRUPT_EVERY == 0) R_CheckUserInterrupt();
        // Count ancestry-genotype combinations for this variant
        VariantCounts cnt = count_variant_combinations(gt_ptr, an_ptr, nrow, ncol, i);

        // Raw (unshrunk) p1/p2:
        // p1 = (2*N1 + N2 + N4) / (2*N1 + N2 + 2*N4 + 2*N7 + N8)
        // p2 = (N4 + 2*N7 + N8) / (2*N1 + N2 + 2*N4 + 2*N7 + N8)
        double total_alt = 2 * cnt.N1 + cnt.N2 + 2 * cnt.N4 + cnt.N5 + 2 * cnt.N7 + cnt.N8;
        double denominator = total_alt - cnt.N5;   // excludes ambiguous mixed hets (N5)
        double p1_raw, p2_raw;
        if (denominator > 0.0) {
            p1_raw = (2.0 * cnt.N1 + cnt.N2 + cnt.N4) / denominator;
            p2_raw = (cnt.N4 + 2.0 * cnt.N7 + cnt.N8) / denominator;
        } else {
            p1_raw = 0.0;
            p2_raw = 0.0;
        }

        double p1, p2;
        if (!have_gla) {
            // No GLA supplied: reproduce the original (pre-shrinkage) behavior
            // exactly, including its hardcoded 0.5/0.5 singleton fallback.
            if (cnt.N5 > 0 && total_alt == (double)cnt.N5) {
                p1 = 0.5;
                p2 = 0.5;
            } else if (cnt.N5 != 0 && denominator > 0.0) {
                p1 = p1_raw;
                p2 = p2_raw;
            } else {
                p1 = 0.0;
                p2 = 0.0;
            }
        } else {
            // GLA-shrinkage: w is the fraction of alt-allele confidence
            // (not carrier count) that is ambiguous. Hom-alt carriers
            // (N1, N4, N7) carry 2 alt alleles' worth of ancestry evidence
            // vs. 1 for het carriers, so they count 2x here too -- this
            // matches the 2x weighting already applied to them in the
            // p1_raw/p2_raw numerator above. total_alt is already exactly
            // this weighted total (2*N1+N2+2*N4+N5+2*N7+N8). The more
            // ambiguous evidence dominates, the more p1/p2 shrinks toward
            // this variant's arm-level global ancestry proportion. w == 1
            // (all evidence ambiguous) reduces to p1 = GLA exactly,
            // subsuming the old singleton special case.
            if (total_alt > 0.0) {
                double w = (double) cnt.N5 / total_alt;
                int arm = arm_ptr[i];
                double gla_p1 = gla_ptr[arm + n_arms * 0];
                double gla_p2 = gla_ptr[arm + n_arms * 1];
                p1 = (1.0 - w) * p1_raw + w * gla_p1;
                p2 = (1.0 - w) * p2_raw + w * gla_p2;
            } else {
                p1 = 0.0;
                p2 = 0.0;
            }
        }

        // Apply splitting to each sample
        for (int j = 0; j < ncol; j++) {
            int gt = gt_ptr[i + nrow * j];
            int an = an_ptr[i + nrow * j];
            int idx = i + nrow * j;
            
            if (an == 3) {
                // Pure African (EUR/EUR = 0)
                afr_ptr[idx] = gt;  // All alt alleles from African
                eur_ptr[idx] = 0.0;
            } else if (an == 1) {
                // Pure European
                afr_ptr[idx] = 0.0;
                eur_ptr[idx] = gt;  // All alt alleles from European
            } else if (an == 2) {
                // Mixed ancestry: split based on p1/p2
                if (gt == 2) {
                    // Homozygous alt: 1 allele to each ancestry
                    afr_ptr[idx] = 1.0;
                    eur_ptr[idx] = 1.0;
                } else if (gt == 1) {
                    // Heterozygous: split by p1/p2
                    afr_ptr[idx] = p1;
                    eur_ptr[idx] = p2;
                } else {
                    // Homozygous ref
                    afr_ptr[idx] = 0.0;
                    eur_ptr[idx] = 0.0;
                }
            } else {
                // Invalid/missing ancestry code
                afr_ptr[idx] = 0.0;
                eur_ptr[idx] = 0.0;
            }
        }
    }
    
    // Create named list result
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, african);
    SET_VECTOR_ELT(result, 1, european);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("african"));
    SET_STRING_ELT(names, 1, mkChar("european"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(9);
    return result;
}

// ============================================================================
// read_bed_file: Read PLINK binary files
// ============================================================================
static SEXP read_bed_file_c(SEXP bed_path, SEXP bim_path, SEXP fam_path) {
    const char *bed = CHAR(STRING_ELT(bed_path, 0));
    const char *bim = CHAR(STRING_ELT(bim_path, 0));
    const char *fam = CHAR(STRING_ELT(fam_path, 0));
    
    FILE *fp = fopen(bed, "rb");
    if (!fp) error("Cannot open bed file: %s", bed);
    
    unsigned char magic[3];
    if (fread(magic, 1, 3, fp) != 3) {
        fclose(fp);
        error("Cannot read bed file header");
    }
    if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
        fclose(fp);
        error("Invalid bed file format");
    }
    
    int n_samples = 0;
    fp = fopen(fam, "r");
    if (!fp) error("Cannot open fam file: %s", fam);
    char line[1000];
    while (fgets(line, sizeof(line), fp)) n_samples++;
    fclose(fp);
    
    int n_variants = 0;
    fp = fopen(bim, "r");
    if (!fp) error("Cannot open bim file: %s", bim);
    while (fgets(line, sizeof(line), fp)) n_variants++;
    fclose(fp);
    
    int store_size = (n_samples + 3) / 4;

    // Rows = variants, columns = samples (matches split_by_ancestry_C's
    // expected orientation and count_ancestry_codes_C's "rows = regions"
    // convention elsewhere in this file).
    SEXP result = PROTECT(allocVector(INTSXP, n_samples * n_variants));
    int *gt = INTEGER(result);
    memset(gt, 0, n_samples * n_variants * sizeof(int));

    fp = fopen(bed, "rb");
    fseek(fp, 3, SEEK_SET);

    // PLINK .bed 2-bit codes (low bits = first sample in the byte):
    //   00 = homozygous first allele, 01 = missing,
    //   10 = heterozygous,            11 = homozygous second allele.
    // Remapped here to match snpStats::read.plink()'s raw numeric scheme
    // (0 = missing, 1 = hom. first allele, 2 = het, 3 = hom. second allele)
    // -- the convention this project's ancestry-encoded .bed files (pt
    // code 1=EUR/EUR, 2=AFR/EUR, 3=AFR/AFR, 0=no call) were built against.
    static const int code_map[4] = {1, 0, 2, 3};

    for (int v = 0; v < n_variants; v++) {
        if (v % INTERRUPT_EVERY == 0 && interrupt_pending()) {
            fclose(fp);
            error("Interrupted while reading bed file");
        }
        unsigned char buffer[store_size];
        if (fread(buffer, 1, store_size, fp) != (size_t)store_size) {
            fclose(fp);
            error("Truncated bed file at variant %d", v);
        }

        for (int s = 0; s < n_samples; s++) {
            int byte_idx = s / 4;
            int bit_idx = (s % 4) * 2;
            int raw = (buffer[byte_idx] >> bit_idx) & 0x03;

            gt[v + n_variants * s] = code_map[raw];
        }
    }
    fclose(fp);

    SEXP dim = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim)[0] = n_variants;
    INTEGER(dim)[1] = n_samples;
    setAttrib(result, R_DimSymbol, dim);

    UNPROTECT(2);
    return result;
}

// ============================================================================
// split_by_ancestry_multi: K-population unphased dosage splitting
// ============================================================================
//
// Algorithm (per variant row):
//   Step 1 — estimate population proportions p[k] from unambiguous observations:
//     * pure pop k, gt=2 → num[k] += 2, D += 2
//     * pure pop k, gt=1 → num[k] += 1, D += 1
//     * mixed pair (a,b), gt=2 → num[a] += 1, num[b] += 1, D += 2
//     * mixed pair (a,b), gt=1 → excluded (ambiguous)
//     p[k] = num[k] / D    (if D==0, all p[k]=0; handled per-pair below)
//
//   Step 2 — assign dosages:
//     * pure pop k, gt=g  → out[k] = g
//     * mixed (a,b), gt=2 → out[a] = 1, out[b] = 1
//     * mixed (a,b), gt=1 → out[a] = p[a]/(p[a]+p[b]), out[b] = p[b]/(p[a]+p[b])
//                           (singleton fallback: 0.5/0.5 when p[a]+p[b]==0)
//     * unknown code or gt=0 → 0
//
// Arguments:
//   gt          int matrix nrow×ncol, values 0/1/2
//   ancestry    int matrix nrow×ncol, same dims
//   pure_codes  int vector length K: pure-ancestry code for each population
//   m_code      int vector length M: mixed-ancestry codes
//   m_pop1      int vector length M: 0-based index of parent pop 1 in pure_codes
//   m_pop2      int vector length M: 0-based index of parent pop 2 in pure_codes
//   gla         n_arms x K numeric matrix of per-arm global local ancestry
//               proportions, or zero-length to disable GLA shrinkage
//               (reproduces the pre-shrinkage behavior bit-for-bit)
//   arm_id      int vector length nrow: each variant's arm index into gla's
//               rows; ignored when gla is empty
static SEXP split_by_ancestry_multi_c(SEXP gt, SEXP ancestry,
                                        SEXP pure_codes,
                                        SEXP m_code, SEXP m_pop1, SEXP m_pop2,
                                        SEXP gla, SEXP arm_id) {
    int K = length(pure_codes);
    int M = length(m_code);
    int *pure = INTEGER(pure_codes);
    int *mc   = INTEGER(m_code);
    int *mp1  = INTEGER(m_pop1);
    int *mp2  = INTEGER(m_pop2);

    SEXP dim  = PROTECT(getAttrib(gt, R_DimSymbol));
    int nrow  = INTEGER(dim)[0];   /* variants */
    int ncol  = INTEGER(dim)[1];   /* samples  */
    int n     = nrow * ncol;

    SEXP gt_i  = PROTECT(coerceVector(gt,       INTSXP));
    SEXP anc_i = PROTECT(coerceVector(ancestry, INTSXP));
    int *gp    = INTEGER(gt_i);
    int *ap    = INTEGER(anc_i);

    int have_gla = length(gla) > 0;
    SEXP gla_real = PROTECT(coerceVector(gla, REALSXP));
    SEXP arm_int  = PROTECT(coerceVector(arm_id, INTSXP));
    double *gla_ptr = REAL(gla_real);
    int *arm_ptr    = INTEGER(arm_int);
    int n_arms      = have_gla ? INTEGER(getAttrib(gla_real, R_DimSymbol))[0] : 0;

    /* Allocate K output matrices */
    SEXP result = PROTECT(allocVector(VECSXP, K));
    for (int p = 0; p < K; p++) {
        SEXP mat = PROTECT(allocMatrix(REALSXP, nrow, ncol));
        memset(REAL(mat), 0, n * sizeof(double));
        SET_VECTOR_ELT(result, p, mat);
        UNPROTECT(1);   /* protected via result */
    }
    double **out = (double **) R_alloc(K, sizeof(double *));
    for (int p = 0; p < K; p++)
        out[p] = REAL(VECTOR_ELT(result, p));

    /* Per-variant working buffers (freed automatically by R_alloc) */
    double *num = (double *) R_alloc(K, sizeof(double));
    double *pk  = (double *) R_alloc(K, sizeof(double));
    int    *amb = (int *)    R_alloc(M, sizeof(int));   /* ambiguous-het count per mixed pair */

    for (int i = 0; i < nrow; i++) {
        if (i % INTERRUPT_EVERY == 0) R_CheckUserInterrupt();

        /* -- Step 1: accumulate unambiguous alt allele counts -- */
        double D = 0.0;
        memset(num, 0, K * sizeof(double));
        memset(amb, 0, M * sizeof(int));

        for (int j = 0; j < ncol; j++) {
            int g = gp[i + nrow * j];
            int a = ap[i + nrow * j];
            if (g == 0) continue;

            /* pure ancestry */
            int hit = 0;
            for (int p = 0; p < K; p++) {
                if (a == pure[p]) {
                    num[p] += g;  D += g;
                    hit = 1;  break;
                }
            }
            if (hit) continue;

            /* mixed ancestry — only hom-alts inform p[k] */
            if (g == 2) {
                for (int m = 0; m < M; m++) {
                    if (a == mc[m]) {
                        num[mp1[m]] += 1.0;
                        num[mp2[m]] += 1.0;
                        D           += 2.0;
                        break;
                    }
                }
            } else {
                /* mixed gt==1: ambiguous, excluded from Step 1's num/D but
                   tracked per-pair for the GLA shrinkage weight below */
                for (int m = 0; m < M; m++) {
                    if (a == mc[m]) { amb[m]++; break; }
                }
            }
        }

        /* Allele-weighted confidence total for the GLA shrinkage weight
           below: D already counts hom-alt carriers 2x vs. het (matching
           the p[k] numerator above); each ambiguous mixed-het contributes
           1 unresolved allele of confidence, same as N5 in the
           2-population case (split_by_ancestry_c). */
        double total_amb = 0.0;
        for (int m = 0; m < M; m++) total_amb += amb[m];
        double total_weighted = D + total_amb;

        /* -- Compute p[k] -- */
        if (D > 0.0)
            for (int p = 0; p < K; p++) pk[p] = num[p] / D;
        else
            memset(pk, 0, K * sizeof(double));  /* singleton case handled per-pair */

        int arm = have_gla ? arm_ptr[i] : 0;

        /* -- Step 2: assign dosages -- */
        for (int j = 0; j < ncol; j++) {
            int g   = gp[i + nrow * j];
            int a   = ap[i + nrow * j];
            int idx = i + nrow * j;
            if (g == 0) continue;

            /* pure ancestry */
            int hit = 0;
            for (int p = 0; p < K; p++) {
                if (a == pure[p]) {
                    out[p][idx] = (double) g;
                    hit = 1;  break;
                }
            }
            if (hit) continue;

            /* mixed ancestry */
            for (int m = 0; m < M; m++) {
                if (a == mc[m]) {
                    int pa = mp1[m], pb = mp2[m];
                    if (g == 2) {
                        out[pa][idx] = 1.0;
                        out[pb][idx] = 1.0;
                    } else if (!have_gla) {
                        /* No GLA supplied: reproduce original behavior
                           exactly, including its 0.5/0.5 singleton fallback. */
                        double sum_ab = pk[pa] + pk[pb];
                        if (sum_ab > 0.0) {
                            out[pa][idx] = pk[pa] / sum_ab;
                            out[pb][idx] = pk[pb] / sum_ab;
                        } else {
                            out[pa][idx] = 0.5;
                            out[pb][idx] = 0.5;
                        }
                    } else {
                        /* GLA shrinkage: w is the fraction of alt-allele
                           confidence (not carrier count) that is ambiguous
                           for this specific pair -- hom-alt carriers count
                           2x vs. het, matching total_weighted's derivation
                           from D above. Raw pair ratio collapses to the GLA
                           conditional fraction when this pair has no
                           unambiguous hom-alt evidence of its own, so
                           blending is well-defined even at w < 1. */
                        double gla_pa = gla_ptr[arm + n_arms * pa];
                        double gla_pb = gla_ptr[arm + n_arms * pb];
                        double gla_sum = gla_pa + gla_pb;
                        double gla_frac_pa = (gla_sum > 0.0) ? gla_pa / gla_sum : 0.5;

                        double sum_ab = pk[pa] + pk[pb];
                        double pa_raw = (sum_ab > 0.0) ? pk[pa] / sum_ab : gla_frac_pa;

                        double w = (double) amb[m] / total_weighted;
                        double blended_pa = (1.0 - w) * pa_raw + w * gla_frac_pa;
                        out[pa][idx] = blended_pa;
                        out[pb][idx] = 1.0 - blended_pa;
                    }
                    break;
                }
            }
            /* unknown ancestry code: leave as 0 */
        }
    }

    UNPROTECT(6);   /* dim, gt_i, anc_i, gla_real, arm_int, result */
    return result;
}

// ============================================================================
// split_phased_by_ancestry: Deterministic split using per-haplotype ancestry
// ============================================================================
// Each haplotype contributes its allele to whichever ancestry population it
// was called as.  NA genotypes are treated as 0 (reference).  Any ancestry
// code that is neither AFR nor EUR contributes nothing to either output.
//
// pop_codes: integer vector of length 2 — pop_codes[0] = AFR code,
//            pop_codes[1] = EUR code (RFMix default: 0 and 1).
static SEXP split_phased_by_ancestry_c(SEXP gt_hap0, SEXP gt_hap1,
                                        SEXP anc_hap0, SEXP anc_hap1,
                                        SEXP pop_codes) {
    SEXP dim = PROTECT(getAttrib(gt_hap0, R_DimSymbol));
    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];
    int afr_code = INTEGER(pop_codes)[0];
    int eur_code = INTEGER(pop_codes)[1];

    SEXP gh0 = PROTECT(coerceVector(gt_hap0, INTSXP));
    SEXP gh1 = PROTECT(coerceVector(gt_hap1, INTSXP));
    SEXP ah0 = PROTECT(coerceVector(anc_hap0, INTSXP));
    SEXP ah1 = PROTECT(coerceVector(anc_hap1, INTSXP));

    SEXP african  = PROTECT(allocMatrix(REALSXP, nrow, ncol));
    SEXP european = PROTECT(allocMatrix(REALSXP, nrow, ncol));

    int *gh0_ptr = INTEGER(gh0);
    int *gh1_ptr = INTEGER(gh1);
    int *ah0_ptr = INTEGER(ah0);
    int *ah1_ptr = INTEGER(ah1);
    double *afr_ptr = REAL(african);
    double *eur_ptr = REAL(european);

    int n = nrow * ncol;
    for (int k = 0; k < n; k++) {
        if (k % INTERRUPT_EVERY_CELLS == 0) R_CheckUserInterrupt();
        int g0 = (gh0_ptr[k] == NA_INTEGER) ? 0 : gh0_ptr[k];
        int g1 = (gh1_ptr[k] == NA_INTEGER) ? 0 : gh1_ptr[k];
        int a0 = ah0_ptr[k];
        int a1 = ah1_ptr[k];

        double afr = 0.0, eur = 0.0;
        if (a0 == afr_code)      afr += g0;
        else if (a0 == eur_code) eur += g0;

        if (a1 == afr_code)      afr += g1;
        else if (a1 == eur_code) eur += g1;

        afr_ptr[k] = afr;
        eur_ptr[k] = eur;
    }

    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, african);
    SET_VECTOR_ELT(result, 1, european);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("african"));
    SET_STRING_ELT(names, 1, mkChar("european"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(9);
    return result;
}

// ============================================================================
// split_phased_multi: K-population phased dosage splitting
// ============================================================================
// pop_codes: named integer vector of length K.
//   pop_codes[p] is the ancestry code for population p.
//   The output list has K matrices, named after pop_codes.
// Each haplotype allele is routed to the pool whose code matches the
// haplotype's local ancestry call.  Unrecognised codes contribute 0
// to all pools.  NA genotypes are treated as 0 (reference allele).
static SEXP split_phased_multi_c(SEXP gt_hap0, SEXP gt_hap1,
                                   SEXP anc_hap0, SEXP anc_hap1,
                                   SEXP pop_codes) {
    int K = length(pop_codes);
    int *codes = INTEGER(pop_codes);
    SEXP code_names = getAttrib(pop_codes, R_NamesSymbol);

    SEXP dim = PROTECT(getAttrib(gt_hap0, R_DimSymbol));
    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];
    int n = nrow * ncol;

    SEXP gh0 = PROTECT(coerceVector(gt_hap0, INTSXP));
    SEXP gh1 = PROTECT(coerceVector(gt_hap1, INTSXP));
    SEXP ah0 = PROTECT(coerceVector(anc_hap0, INTSXP));
    SEXP ah1 = PROTECT(coerceVector(anc_hap1, INTSXP));

    int *gh0_ptr = INTEGER(gh0);
    int *gh1_ptr = INTEGER(gh1);
    int *ah0_ptr = INTEGER(ah0);
    int *ah1_ptr = INTEGER(ah1);

    /* Allocate result list and K output matrices */
    SEXP result = PROTECT(allocVector(VECSXP, K));
    for (int p = 0; p < K; p++) {
        SEXP mat = PROTECT(allocMatrix(REALSXP, nrow, ncol));
        memset(REAL(mat), 0, n * sizeof(double));
        SET_VECTOR_ELT(result, p, mat);
        UNPROTECT(1);   /* mat now protected via result */
    }

    /* Cache per-population output pointers to avoid repeated VECTOR_ELT calls */
    double **out = (double **) R_alloc(K, sizeof(double *));
    for (int p = 0; p < K; p++)
        out[p] = REAL(VECTOR_ELT(result, p));

    /* Route each haplotype allele to its population pool */
    for (int k = 0; k < n; k++) {
        if (k % INTERRUPT_EVERY_CELLS == 0) R_CheckUserInterrupt();
        int g0 = (gh0_ptr[k] == NA_INTEGER) ? 0 : gh0_ptr[k];
        int g1 = (gh1_ptr[k] == NA_INTEGER) ? 0 : gh1_ptr[k];
        int a0 = ah0_ptr[k];
        int a1 = ah1_ptr[k];

        for (int p = 0; p < K; p++) {
            if (a0 == codes[p]) { out[p][k] += g0; break; }
        }
        for (int p = 0; p < K; p++) {
            if (a1 == codes[p]) { out[p][k] += g1; break; }
        }
    }

    /* Name the output list after the pop_codes names */
    if (code_names != R_NilValue) {
        SEXP names = PROTECT(duplicate(code_names));
        setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1);
    }

    UNPROTECT(6);   /* dim, gh0, gh1, ah0, ah1, result */
    return result;
}

// ============================================================================
// Exported functions (registered in init.c)
// ============================================================================
SEXP count_ancestry_codes(SEXP mat, SEXP code) {
    return count_ancestry_codes_c(mat, code);
}

SEXP split_by_ancestry(SEXP gt_genotype, SEXP ancestry, SEXP gla, SEXP arm_id) {
    return split_by_ancestry_c(gt_genotype, ancestry, gla, arm_id);
}

SEXP read_bed_file(SEXP bed_path, SEXP bim_path, SEXP fam_path) {
    return read_bed_file_c(bed_path, bim_path, fam_path);
}

SEXP split_phased_multi(SEXP gt_hap0, SEXP gt_hap1,
                         SEXP anc_hap0, SEXP anc_hap1, SEXP pop_codes) {
    return split_phased_multi_c(gt_hap0, gt_hap1, anc_hap0, anc_hap1, pop_codes);
}

SEXP split_by_ancestry_multi(SEXP gt, SEXP ancestry, SEXP pure_codes,
                               SEXP m_code, SEXP m_pop1, SEXP m_pop2,
                               SEXP gla, SEXP arm_id) {
    return split_by_ancestry_multi_c(gt, ancestry, pure_codes,
                                     m_code, m_pop1, m_pop2, gla, arm_id);
}

SEXP split_phased_by_ancestry(SEXP gt_hap0, SEXP gt_hap1,
                               SEXP anc_hap0, SEXP anc_hap1,
                               SEXP pop_codes) {
    return split_phased_by_ancestry_c(gt_hap0, gt_hap1, anc_hap0, anc_hap1, pop_codes);
}
