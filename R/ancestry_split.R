# ============================================================================
# ancestry_split (Step 1): split a VCF by local ancestry
# ============================================================================
#
# Low-level C-backed split primitives (diploid/dosage + phased/haplotype),
# MSP/VCF parsing helpers, and the ancestry_split() entry point plus its
# backward-compatible ancestry_split_dosage() / ancestry_split_phased()
# wrappers.

# ============================================================================
# Internal C wrappers
# ============================================================================

.count_ancestry_codes <- function(mat, code) {
    storage.mode(mat) <- "integer"
    code <- as.integer(code)
    .Call("count_ancestry_codes_C", mat, code, PACKAGE = "lantern")
}

#' Normalize gla/arm_id into .Call-ready sentinels
#'
#' \code{gla = NULL} (either argument NULL) disables shrinkage: an empty
#' matrix/vector signal the C side to reproduce the pre-shrinkage behavior
#' exactly (see \code{ancestry.c}).
#' @keywords internal
.gla_call_args <- function(gla, arm_id) {
  if (is.null(gla) || is.null(arm_id))
    return(list(gla = matrix(numeric(0), nrow = 0, ncol = 0), arm_id = integer(0)))
  storage.mode(gla) <- "double"
  list(gla = gla, arm_id = as.integer(arm_id))
}

#' Uniform shrinkage target: 1/K for each of K populations
#'
#' One-row target matrix for \code{split_diploid()} /
#' \code{split_diploid_multi()}'s \code{gla} argument (use with an all-zero
#' \code{arm_id}). For every ambiguous pair the C code shrinks toward the
#' pair's target share \code{t_a / (t_a + t_b)}, which is exactly 1/2 here.
#' @keywords internal
.uniform_shrink_target <- function(pop_names) {
  matrix(1 / length(pop_names), nrow = 1L, ncol = length(pop_names),
         dimnames = list("all", pop_names))
}

.split_diploid <- function(gt_genotype, ancestry, gla = NULL, arm_id = NULL) {
    storage.mode(gt_genotype) <- "integer"
    storage.mode(ancestry) <- "integer"
    g <- .gla_call_args(gla, arm_id)
    .Call("split_by_ancestry_C", gt_genotype, ancestry, g$gla, g$arm_id, PACKAGE = "lantern")
}

.split_haplotype <- function(gt_hap0, gt_hap1, anc_hap0, anc_hap1,
                           pop_codes = c(AFR = 0L, EUR = 1L)) {
    storage.mode(gt_hap0)  <- "integer"
    storage.mode(gt_hap1)  <- "integer"
    storage.mode(anc_hap0) <- "integer"
    storage.mode(anc_hap1) <- "integer"
    pop_codes <- as.integer(pop_codes)
    .Call("split_phased_by_ancestry_C", gt_hap0, gt_hap1, anc_hap0, anc_hap1,
          pop_codes, PACKAGE = "lantern")
}

# ============================================================================
# Exported split primitives
# ============================================================================

#' Count ancestry codes in matrix
#'
#' Count occurrences of a specific ancestry code in each row of a matrix.
#' Internal building block used by \code{\link{ancestry_split}} (to compute
#' \code{ancestry_counts}) and \code{\link{ancestry_split_dosage}}.
#'
#' @param mat Integer matrix with ancestry codes (1, 2, 3)
#'   - Rows represent genomic regions/windows
#'   - Columns represent samples
#' @param code Ancestry code to count:
#'   - 1 = EUR/EUR (Pure European)
#'   - 2 = AFR/EUR (Mixed)
#'   - 3 = AFR/AFR (Pure African)
#'
#' @return Integer vector with count of \code{code} per row
#'
#' @keywords internal
count_ancestry_codes <- function(mat, code) {
    .count_ancestry_codes(mat, code)
}

#' Split genotype matrix by ancestry
#'
#' Split a genotype matrix into African and European ancestry-specific
#' dosage matrices based on a parent-of-origin ancestry matrix.
#'
#' @param gt_genotype Integer matrix of genotypes (0, 1, 2)
#'   - Rows represent variants
#'   - Columns represent samples
#' @param ancestry Integer matrix of ancestry codes (1, 2, 3)
#'   - Same dimensions as gt_genotype
#'   - Values: 1=EUR/EUR, 2=AFR/EUR (mixed), 3=AFR/AFR
#'
#' @return List with two elements:
#'   \item{african}{Matrix of African ancestry-specific dosages}
#'   \item{european}{Matrix of European ancestry-specific dosages}
#'
#' @examples
#' # Genotype matrix: 5 variants x 4 samples
#' gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1, 0,
#'                1, 1, 1, 0, 2, 1, 0, 1, 0, 1), nrow = 5, ncol = 4)
#'
#' # Ancestry matrix: same dimensions
#' ancestry <- matrix(c(3, 2, 1, 3, 2, 1, 2, 2, 1, 1,
#'                      3, 1, 2, 1, 3, 2, 1, 2, 1, 3), nrow = 5, ncol = 4)
#'
#' result <- split_diploid(gt, ancestry)
#' result$african
#' result$european
#'
#' @param gla Optional shrinkage target: numeric matrix with one row per
#'   target group and two columns (AFR, EUR) of target proportions. A mixed
#'   heterozygote's split is shrunk toward its row's AFR share by the
#'   ambiguous share \eqn{w} of the variant's evidence.
#'   \code{\link{ancestry_split}} uses a single row of 1/2. When \code{NULL}
#'   (default), no shrinkage: the raw p1/p2 ratio is used, with a 0.5/0.5
#'   fallback only when every carrier is ambiguous.
#' @param arm_id Optional integer vector (length = \code{nrow(gt_genotype)})
#'   giving each variant's row index into \code{gla} (0-based; all zeros for
#'   a one-row target). Required together with \code{gla}.
#'
#' @export
split_diploid <- function(gt_genotype, ancestry, gla = NULL, arm_id = NULL) {
    .split_diploid(gt_genotype, ancestry, gla, arm_id)
}

#' Split unphased genotype matrix by ancestry for K populations
#'
#' Generalises \code{\link{split_diploid}} to an arbitrary number of
#' ancestries.  For each variant, population allele proportions
#' \eqn{p_1, \ldots, p_K} (summing to 1) are estimated from samples with
#' unambiguous ancestry (pure-ancestry homozygotes and heterozygotes, plus
#' mixed hom-alts which each contribute one allele to each parent pool).
#' Ambiguous heterozygotes in each mixed pair are then split by the
#' conditional pairwise ratio \eqn{p_i / (p_i + p_j)}.
#' When no unambiguous alt alleles exist for a pair (singleton edge case),
#' the split defaults to 0.5 / 0.5.
#'
#' @param gt_genotype Integer matrix (variants x samples) of genotypes (0/1/2).
#' @param ancestry Integer matrix (variants x samples) of ancestry codes,
#'   same dimensions as \code{gt_genotype}.
#' @param pure_codes Named integer vector mapping each population label to its
#'   pure-ancestry diploid code.
#'   Example: \code{c(AFR = 3L, EUR = 1L, NAT = 4L)}.
#' @param mixed_codes Data frame with three columns:
#'   \describe{
#'     \item{code}{Integer ancestry code for this mixed pair.}
#'     \item{pop1}{Character name of the first parent population
#'       (must be a name in \code{pure_codes}).}
#'     \item{pop2}{Character name of the second parent population
#'       (must be a name in \code{pure_codes}).}
#'   }
#'   One row per mixed-ancestry diploid type.
#'   Example for 3 populations:
#'   \code{data.frame(code=c(2L,5L,6L), pop1=c("AFR","AFR","EUR"), pop2=c("EUR","NAT","NAT"))}
#'
#' @return Named list of K numeric matrices (variants x samples), one per
#'   population in \code{pure_codes}, in the same order.  Each element is
#'   named after the corresponding entry in \code{pure_codes}.
#'
#' @seealso \code{\link{split_diploid}} for the hardcoded 2-population
#'   (AFR + EUR) wrapper; \code{\link{split_haplotype_multi}} for the phased
#'   K-population analogue.
#'
#' @examples
#' # 3-population example: AFR=3, EUR=1, NAT=4; mixed codes 2/5/6
#' set.seed(1)
#' gt  <- matrix(sample(0:2, 20, replace = TRUE, prob = c(.8,.15,.05)),
#'               nrow = 4, ncol = 5)
#' # ancestry codes: 3=AFR/AFR, 1=EUR/EUR, 4=NAT/NAT,
#' #                 2=AFR/EUR, 5=AFR/NAT, 6=EUR/NAT
#' anc <- matrix(sample(c(1L,2L,3L,4L,5L,6L), 20, replace = TRUE),
#'               nrow = 4, ncol = 5)
#' pure  <- c(AFR = 3L, EUR = 1L, NAT = 4L)
#' mixed <- data.frame(code = c(2L, 5L, 6L),
#'                     pop1 = c("AFR", "AFR", "EUR"),
#'                     pop2 = c("EUR", "NAT", "NAT"))
#' out <- split_diploid_multi(gt, anc, pure, mixed)
#' names(out)   # "AFR" "EUR" "NAT"
#'
#' @param gla Optional shrinkage target: numeric matrix with one row per
#'   target group and columns matching \code{names(pure_codes)}. Each
#'   ambiguous pair (a, b) is shrunk toward \eqn{t_a / (t_a + t_b)} by the
#'   pair's ambiguous share \eqn{w} of the evidence.
#'   \code{\link{ancestry_split}} uses a single row of 1/K (i.e. 1/2 for
#'   every pair). When \code{NULL} (default), no shrinkage: the raw p[k]
#'   ratio is used, with a 0.5/0.5 fallback only when a pair has no
#'   unambiguous evidence.
#' @param arm_id Optional integer vector (length = \code{nrow(gt_genotype)})
#'   giving each variant's row index into \code{gla} (0-based; all zeros for
#'   a one-row target). Required together with \code{gla}.
#'
#' @export
split_diploid_multi <- function(gt_genotype, ancestry,
                                     pure_codes, mixed_codes,
                                     gla = NULL, arm_id = NULL) {
    if (is.null(names(pure_codes)))
        stop("pure_codes must be a named integer vector, e.g. c(AFR=3L, EUR=1L)")
    if (!is.data.frame(mixed_codes) ||
        !all(c("code", "pop1", "pop2") %in% names(mixed_codes)))
        stop("mixed_codes must be a data.frame with columns: code, pop1, pop2")
    pop_names <- names(pure_codes)
    bad <- setdiff(c(mixed_codes$pop1, mixed_codes$pop2), pop_names)
    if (length(bad))
        stop("mixed_codes references populations not in pure_codes: ",
             paste(bad, collapse = ", "))
    if (!is.null(gla) && !identical(colnames(gla), pop_names))
        stop("gla's columns must match names(pure_codes) exactly, in order: ",
             paste(pop_names, collapse = ", "))

    storage.mode(gt_genotype) <- "integer"
    storage.mode(ancestry)    <- "integer"
    storage.mode(pure_codes)  <- "integer"

    m_code <- as.integer(mixed_codes$code)
    m_pop1 <- match(mixed_codes$pop1, pop_names) - 1L   # 0-based index
    m_pop2 <- match(mixed_codes$pop2, pop_names) - 1L
    g <- .gla_call_args(gla, arm_id)

    result <- .Call("split_by_ancestry_multi_C",
                    gt_genotype, ancestry,
                    pure_codes, m_code, m_pop1, m_pop2,
                    g$gla, g$arm_id,
                    PACKAGE = "lantern")
    names(result) <- pop_names
    result
}

#' Split phased haplotypes by ancestry
#'
#' Split phased haplotypes deterministically into African and European
#' ancestry-specific dosage matrices.  Each haplotype's allele is added
#' to whichever ancestry pool matches its local ancestry call; haplotypes
#' with unrecognised ancestry codes contribute nothing to either pool.
#'
#' @param gt_hap0 Integer matrix (variants x samples) of haplotype-0 alleles (0/1).
#' @param gt_hap1 Integer matrix (variants x samples) of haplotype-1 alleles (0/1).
#' @param anc_hap0 Integer matrix (variants x samples) of haplotype-0 ancestry codes.
#' @param anc_hap1 Integer matrix (variants x samples) of haplotype-1 ancestry codes.
#' @param pop_codes Named integer vector of length 2 giving the AFR and EUR
#'   ancestry codes as used in the MSP file.  Defaults to
#'   \code{c(AFR = 0L, EUR = 1L)} (RFMix convention).
#'
#' @return List with two numeric matrices:
#'   \item{african}{African ancestry-specific dosage (variants x samples)}
#'   \item{european}{European ancestry-specific dosage (variants x samples)}
#'
#' @examples
#' gt0 <- matrix(c(1L, 0L, 0L, 1L), 2, 2)
#' gt1 <- matrix(c(0L, 1L, 1L, 0L), 2, 2)
#' a0  <- matrix(c(0L, 1L, 0L, 1L), 2, 2)
#' a1  <- matrix(c(1L, 0L, 1L, 0L), 2, 2)
#' split_haplotype(gt0, gt1, a0, a1)
#'
#' @export
split_haplotype <- function(gt_hap0, gt_hap1, anc_hap0, anc_hap1,
                         pop_codes = c(AFR = 0L, EUR = 1L)) {
    .split_haplotype(gt_hap0, gt_hap1, anc_hap0, anc_hap1, pop_codes)
}

#' Split phased haplotypes into K population-specific dosage matrices
#'
#' Generalises \code{\link{split_haplotype}} to an arbitrary number of ancestries.
#' Each haplotype's allele is routed to the dosage pool whose code matches its
#' local ancestry call.  Haplotypes with unrecognised codes contribute 0 to
#' all pools.  NA genotypes are treated as 0 (reference allele).
#'
#' For the two-population case \code{split_haplotype_multi} and \code{split_haplotype}
#' are equivalent; prefer \code{split_haplotype_multi} for new code that may later
#' be extended to three or more ancestries.
#'
#' @param gt_hap0 Integer matrix (variants x samples) of haplotype-0 alleles (0/1).
#' @param gt_hap1 Integer matrix (variants x samples) of haplotype-1 alleles (0/1).
#' @param anc_hap0 Integer matrix (variants x samples) of haplotype-0 ancestry codes.
#' @param anc_hap1 Integer matrix (variants x samples) of haplotype-1 ancestry codes.
#' @param pop_codes Named integer vector mapping population label to ancestry code,
#'   e.g. \code{c(AFR = 0L, EUR = 1L, NAT = 2L)}.  Codes must match those used
#'   in the MSP file (RFMix default: 0-based integers).
#'
#' @return Named list of K numeric matrices (variants x samples), one per
#'   entry in \code{pop_codes}, in the same order.  List names equal the
#'   names of \code{pop_codes}.
#'
#' @seealso \code{\link{split_haplotype}} for the two-population convenience wrapper,
#'   \code{\link{split_diploid}} for unphased (diploid-code) splitting.
#'
#' @examples
#' # Three-population panel: AFR=0, EUR=1, NAT=2
#' gt0 <- matrix(c(1L, 0L, 0L, 1L, 1L, 0L), nrow = 3, ncol = 2)
#' gt1 <- matrix(c(0L, 1L, 1L, 0L, 0L, 1L), nrow = 3, ncol = 2)
#' a0  <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 3, ncol = 2)
#' a1  <- matrix(c(1L, 2L, 0L, 2L, 0L, 1L), nrow = 3, ncol = 2)
#' out <- split_haplotype_multi(gt0, gt1, a0, a1, c(AFR = 0L, EUR = 1L, NAT = 2L))
#' names(out)   # "AFR" "EUR" "NAT"
#' out$NAT      # native-ancestry dosage matrix
#'
#' @export
split_haplotype_multi <- function(gt_hap0, gt_hap1, anc_hap0, anc_hap1,
                                pop_codes = c(AFR = 0L, EUR = 1L)) {
    if (is.null(names(pop_codes)))
        stop("pop_codes must be a named integer vector, e.g. c(AFR=0L, EUR=1L, NAT=2L)")
    if (length(pop_codes) < 2L)
        stop("pop_codes must contain at least 2 populations")
    pop_names <- names(pop_codes)
    storage.mode(gt_hap0)  <- "integer"
    storage.mode(gt_hap1)  <- "integer"
    storage.mode(anc_hap0) <- "integer"
    storage.mode(anc_hap1) <- "integer"
    storage.mode(pop_codes) <- "integer"   # preserves names unlike as.integer()
    result <- .Call("split_phased_multi_C", gt_hap0, gt_hap1, anc_hap0, anc_hap1,
                    pop_codes, PACKAGE = "lantern")
    names(result) <- pop_names             # guarantee names even if C skips them
    result
}

# ============================================================================
# PLINK BED reader (Step 1 input helper for BED-encoded local ancestry)
# ============================================================================

#' Read a PLINK .bed file of BED-encoded local ancestry calls
#'
#' Reads a PLINK binary genotype file (\code{.bed}/\code{.bim}/\code{.fam}
#' trio) into an integer matrix, for the case where local ancestry has been
#' encoded as PLINK genotype calls at ancestry-tract "SNPs" (as consumed by
#' \code{src/step1_vcf_split_by_ancestry.R}): homozygous first allele = pure
#' EUR/EUR (code 1), heterozygous = mixed AFR/EUR (code 2), homozygous
#' second allele = pure AFR/AFR (code 3), missing = no call (code 0). This
#' matches \code{snpStats::read.plink()}'s raw numeric convention (values
#' feed directly into \code{\link{split_diploid}}/\code{\link{ancestry_split_dosage}}
#' as \code{pt_matrix}/\code{ancestry} without \code{snpStats} as a
#' dependency.
#'
#' @param bed Path to the \code{.bed} file.
#' @param bim Path to the \code{.bim} file (variant IDs are taken from
#'   column 2).
#' @param fam Path to the \code{.fam} file (sample IDs are taken from
#'   column 2, i.e. IID).
#'
#' @return Integer matrix (variants x samples), dimnamed from the
#'   \code{.bim}/\code{.fam} files. Values: 0 = no call, 1 = homozygous
#'   first allele, 2 = heterozygous, 3 = homozygous second allele.
#'
#' @seealso \code{\link{split_diploid}}, \code{\link{ancestry_split_dosage}}
#'
#' @examples
#' ext <- system.file("extdata", package = "lantern")
#' pt <- read_bed_file(file.path(ext, "toy_ancestry.bed"),
#'                     file.path(ext, "toy_ancestry.bim"),
#'                     file.path(ext, "toy_ancestry.fam"))
#' pt[, 1:6]   # ancestry tracts x samples, codes 1/2/3
#'
#' @export
read_bed_file <- function(bed, bim, fam) {
  bim_dt <- data.table::fread(bim, header = FALSE)
  fam_dt <- data.table::fread(fam, header = FALSE)
  variant_ids <- as.character(bim_dt[[2]])
  sample_ids  <- as.character(fam_dt[[2]])

  mat <- .Call("read_bed_file_C", bed, bim, fam, PACKAGE = "lantern")
  dimnames(mat) <- list(variant_ids, sample_ids)
  mat
}

# ============================================================================
# MSP / VCF parsing helpers
# ============================================================================

#' Parse RFMix MSP file
#'
#' Parse the RFMix MSP (local ancestry tract) file into ancestry matrices.
#' MSP format: two header lines followed by data rows.
#' Line 1: population codes (e.g., "#Subpopulation order/codes: AFR=0  EUR=1").
#' Line 2: column headers.
#' Data rows: ancestry calls per haplotype per tract.
#'
#' @param msp_path Path to MSP file (plain text or gzipped).
#' @param verbose Print progress messages.
#' @return List with: \code{pop_codes} (named integer vector),
#'   \code{sample_ids} (character), \code{tract_df} (data.frame sorted by
#'   spos), \code{anc_hap0} (integer matrix: samples x tracts),
#'   \code{anc_hap1} (integer matrix: samples x tracts).
#' @keywords internal
.parse_msp <- function(msp_path, verbose = TRUE) {
  if (verbose) message("  Parsing MSP file: ", msp_path)

  con <- if (grepl("\\.gz$", msp_path)) gzfile(msp_path, "rt") else file(msp_path, "r")
  on.exit(close(con))

  line1 <- readLines(con, n = 1)
  line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
  code_parts <- strsplit(line1, if (grepl("\t", line1)) "\t" else "\\s+")[[1]]
  pop_codes <- integer()
  for (part in code_parts) {
    m <- regmatches(part, regexec("([A-Za-z0-9_]+)\\s*=\\s*([0-9]+)", part))[[1]]
    if (length(m) == 3) pop_codes[m[2]] <- as.integer(m[3])
  }
  if (length(pop_codes) == 0) {
    warning("Could not parse pop codes from MSP; using default AFR=0, EUR=1")
    pop_codes <- c(AFR = 0L, EUR = 1L)
  }

  line2      <- readLines(con, n = 1)
  headers    <- strsplit(line2, "\t")[[1]]
  sample_hap_cols <- headers[7:length(headers)]
  sample_ids <- unique(sub("\\.[01]$", "", sample_hap_cols))

  if (verbose) {
    message("    Found ", length(sample_ids), " samples, ", length(pop_codes), " populations")
    message("    Pop codes: ", paste(names(pop_codes), "=", pop_codes, collapse = ", "))
  }

  data_lines <- readLines(con)
  if (length(data_lines) == 0) stop("MSP file has no data rows")

  data_mat <- do.call(rbind, lapply(data_lines, function(l) strsplit(l, "\t")[[1]]))

  tract_df <- data.frame(
    chrom = data_mat[, 1],
    spos  = as.integer(data_mat[, 2]),
    epos  = as.integer(data_mat[, 3]),
    stringsAsFactors = FALSE
  )
  tract_order <- order(tract_df$spos)
  tract_df    <- tract_df[tract_order, ]
  rownames(tract_df) <- NULL

  # Ancestry columns start at column 7 in data rows (6 metadata fields: chm spos epos sgpos egpos "n snps")
  anc_values   <- data_mat[, 7:ncol(data_mat), drop = FALSE]
  n_tracts_orig <- nrow(tract_df)
  n_samples     <- length(sample_ids)

  anc_hap0_mat <- matrix(0L, nrow = n_samples, ncol = n_tracts_orig)
  anc_hap1_mat <- matrix(0L, nrow = n_samples, ncol = n_tracts_orig)
  rownames(anc_hap0_mat) <- rownames(anc_hap1_mat) <- sample_ids

  for (j in seq_len(n_tracts_orig)) {
    vals <- anc_values[tract_order[j], ]
    for (s in seq_len(n_samples)) {
      anc_hap0_mat[s, j] <- as.integer(vals[2 * s - 1])
      anc_hap1_mat[s, j] <- as.integer(vals[2 * s])
    }
  }

  list(pop_codes  = pop_codes,
       sample_ids = sample_ids,
       tract_df   = tract_df,
       anc_hap0   = anc_hap0_mat,
       anc_hap1   = anc_hap1_mat)
}

#' Parse phased GT matrix into haplotype integer matrices (vectorised)
#'
#' @param gt_mat Character matrix (variants x samples).
#' @param n_variants Number of variants.
#' @param n_samples Number of samples.
#' @param sample_order Column names for output matrices.
#' @param verbose Print progress.
#' @return List: \code{hap0}, \code{hap1} (integer matrices), \code{missing}
#'   (logical matrix).
#' @keywords internal
.parse_phased_gt_matrix <- function(gt_mat, n_variants, n_samples,
                                    sample_order, verbose = FALSE) {
  gt_vec  <- as.character(gt_mat)
  missing_vec <- is.na(gt_vec) | gt_vec %in% c(".", "./.", ".|.")

  unphased <- grepl("/", gt_vec, fixed = TRUE) & !grepl("|", gt_vec, fixed = TRUE)
  if (any(unphased))
    warning("Unphased GT (\"/\") encountered for ", sum(unphased),
            " cells. Treating first allele as hap0.")

  gt_norm <- gt_vec
  gt_norm[unphased] <- sub("/", "|", gt_norm[unphased], fixed = TRUE)

  nch      <- nchar(gt_norm)
  hap0_str <- ifelse(nch >= 1, substr(gt_norm, 1, 1), "0")
  hap1_str <- ifelse(nch >= 3, substr(gt_norm, 3, 3),
                     ifelse(nch >= 1, substr(gt_norm, 1, 1), "0"))
  hap0_str[missing_vec] <- "0"
  hap1_str[missing_vec] <- "0"

  hap0_vec <- suppressWarnings(as.integer(hap0_str))
  hap1_vec <- suppressWarnings(as.integer(hap1_str))
  hap0_vec[is.na(hap0_vec)] <- 0L
  hap1_vec[is.na(hap1_vec)] <- 0L

  multi <- hap0_vec > 1L | hap1_vec > 1L
  if (any(multi)) {
    warning("GT with allele >1 (multiallelic?) for ", sum(multi),
            " cells. Treating as missing (0/0).")
    hap0_vec[multi] <- 0L; hap1_vec[multi] <- 0L; missing_vec[multi] <- TRUE
  }

  hap0_mat <- matrix(hap0_vec, nrow = n_variants, ncol = n_samples)
  hap1_mat <- matrix(hap1_vec, nrow = n_variants, ncol = n_samples)
  missing_mat <- matrix(missing_vec, nrow = n_variants, ncol = n_samples)
  colnames(hap0_mat) <- colnames(hap1_mat) <- sample_order

  list(hap0 = hap0_mat, hap1 = hap1_mat, missing = missing_mat)
}

# ============================================================================
# Shared VCF + MSP parsing (steps 1-6, used by ancestry_split() below)
# ============================================================================

#' Resolve a chromosome name against a VCF/BCF's declared contigs
#'
#' Tries \code{chrom} as given, with a leading "chr" stripped, and with a
#' leading "chr" added, against the \code{##contig=<ID=...>} lines in the
#' file's header. Used to restrict \code{bcftools query} to one chromosome
#' via \code{-t} without guessing the file's naming convention. Returns
#' \code{NULL} (triggering an unrestricted query, filtered in R instead) if
#' the header has no contig lines or none of the candidates match.
#'
#' @keywords internal
.resolve_vcf_contig <- function(vcf_path, chrom) {
  header <- tryCatch(
    system2("bcftools", c("view", "-h", shQuote(vcf_path)),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character(0))
  contig_ids <- sub("^##contig=<ID=([^,>]+).*$", "\\1",
                     grep("^##contig=<ID=", header, value = TRUE))
  if (length(contig_ids) == 0) return(NULL)
  chrom_clean <- sub("^chr", "", chrom)
  candidates  <- c(chrom, chrom_clean, paste0("chr", chrom_clean))
  match_id    <- candidates[candidates %in% contig_ids]
  if (length(match_id) == 0) return(NULL)
  match_id[1]
}

#' Choose the VCF-reading backend
#'
#' \code{"bcftools"} when the \code{bcftools} binary is on \code{PATH}
#' (fastest; also reads BCF), otherwise \code{"seqarray"} (no external
#' binary). \code{options(lantern.vcf_reader = "bcftools"/"seqarray")}
#' forces one; the default \code{"auto"} picks as above.
#'
#' @keywords internal
.vcf_reader <- function() {
  reader <- match.arg(getOption("lantern.vcf_reader", "auto"),
                      c("auto", "bcftools", "seqarray"))
  if (reader == "auto")
    reader <- if (nzchar(Sys.which("bcftools"))) "bcftools" else "seqarray"
  reader
}

#' Read the sample IDs of a VCF
#'
#' @keywords internal
.vcf_sample_ids <- function(vcf_path, reader) {
  if (reader == "seqarray") return(SeqArray::seqVCF_SampID(vcf_path))
  res <- system2("bcftools", c("query", "-l", shQuote(vcf_path)),
                 stdout = TRUE, stderr = FALSE)
  if (length(res) == 0) stop("bcftools query -l returned nothing")
  res
}

#' Read variant sites and GT strings from a VCF
#'
#' Both backends return the same thing: \code{chrom}, \code{pos},
#' \code{ref}, \code{alt} vectors and \code{gt}, a character matrix
#' (variants x all VCF samples) of GT strings exactly as
#' \code{bcftools query -f '[\%GT]'} prints them (e.g. \code{"0|1"},
#' \code{"./."}), so everything downstream is backend-agnostic.
#' \code{chrom} restricts the read where the backend can do so cheaply;
#' callers still filter by chromosome afterward.
#'
#' @keywords internal
.vcf_read_gt <- function(vcf_path, chrom, reader, verbose = TRUE) {
  if (reader == "seqarray")
    return(.vcf_read_gt_seqarray(vcf_path, chrom, verbose = verbose))

  target_arg <- character(0)
  if (!is.null(chrom)) {
    resolved <- .resolve_vcf_contig(vcf_path, chrom)
    if (!is.null(resolved)) {
      target_arg <- c("-t", shQuote(resolved))
      if (verbose) message("  Restricting bcftools query to contig '", resolved,
                           "' (other chromosomes are never read into R)")
    } else if (verbose) {
      message("  Could not match chromosome '", chrom, "' against the VCF's ",
              "##contig header; querying all variants and filtering ",
              "afterward instead (slower, more memory)")
    }
  }
  args_gt <- c("query", target_arg,
               "-f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"),
               shQuote(vcf_path))
  gt_output <- tryCatch(
    system2("bcftools", args_gt, stdout = TRUE, stderr = FALSE),
    error = function(e) stop("bcftools query failed: ", e$message))
  if (length(gt_output) == 0) stop("bcftools query returned no variants")
  if (verbose) message("  ", length(gt_output), " variant lines read")

  tmp_gt <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_gt), add = TRUE)
  writeLines(gt_output, tmp_gt)
  gt_dt <- data.table::fread(tmp_gt, header = FALSE, sep = "\t")

  list(chrom = gt_dt[[1]],
       pos   = as.integer(gt_dt[[2]]),
       ref   = as.character(gt_dt[[3]]),
       alt   = as.character(gt_dt[[4]]),
       gt    = as.matrix(gt_dt[, 5:ncol(gt_dt), drop = FALSE]))
}

#' @rdname dot-vcf_read_gt
#' @keywords internal
.vcf_read_gt_seqarray <- function(vcf_path, chrom, verbose = TRUE) {
  if (grepl("\\.bcf$", vcf_path, ignore.case = TRUE))
    stop("Reading BCF input requires bcftools on PATH; ",
         "convert it to VCF (.vcf / .vcf.gz) or install bcftools")
  if (verbose) message("  bcftools not found; reading the VCF via SeqArray")

  gds_tmp <- tempfile(fileext = ".gds")
  on.exit(unlink(gds_tmp), add = TRUE)
  # GT only; keep chromosome names verbatim (SeqArray strips "chr" by
  # default, bcftools does not) so variant_info matches the bcftools path.
  SeqArray::seqVCF2GDS(vcf_path, gds_tmp,
                       info.import = character(0), fmt.import = character(0),
                       ignore.chr.prefix = "", storage.option = "LZ4_RA.fast",
                       optimize = FALSE, digest = FALSE, verbose = FALSE)
  gds <- SeqArray::seqOpen(gds_tmp)
  on.exit(SeqArray::seqClose(gds), add = TRUE, after = FALSE)

  if (!is.null(chrom)) {
    keep <- sub("^chr", "", SeqArray::seqGetData(gds, "chromosome")) ==
      sub("^chr", "", chrom)
    if (!any(keep)) stop("No variants on chromosome ", chrom)
    SeqArray::seqSetFilter(gds, variant.sel = keep, verbose = FALSE)
  }

  sample_ids <- SeqArray::seqGetData(gds, "sample.id")
  n_samples  <- length(sample_ids)
  vcf_chrom  <- SeqArray::seqGetData(gds, "chromosome")
  n_variants <- length(vcf_chrom)
  if (n_variants == 0) stop("VCF contains no variants")
  if (verbose) message("  ", n_variants, " variants read")

  alt <- SeqArray::seqGetData(gds, "$alt")
  alt[alt == ""] <- "."   # bcftools prints a missing ALT as "."

  # genotype: ploidy x sample x variant allele indices (NA = missing);
  # phase: sample x variant, 1 = phased ("|"), 0 = unphased ("/").
  geno  <- SeqArray::seqGetData(gds, "genotype")
  geno  <- array(geno, dim = c(length(geno) / (n_samples * n_variants),
                               n_samples, n_variants))
  allele_str <- function(a) {
    s <- as.character(a)
    s[is.na(a)] <- "."
    s
  }
  if (dim(geno)[1] == 1) {
    gt <- allele_str(geno[1, , ])
  } else {
    phase <- SeqArray::seqGetData(gds, "phase")
    gt <- paste0(allele_str(geno[1, , ]),
                 ifelse(as.vector(phase) == 1, "|", "/"),
                 allele_str(geno[2, , ]))
  }
  gt <- t(matrix(gt, nrow = n_samples, ncol = n_variants))
  colnames(gt) <- sample_ids

  list(chrom = vcf_chrom,
       pos   = SeqArray::seqGetData(gds, "position"),
       ref   = SeqArray::seqGetData(gds, "$ref"),
       alt   = alt,
       gt    = gt)
}

#' Parse an MSP file and a phased VCF into haplotype/ancestry matrices
#'
#' Internal helper shared by \code{\link{ancestry_split}}: parses the RFMix
#' MSP file, queries VCF sample IDs, intersects samples, parses phased GT
#' into haplotype matrices, maps each variant to its ancestry tract, and
#' broadcasts tract-level ancestry to variant-level ancestry.
#'
#' @keywords internal
.parse_vcf_msp_common <- function(vcf_path, msp_path, chrom = NULL, verbose = TRUE) {

  # ---- Step 1: Parse MSP ----
  if (verbose) message("Step 1: Parsing RFMix MSP file...")
  msp_data     <- .parse_msp(msp_path, verbose = verbose)
  msp_samples  <- msp_data$sample_ids
  tract_df     <- msp_data$tract_df
  anc_hap0_mat <- msp_data$anc_hap0
  anc_hap1_mat <- msp_data$anc_hap1
  pop_codes    <- msp_data$pop_codes

  # Restrict tracts to the requested chromosome up front. Besides avoiding
  # wasted work on other chromosomes' tracts, this fixes a correctness gap:
  # without it, a multi-chromosome MSP file's tract_df is sorted by spos
  # alone (chromosome-agnostic), so findInterval() below could match a
  # variant to the wrong chromosome's tract purely by position overlap.
  if (!is.null(chrom)) {
    chrom_clean <- sub("^chr", "", chrom)
    tract_keep  <- sub("^chr", "", tract_df$chrom) == chrom_clean
    if (!any(tract_keep))
      stop("No ancestry tracts on chromosome ", chrom, " in MSP file")
    tract_df     <- tract_df[tract_keep, , drop = FALSE]
    anc_hap0_mat <- anc_hap0_mat[, tract_keep, drop = FALSE]
    anc_hap1_mat <- anc_hap1_mat[, tract_keep, drop = FALSE]
    rownames(tract_df) <- NULL
  }
  n_tracts <- nrow(tract_df)
  if (verbose) message("  MSP samples: ", length(msp_samples),
                       "  Tracts: ", n_tracts,
                       if (!is.null(chrom)) c(" (chr ", chrom, " only)"))

  # ---- Step 2: VCF sample IDs ----
  if (verbose) message("\nStep 2: Querying VCF sample IDs...")
  reader <- .vcf_reader()
  vcf_samples <- tryCatch(
    .vcf_sample_ids(vcf_path, reader),
    error = function(e) stop("Could not read VCF sample IDs: ", e$message))
  if (verbose) message("  VCF samples: ", length(vcf_samples))

  # ---- Step 3: Intersect samples ----
  if (verbose) message("\nStep 3: Intersecting samples...")
  common_samples <- intersect(vcf_samples, msp_samples)
  if (length(common_samples) == 0)
    stop("No common samples between VCF and MSP")
  dropped_vcf <- setdiff(vcf_samples, common_samples)
  dropped_msp <- setdiff(msp_samples, common_samples)
  if (verbose) message("  Common: ", length(common_samples))

  anc_hap0_common <- anc_hap0_mat[common_samples, , drop = FALSE]
  anc_hap1_common <- anc_hap1_mat[common_samples, , drop = FALSE]

  # ---- Step 4: Parse VCF GT ----
  if (verbose) message("\nStep 4: Parsing VCF genotypes...")
  vcf <- .vcf_read_gt(vcf_path, chrom, reader, verbose = verbose)
  vcf_chrom <- vcf$chrom
  vcf_pos   <- vcf$pos
  vcf_ref   <- vcf$ref
  vcf_alt   <- vcf$alt

  is_biallelic <- !grepl(",", vcf_alt, fixed = TRUE)
  n_multi <- sum(!is_biallelic)
  if (n_multi > 0) {
    if (verbose) message("  Skipping ", n_multi, " multiallelic variants")
    vcf_chrom <- vcf_chrom[is_biallelic]; vcf_pos <- vcf_pos[is_biallelic]
    vcf_ref   <- vcf_ref[is_biallelic];   vcf_alt <- vcf_alt[is_biallelic]
  }
  n_variants <- length(vcf_pos)

  sample_idx <- match(common_samples, vcf_samples)
  gt_mat     <- vcf$gt[is_biallelic, sample_idx, drop = FALSE]

  gt_parsed  <- .parse_phased_gt_matrix(gt_mat, n_variants,
                                         length(common_samples),
                                         common_samples, verbose = verbose)
  gt_hap0_mat <- gt_parsed$hap0
  gt_hap1_mat <- gt_parsed$hap1

  n_missing <- sum(gt_parsed$missing)
  if (n_missing > 0 && verbose)
    message("  Missing GT: ", n_missing, " (treated as 0/0)")

  # ---- Step 5: Map variants to tracts ----
  if (verbose) message("\nStep 5: Mapping variants to ancestry tracts...")
  vcf_chrom_clean   <- sub("^chr", "", vcf_chrom)
  tract_chrom_clean <- sub("^chr", "", tract_df$chrom)

  if (!is.null(chrom)) {
    chrom_clean <- sub("^chr", "", chrom)
    keep <- vcf_chrom_clean == chrom_clean
    if (!any(keep)) stop("No variants on chromosome ", chrom)
    vcf_chrom <- vcf_chrom[keep]; vcf_pos <- vcf_pos[keep]
    vcf_ref   <- vcf_ref[keep];   vcf_alt <- vcf_alt[keep]
    gt_hap0_mat <- gt_hap0_mat[keep, , drop = FALSE]
    gt_hap1_mat <- gt_hap1_mat[keep, , drop = FALSE]
    vcf_chrom_clean <- vcf_chrom_clean[keep]
    n_variants <- sum(keep)
    if (verbose) message("  Filtered to chrom ", chrom, ": ", n_variants, " variants")
  }

  tract_idx   <- findInterval(vcf_pos, tract_df$spos)
  valid_tract <- tract_idx > 0
  if (any(valid_tract))
    valid_tract[valid_tract] <-
      vcf_pos[valid_tract] <= tract_df$epos[tract_idx[valid_tract]]
  # tract_df is already restricted to `chrom` (Step 1) whenever chrom is
  # supplied, so this check is cheap and always correct either way -- no
  # special-cased shortcut needed.
  valid_chrom <- vapply(seq_along(tract_idx), function(i) {
    tract_idx[i] > 0 &&
      vcf_chrom_clean[i] == tract_chrom_clean[tract_idx[i]]
  }, logical(1))
  keep_var <- valid_tract & valid_chrom
  n_no_tract <- sum(!keep_var)
  if (n_no_tract > 0 && verbose)
    message("  Variants with no tract: ", n_no_tract, " (dropped)")

  if (n_no_tract > 0) {
    vcf_chrom   <- vcf_chrom[keep_var]; vcf_pos <- vcf_pos[keep_var]
    vcf_ref     <- vcf_ref[keep_var];   vcf_alt <- vcf_alt[keep_var]
    gt_hap0_mat <- gt_hap0_mat[keep_var, , drop = FALSE]
    gt_hap1_mat <- gt_hap1_mat[keep_var, , drop = FALSE]
    tract_idx   <- tract_idx[keep_var]
    n_variants  <- sum(keep_var)
  }
  if (n_variants == 0) stop("No variants with valid ancestry tracts")

  # ---- Step 6: Broadcast tract ancestry to variants ----
  if (verbose) message("\nStep 6: Broadcasting ancestry to variants...")
  anc_hap0_var <- matrix(0L, nrow = n_variants, ncol = length(common_samples))
  anc_hap1_var <- matrix(0L, nrow = n_variants, ncol = length(common_samples))
  for (i in seq_len(n_variants)) {
    t <- tract_idx[i]
    anc_hap0_var[i, ] <- anc_hap0_common[, t]
    anc_hap1_var[i, ] <- anc_hap1_common[, t]
  }

  variant_info <- data.frame(chrom = vcf_chrom, pos = vcf_pos,
                              ref = vcf_ref, alt = vcf_alt,
                              stringsAsFactors = FALSE)

  list(
    gt_hap0_mat     = gt_hap0_mat,
    gt_hap1_mat     = gt_hap1_mat,
    anc_hap0_var    = anc_hap0_var,
    anc_hap1_var    = anc_hap1_var,
    anc_hap0_tract  = anc_hap0_common,
    anc_hap1_tract  = anc_hap1_common,
    variant_info    = variant_info,
    common_samples  = common_samples,
    pop_codes       = pop_codes,
    tract_df        = tract_df,
    overlap = list(
      n_vcf_samples           = length(vcf_samples),
      n_msp_samples           = length(msp_samples),
      n_common                = length(common_samples),
      n_multiallelic_filtered = n_multi,
      n_no_tract              = n_no_tract,
      dropped_samples_vcf     = dropped_vcf,
      dropped_samples_msp     = dropped_msp
    )
  )
}

# ============================================================================
# ancestry_split (Step 1 entry point)
# ============================================================================

#' Split a phased VCF by local ancestry (Step 1)
#'
#' Parse an RFMix MSP file and a phased VCF, intersect samples, map each
#' variant to its ancestry tract, and split genotypes into per-population
#' dosage matrices using either the proportional (\code{"dosage"}) or
#' deterministic (\code{"haplotype"}) split algorithm. This is the
#' package's Step 1 entry point: it only computes in-memory dosage
#' matrices, with no file I/O. Pass the result to
#' \code{\link{write_ancestry_gds}} (Step 2) to materialize GDS files, and
#' then to \code{\link{ancestry_smmat}} (Step 3) to run SMMAT + Cauchy
#' combination.
#'
#' @param vcf_path Path to phased VCF/BCF file (plain or gzipped). Read
#'   with \code{bcftools} if it is on \code{PATH} (fastest; required for
#'   BCF), otherwise with \pkg{SeqArray}. Set
#'   \code{options(lantern.vcf_reader = "bcftools")} or \code{"seqarray"}
#'   to force one.
#' @param msp_path Path to RFMix MSP file (plain text or gzipped TSV).
#' @param mode Split algorithm: \code{"dosage"} (proportional p1/p2 split
#'   via \code{\link{split_diploid_multi}}, unphased) or \code{"haplotype"}
#'   (deterministic per-haplotype split via
#'   \code{\link{split_haplotype_multi}}, phased). Either/or -- call
#'   \code{ancestry_split()} twice if both are needed.
#' @param chrom Chromosome to process (e.g., \code{"chr19"} or \code{"19"}).
#'   If \code{NULL}, all chromosomes present in the VCF are used. When
#'   supplied, \code{vcf_path} and \code{msp_path} may contain other
#'   chromosomes too: MSP tracts outside \code{chrom} are dropped up front,
#'   and \code{bcftools query} is restricted to \code{chrom} via \code{-t}
#'   when the VCF's declared or inferred contigs allow resolving the exact
#'   name to use (whichever of \code{chrom}, with, or without a leading
#'   "chr" matches) -- other chromosomes' genotypes are never read into R.
#'   If the contig can't be resolved, all variants are queried and filtered
#'   in R instead (slower, more memory, same result).
#' @param verbose Print step-by-step progress messages.
#'
#' @section Shrinkage toward 1/2 (dosage mode): A mixed-ancestry
#'   heterozygote carries one haplotype from each of its two ancestries, so
#'   which one holds the alt allele cannot be read off the genotype. Its
#'   allele is split between the two ancestries in proportion to the
#'   variant's unambiguous carriers (the raw ratio), shrunk toward 1/2 by
#'   \eqn{w}, the share of the variant's alt-allele evidence that comes from
#'   ambiguous heterozygotes:
#'   \deqn{p = (1 - w)\, p_{raw} + w / 2.}
#'   A variant carried only by ambiguous heterozygotes (\eqn{w = 1}) is
#'   split 1/2 : 1/2. With K populations every pair is shrunk toward 1/2
#'   in the same way. \code{mode = "haplotype"} is deterministic and
#'   applies no shrinkage.
#'
#' @return Invisibly, a list with one named numeric matrix (variants x
#'   samples) per population (named after the MSP population codes, e.g.
#'   \code{$AFR}, \code{$EUR}), plus:
#'   \item{variant_info}{data.frame with columns chrom, pos, ref, alt.}
#'   \item{sample_ids}{Character vector of common sample IDs.}
#'   \item{mode}{The split mode used.}
#'   \item{ancestry_counts}{Numeric matrix (variants x populations): for
#'     each variant, the number of samples with pure ancestry in that
#'     population. Feeds \code{\link{ancestry_smmat}}'s per-gene weight
#'     calculation.}
#'   \item{tract_info}{data.frame of ancestry tracts from the MSP file.}
#'   \item{overlap}{List of intersection/filtering statistics.}
#'
#' @seealso \code{\link{write_ancestry_gds}}, \code{\link{ancestry_smmat}}
#'
#' @examples
#' vcf <- system.file("extdata", "toy.vcf.gz", package = "lantern")
#' msp <- system.file("extdata", "toy.msp.tsv", package = "lantern")
#' split <- ancestry_split(vcf, msp, mode = "dosage", chrom = "chr19",
#'                         verbose = FALSE)
#' split$AFR[1:3, 1:5]   # African dosages (variants x samples)
#' split$EUR[1:3, 1:5]   # European dosages
#' head(split$variant_info)
#'
#' @export
ancestry_split <- function(vcf_path, msp_path, mode = c("dosage", "haplotype"),
                            chrom = NULL, verbose = TRUE) {
  mode <- match.arg(mode)

  if (verbose) message("=== LANTERN Ancestry Split (Step 1, mode = ", mode, ") ===\n")

  common <- .parse_vcf_msp_common(vcf_path, msp_path, chrom = chrom, verbose = verbose)

  gt_hap0_mat    <- common$gt_hap0_mat
  gt_hap1_mat    <- common$gt_hap1_mat
  anc_hap0_var   <- common$anc_hap0_var
  anc_hap1_var   <- common$anc_hap1_var
  variant_info   <- common$variant_info
  common_samples <- common$common_samples
  pop_codes      <- common$pop_codes
  n_variants     <- nrow(variant_info)
  K              <- length(pop_codes)
  if (K < 2) stop("ancestry_split() requires at least 2 populations in the MSP")
  pop_names      <- names(pop_codes)

  # ---- Auto-generate diploid ancestry codes for K populations ----
  # pure_codes_named: pop_name -> unique diploid code 1..K (pop index)
  # mixed_codes_df  : data.frame(code, pop1, pop2) for all K*(K-1)/2 pairs
  # lookup[p0_idx, p1_idx] -> diploid code  (vectorised ancestry conversion)
  hap_codes_v      <- as.integer(pop_codes)
  pure_codes_named <- setNames(seq_len(K), pop_names)
  pairs_mat        <- combn(K, 2)
  M                <- ncol(pairs_mat)
  mixed_code_v     <- K + seq_len(M)
  mixed_codes_df   <- data.frame(
    code = as.integer(mixed_code_v),
    pop1 = pop_names[pairs_mat[1L, ]],
    pop2 = pop_names[pairs_mat[2L, ]],
    stringsAsFactors = FALSE
  )
  lookup <- matrix(0L, nrow = K, ncol = K)
  diag(lookup) <- seq_len(K)
  for (m in seq_len(M)) {
    i <- pairs_mat[1L, m]; j <- pairs_mat[2L, m]
    lookup[i, j] <- mixed_code_v[m]
    lookup[j, i] <- mixed_code_v[m]
  }

  p0_mat <- matrix(match(as.vector(anc_hap0_var), hap_codes_v),
                    nrow = n_variants, ncol = length(common_samples))
  p1_mat <- matrix(match(as.vector(anc_hap1_var), hap_codes_v),
                    nrow = n_variants, ncol = length(common_samples))
  anc_diploid <- matrix(lookup[cbind(as.vector(p0_mat), as.vector(p1_mat))],
                         nrow = n_variants, ncol = length(common_samples))

  # ---- Shrinkage target (dosage mode; haplotype splitting is deterministic
  # and has no ambiguous ratio to shrink): uniform 1/K for every population,
  # i.e. each ambiguous pair's split is shrunk toward 1/2 by the share w of
  # its evidence that is ambiguous. One target row shared by all variants,
  # so no genomic position / arm lookup is needed. ----
  shrink_target <- NULL
  arm_id        <- NULL
  if (mode == "dosage") {
    shrink_target <- .uniform_shrink_target(pop_names)
    arm_id        <- integer(n_variants)
  }

  # ---- Split by ancestry ----
  if (verbose) message("\nStep 7: Splitting genotypes by ancestry (", mode, ")...")
  if (mode == "haplotype") {
    res <- split_haplotype_multi(gt_hap0_mat, gt_hap1_mat,
                                  anc_hap0_var, anc_hap1_var,
                                  pop_codes = pop_codes)
  } else {
    gt_diploid <- gt_hap0_mat + gt_hap1_mat
    res <- split_diploid_multi(gt_diploid, anc_diploid,
                                pure_codes_named, mixed_codes_df,
                                gla = shrink_target, arm_id = arm_id)
  }

  # ---- Per-variant pure-ancestry sample counts (feeds Step 3 gene weights) ----
  # matrix() wrap avoids vapply's auto-simplification collapsing to a plain
  # vector when n_variants == 1
  ancestry_counts <- matrix(
    vapply(pure_codes_named, function(code) count_ancestry_codes(anc_diploid, code),
           numeric(n_variants)),
    nrow = n_variants, ncol = K, dimnames = list(NULL, pop_names))

  # ---- Filter monomorphic variants ----
  if (verbose) message("\nStep 8: Filtering monomorphic variants...")
  total   <- Reduce("+", res)
  has_alt <- rowSums(total) > 0
  n_mono  <- sum(!has_alt)
  if (n_mono > 0) {
    if (verbose) message("  Removed ", n_mono, " monomorphic variants")
    res             <- lapply(res, function(m) m[has_alt, , drop = FALSE])
    variant_info    <- variant_info[has_alt, , drop = FALSE]
    ancestry_counts <- ancestry_counts[has_alt, , drop = FALSE]
  }
  n_final <- nrow(variant_info)
  if (n_final == 0) stop("All variants are monomorphic after filtering")
  rownames(variant_info) <- NULL

  if (verbose) message("\n=== Step 1 Complete: ", n_final, " variants, ",
                       length(common_samples), " samples, ", K, " populations ===\n")

  overlap <- c(common$overlap,
               n_variants_kept        = n_final,
               n_monomorphic_filtered = n_mono)

  invisible(c(res, list(
    variant_info    = variant_info,
    sample_ids      = common_samples,
    mode            = mode,
    ancestry_counts = ancestry_counts,
    tract_info      = common$tract_df,
    overlap         = overlap
  )))
}

# ============================================================================
# ancestry_split_dosage (backward-compatible Step 1 wrapper)
# ============================================================================

#' Run ancestry splitting pipeline
#'
#' Complete pipeline: split genotype matrix by ancestry, return dosage matrices.
#' Handles sample and variant mismatches by using only overlapping data.
#' Memory efficient: single-pass subsetting using index vectors.
#'
#' @param gt_matrix Integer matrix of genotypes (rows=variants, cols=samples).
#'   Values: 0=homozygous ref, 1=heterozygous, 2=homozygous alt.
#'   rownames: variant IDs (e.g., \code{"chr:pos"}). colnames: sample IDs.
#'   May be \code{NULL} when \code{vcf_path} and \code{msp_path} are supplied.
#' @param pt_matrix Integer matrix of parent-of-origin ancestry codes
#'   (rows=samples, cols=regions/variants).
#'   Values: 1=EUR/EUR, 2=AFR/EUR (mixed), 3=AFR/AFR.
#'   rownames: sample IDs (must match colnames of \code{gt_matrix}).
#'   colnames: region IDs (e.g., \code{"chr:start-end"}).
#'   May be \code{NULL} when \code{vcf_path} and \code{msp_path} are supplied.
#' @param vcf_path Path to a phased VCF or BCF file.  When supplied together
#'   with \code{msp_path}, the function delegates to
#'   \code{\link{ancestry_split}} (Step 1, \code{mode = "dosage"}), bypassing
#'   the \code{gt_matrix}/\code{pt_matrix} arguments.
#' @param msp_path Path to an RFMix MSP file (plain text or gzipped TSV).
#'   Used together with \code{vcf_path} as a shortcut input.
#' @param chrom Chromosome name to restrict to (e.g. \code{"chr19"}).
#'   Only used when \code{vcf_path}/\code{msp_path} are supplied.
#' @param verbose Print progress messages (default TRUE)
#'
#' @section Shrinkage: Both input paths split ambiguous mixed-ancestry
#'   heterozygotes with the same shrinkage toward 1/2 as
#'   \code{\link{ancestry_split}} (see its "Shrinkage toward 1/2" section).
#'   The target is the same for every variant, so the direct
#'   \code{gt_matrix}/\code{pt_matrix} path needs no position information.
#'
#' @return List with elements:
#'   \item{african}{African ancestry-specific dosage matrix}
#'   \item{european}{European ancestry-specific dosage matrix}
#'   \item{counts}{List of ancestry counts per region/sample}
#'   \item{overlap}{Overlap information}
#'
#' @examples
#' # Example with mismatched dimensions
#' gt <- matrix(c(2, 1, 0, 1, 2, 1), nrow = 2, ncol = 3,
#'              dimnames = list(c("chr1:100", "chr1:200"),
#'                             c("sample_A", "sample_B", "sample_C")))
#'
#' pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
#'              dimnames = list(c("sample_A", "sample_B", "sample_D"),
#'                             c("chr1:50-150", "chr1:150-250")))
#'
#' result <- ancestry_split_dosage(gt, pt)
#' # Only sample_A and sample_B are kept (common to both)
#' # Only chr1:100 and chr1:200 that overlap regions are kept
#'
#' # Shortcut: build matrices directly from VCF + MSP
#' vcf <- system.file("extdata", "toy.vcf.gz", package = "lantern")
#' msp <- system.file("extdata", "toy.msp.tsv", package = "lantern")
#' result <- ancestry_split_dosage(vcf_path = vcf, msp_path = msp,
#'                                 chrom = "chr19", verbose = FALSE)
#'
#' @export
ancestry_split_dosage <- function(gt_matrix = NULL, pt_matrix = NULL,
                                   vcf_path = NULL, msp_path = NULL,
                                   chrom = NULL, verbose = TRUE) {

  # ---- VCF/MSP shortcut path: delegate to ancestry_split() (Step 1) ----
  if (!is.null(vcf_path) && !is.null(msp_path)) {
    split_result <- ancestry_split(vcf_path, msp_path, mode = "dosage",
                                    chrom = chrom, verbose = verbose)
    pop_names <- setdiff(names(split_result),
                          c("variant_info", "sample_ids", "mode",
                            "ancestry_counts", "tract_info", "overlap"))
    if (length(pop_names) != 2)
      stop("ancestry_split_dosage()'s vcf_path/msp_path shortcut only ",
           "supports 2-population (AFR/EUR) MSP files; use ancestry_split() ",
           "directly for K > 2 populations.")
    afr_name <- pop_names[toupper(pop_names) %in% c("AFR", "AFRICAN")]
    eur_name <- pop_names[toupper(pop_names) %in% c("EUR", "EUROPEAN")]
    if (length(afr_name) != 1 || length(eur_name) != 1) {
      afr_name <- pop_names[1]; eur_name <- pop_names[2]
    }
    n_samp     <- length(split_result$sample_ids)
    afr_counts <- split_result$ancestry_counts[, afr_name]
    eur_counts <- split_result$ancestry_counts[, eur_name]
    return(invisible(list(
      african  = split_result[[afr_name]],
      european = split_result[[eur_name]],
      counts   = list(african  = afr_counts,
                      european = eur_counts,
                      mixed    = n_samp - afr_counts - eur_counts),
      overlap  = split_result$overlap
    )))
  }
  if (is.null(gt_matrix) || is.null(pt_matrix))
    stop("Provide either (gt_matrix, pt_matrix) or (vcf_path, msp_path)")

  if (verbose) message("=== LANTERN Ancestry Pipeline ===\n")

  # Record original dimensions
  n_gt_samples_orig <- ncol(gt_matrix)
  n_gt_variants_orig <- nrow(gt_matrix)
  n_pt_samples_orig <- nrow(pt_matrix)
  n_pt_regions_orig <- ncol(pt_matrix)

  # ========================================================================
  # Step 1: Find sample overlap (memory efficient - use indices, not copies)
  # ========================================================================
  if (verbose) message("Step 1: Finding sample overlap...")

  gt_sample_names <- colnames(gt_matrix)
  pt_sample_names <- rownames(pt_matrix)

  # Handle unnamed matrices
  if (is.null(gt_sample_names) || is.null(pt_sample_names)) {
    n_gt <- ncol(gt_matrix)
    n_pt <- nrow(pt_matrix)
    min_n <- min(n_gt, n_pt)
    if (n_gt != n_pt) {
      warning("Matrices have different dimensions (GT: ", n_gt, " cols, PT: ", n_pt,
              " rows). Using first ", min_n, " samples.")
    }
    gt_sample_names <- paste0("sample_", seq_len(min_n))
    pt_sample_names <- paste0("sample_", seq_len(min_n))
    colnames(gt_matrix) <- gt_sample_names
    rownames(pt_matrix) <- pt_sample_names
  }

  # Find common samples by name
  common_samples <- intersect(gt_sample_names, pt_sample_names)

  if (length(common_samples) == 0) {
    stop("No common samples found between GT and PT matrices")
  }

  # Calculate indices for subsetting (avoid copying until needed)
  gt_sample_idx <- match(common_samples, gt_sample_names)
  pt_sample_idx <- match(common_samples, pt_sample_names)

  n_common_samples <- length(common_samples)
  dropped_samples <- c(
    setdiff(gt_sample_names, common_samples),
    setdiff(pt_sample_names, common_samples)
  )

  if (verbose) {
    message("  Samples: ", n_gt_samples_orig, " in GT, ", n_pt_samples_orig, " in PT")
    message("  Common samples: ", n_common_samples)
    if (length(dropped_samples) > 0) {
      message("  Dropped: ", paste(head(dropped_samples, 3), collapse = ", "),
              if (length(dropped_samples) > 3) "..." else "")
    }
  }

  # ========================================================================
  # Step 2: Find variant/region overlap
  # ========================================================================
  if (verbose) message("\nStep 2: Finding variant/region overlap...")

  gt_var_names <- rownames(gt_matrix)
  pt_region_names <- colnames(pt_matrix)

  # Handle unnamed matrices
  if (is.null(gt_var_names) || is.null(pt_region_names)) {
    n_gt <- nrow(gt_matrix)
    n_pt <- ncol(pt_matrix)
    min_n <- min(n_gt, n_pt)
    if (n_gt != n_pt) {
      warning("Matrices have different dimensions (GT: ", n_gt, " rows, PT: ",
              n_pt, " cols). Using first ", min_n, " variants.")
    }
    gt_var_names <- paste0("var_", seq_len(min_n))
    pt_region_names <- paste0("var_", seq_len(min_n))
    rownames(gt_matrix) <- gt_var_names
    colnames(pt_matrix) <- pt_region_names
  }

  # Try exact matching first (fast)
  common_vars <- intersect(gt_var_names, pt_region_names)

  # If no exact match, try coordinate-based matching
  if (length(common_vars) == 0) {
    if (verbose) message("  No exact matches. Trying coordinate-based matching...")

    # Pre-parse GT coordinates: chr:pos -> list(chr, pos)
    parse_gt_coord <- function(name) {
      parts <- strsplit(name, ":")[[1]]
      if (length(parts) >= 2) {
        pos <- suppressWarnings(as.numeric(parts[2]))
        if (is.na(pos)) pos <- NULL
        list(chr = parts[1], pos = pos)
      } else {
        NULL
      }
    }

    # Pre-parse PT coordinates: chr:start-end -> list(chr, start, end)
    parse_pt_coord <- function(name) {
      chr_part <- strsplit(name, ":")[[1]]
      if (length(chr_part) >= 2) {
        coords <- strsplit(chr_part[2], "-")[[1]]
        if (length(coords) == 2) {
          start <- suppressWarnings(as.numeric(coords[1]))
          end <- suppressWarnings(as.numeric(coords[2]))
          if (!is.na(start) && !is.na(end)) {
            return(list(chr = chr_part[1], start = start, end = end))
          }
        }
      }
      NULL
    }

    # Parse all coordinates upfront (vectorized approach)
    gt_coords <- lapply(gt_var_names, parse_gt_coord)
    pt_coords <- lapply(pt_region_names, parse_pt_coord)

    # Find valid matches using matrix operations where possible
    matched_gt_idx <- integer()
    matched_pt_idx <- integer()

    # Build match matrix for O(n*m) scan
    for (i in seq_along(gt_coords)) {
      if (is.null(gt_coords[[i]]$pos)) next
      for (j in seq_along(pt_coords)) {
        pt <- pt_coords[[j]]
        if (is.null(pt$start)) next
        if (gt_coords[[i]]$chr == pt$chr &&
            gt_coords[[i]]$pos >= pt$start &&
            gt_coords[[i]]$pos <= pt$end) {
          matched_gt_idx <- c(matched_gt_idx, i)
          matched_pt_idx <- c(matched_pt_idx, j)
          break  # Take first matching region
        }
      }
    }

    # Clean up intermediate objects
    rm(gt_coords, pt_coords)
    gc()

    if (length(matched_gt_idx) > 0) {
      common_vars <- unique(c(gt_var_names[matched_gt_idx], pt_region_names[matched_pt_idx]))
      if (verbose) {
        message("  Found ", length(matched_gt_idx), " variants in ",
                length(unique(matched_pt_idx)), " regions")
      }
    } else {
      # Fallback to positional match by order
      warning("No variant/region name or coordinate overlap. Matching by order.")
      min_n <- min(length(gt_var_names), length(pt_region_names))
      rownames(gt_matrix) <- paste0("var_", seq_len(nrow(gt_matrix)))
      colnames(pt_matrix) <- paste0("region_", seq_len(ncol(pt_matrix)))
      common_vars <- paste0("shared_", seq_len(min_n))
    }
  }

  # Calculate indices for variant/region subsetting
  gt_var_idx <- match(intersect(gt_var_names, common_vars), gt_var_names)
  pt_region_idx <- match(intersect(pt_region_names, common_vars), pt_region_names)

  n_common_vars <- length(common_vars)
  dropped_variants <- c(
    setdiff(gt_var_names, common_vars),
    setdiff(pt_region_names, common_vars)
  )

  if (verbose) {
    message("  Variants in GT: ", n_gt_variants_orig)
    message("  Regions in PT: ", n_pt_regions_orig)
    message("  Common variants/regions: ", n_common_vars)
    if (length(dropped_variants) > 0) {
      message("  Dropped: ", paste(head(dropped_variants, 3), collapse = ", "),
              if (length(dropped_variants) > 3) "..." else "")
    }
  }

  # ========================================================================
  # Step 3: Single-pass subsetting (memory efficient)
  # ========================================================================
  if (verbose) message("\nStep 3: Subsetting matrices...")

  # Subset both matrices in single operation using calculated indices
  # GT: rows = variants we want, cols = samples we want
  # PT: rows = samples we want, cols = regions we want
  gt_subset <- gt_matrix[gt_var_idx, gt_sample_idx, drop = FALSE]
  pt_subset <- pt_matrix[pt_sample_idx, pt_region_idx, drop = FALSE]

  # Clean up intermediates immediately
  rm(gt_matrix, pt_matrix)
  gc()

  # ========================================================================
  # Step 4: Validate dimensions
  # ========================================================================
  if (ncol(gt_subset) != nrow(pt_subset)) {
    stop("After overlap filtering: GT cols (", ncol(gt_subset),
         ") != PT rows (", nrow(pt_subset), ")")
  }
  if (nrow(gt_subset) != ncol(pt_subset)) {
    stop("After overlap filtering: GT rows (", nrow(gt_subset),
         ") != PT cols (", ncol(pt_subset), "). ",
         "Each GT variant should match to exactly one PT region.")
  }

  # ========================================================================
  # Step 4: Filter monomorphic variants (no alt alleles)
  # ========================================================================
  if (verbose) message("\nStep 3b: Filtering monomorphic variants...")

  variant_has_alt <- rowSums(gt_subset > 0) > 0
  n_monomorphic <- sum(!variant_has_alt)

  if (n_monomorphic > 0) {
    if (verbose) message("  Removed ", n_monomorphic, " monomorphic variants (no alt alleles)")

    dropped_variants <- c(dropped_variants, rownames(gt_subset)[!variant_has_alt])

    gt_subset <- gt_subset[variant_has_alt, , drop = FALSE]
    pt_subset <- pt_subset[, variant_has_alt, drop = FALSE]
  }

  if (nrow(gt_subset) == 0) {
    stop("All variants are monomorphic (no alt alleles)")
  }

  # ========================================================================
  # Step 5: Split genotypes by ancestry (C backend)
  # ========================================================================
  if (verbose) message("\nStep 4: Splitting genotypes by ancestry...")

  # Same shrinkage toward 1/2 as ancestry_split(); the target is uniform, so
  # it needs no variant positions and stays aligned after the filtering above.
  result <- split_diploid(gt_subset, pt_subset,
                          gla = .uniform_shrink_target(c("AFR", "EUR")),
                          arm_id = integer(nrow(gt_subset)))

  # Capture post-filter variant count BEFORE removing gt_subset
  n_variants_kept_final <- nrow(gt_subset)

  # Clean up gt_subset after C call (keep pt_subset for counts - it's small)
  rm(gt_subset)
  gc()

  if (verbose) {
    message("  -> African dosage matrix: ", nrow(result$african), " x ", ncol(result$african))
    message("  -> European dosage matrix: ", nrow(result$european), " x ", ncol(result$european))
  }

  # ========================================================================
  # Step 6: Count ancestries per region
  # ========================================================================
  if (verbose) message("\nStep 5: Counting ancestries per region...")

  counts <- list(
    african = count_ancestry_codes(pt_subset, 3),
    european = count_ancestry_codes(pt_subset, 1),
    mixed = count_ancestry_codes(pt_subset, 2)
  )

  if (verbose) {
    message("  -> African (3): median = ", median(counts$african))
    message("  -> European (1): median = ", median(counts$european))
    message("  -> Mixed (2): median = ", median(counts$mixed))
  }

  if (verbose) message("\n=== Pipeline Complete ===\n")

  invisible(list(
    african = result$african,
    european = result$european,
    counts = counts,
    overlap = list(
      n_samples_total = n_gt_samples_orig + n_pt_samples_orig - n_common_samples,
      n_samples_kept = n_common_samples,
      n_variants_total = n_gt_variants_orig,
      n_variants_kept = n_variants_kept_final,
      n_monomorphic_filtered = n_monomorphic,
      dropped_samples = unique(dropped_samples),
      dropped_variants = unique(dropped_variants)
    )
  ))
}

# ============================================================================
# ancestry_split_phased (backward-compatible Step 1 + Step 2 wrapper)
# ============================================================================

#' Run phased ancestry splitting pipeline
#'
#' Backward-compatible 2-population (AFR/EUR) wrapper around
#' \code{\link{ancestry_split}} (Step 1, \code{mode = "haplotype"}) and,
#' when \code{write_vcf = TRUE}, \code{\link{write_ancestry_gds}} (Step 2).
#' For K > 2 populations or to keep the split and GDS-writing steps
#' separate, call \code{ancestry_split()} and \code{write_ancestry_gds()}
#' directly.
#'
#' @param vcf_path Path to phased VCF/BCF file (plain or gzipped). Read
#'   with \code{bcftools} if it is on \code{PATH} (fastest; required for
#'   BCF), otherwise with \pkg{SeqArray}. Set
#'   \code{options(lantern.vcf_reader = "bcftools")} or \code{"seqarray"}
#'   to force one.
#' @param msp_path Path to RFMix MSP file (plain text or gzipped TSV).
#' @param out_path Output directory for GDS and cache files.
#' @param chrom Chromosome to process (e.g., \code{"chr19"} or \code{"19"}).
#'   If \code{NULL}, all chromosomes present in the VCF are used.
#' @param write_vcf Logical; if \code{TRUE}, materialize ancestry-specific
#'   GDS files via \code{\link{write_ancestry_gds}} (Step 2).
#' @param verbose Print step-by-step progress messages.
#'
#' @return Invisibly, a list with elements:
#'   \item{african}{Numeric matrix (variants x samples) of African dosages.}
#'   \item{european}{Numeric matrix (variants x samples) of European dosages.}
#'   \item{variant_info}{data.frame with columns chrom, pos, ref, alt.}
#'   \item{sample_ids}{Character vector of common sample IDs.}
#'   \item{tract_info}{data.frame of ancestry tracts from the MSP file.}
#'   \item{overlap}{List of intersection statistics.}
#'   \item{vcf_paths}{(when \code{write_vcf = TRUE}) Named list with
#'     \code{african_gds}/\code{european_gds} paths.}
#'
#' @examples
#' vcf <- system.file("extdata", "toy.vcf.gz", package = "lantern")
#' msp <- system.file("extdata", "toy.msp.tsv", package = "lantern")
#' out_dir <- tempfile()
#' result <- ancestry_split_phased(vcf, msp, out_path = out_dir,
#'                                 chrom = "chr19", write_vcf = FALSE,
#'                                 verbose = FALSE)
#' head(result$variant_info)
#' unlink(out_dir, recursive = TRUE)
#'
#' @export
ancestry_split_phased <- function(vcf_path, msp_path, out_path,
                                chrom = NULL, write_vcf = TRUE,
                                verbose = TRUE) {

  if (verbose) message("=== LANTERN Phased Ancestry Pipeline ===\n")

  split_result <- ancestry_split(vcf_path, msp_path, mode = "haplotype",
                                  chrom = chrom, verbose = verbose)

  pop_names <- setdiff(names(split_result),
                        c("variant_info", "sample_ids", "mode",
                          "ancestry_counts", "tract_info", "overlap"))
  if (length(pop_names) != 2)
    stop("ancestry_split_phased() only supports 2-population (AFR/EUR) MSP ",
         "files; use ancestry_split() directly for K > 2 populations.")

  afr_name <- pop_names[toupper(pop_names) %in% c("AFR", "AFRICAN")]
  eur_name <- pop_names[toupper(pop_names) %in% c("EUR", "EUROPEAN")]
  if (length(afr_name) != 1 || length(eur_name) != 1) {
    afr_name <- pop_names[1]; eur_name <- pop_names[2]
  }
  african_mat  <- split_result[[afr_name]]
  european_mat <- split_result[[eur_name]]

  vcf_paths <- NULL
  if (write_vcf) {
    if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
    gds_paths <- write_ancestry_gds(split_result, out_path, verbose = verbose)
    vcf_paths <- list(african_gds  = unname(gds_paths[[afr_name]]),
                      european_gds = unname(gds_paths[[eur_name]]))
  }

  if (verbose) message("\n=== Pipeline Complete ===\n")

  invisible(list(african      = african_mat,
                 european     = european_mat,
                 variant_info = split_result$variant_info,
                 sample_ids   = split_result$sample_ids,
                 tract_info   = split_result$tract_info,
                 overlap      = split_result$overlap,
                 vcf_paths    = vcf_paths))
}
