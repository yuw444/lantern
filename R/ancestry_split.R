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

# Checks for the optional shrinkage target of split_diploid()/
# split_diploid_multi(): both or neither, one target row per arm_id value.
.check_shrink_target <- function(gla, arm_id, n_variants, n_pops) {
  if (is.null(gla) != is.null(arm_id))
    stop("`gla` and `arm_id` must be supplied together", call. = FALSE)
  if (is.null(gla)) return(invisible(NULL))
  .check_matrix(gla, "gla")
  if (ncol(gla) != n_pops)
    stop("`gla` must have one column per population (", n_pops, ")",
         call. = FALSE)
  if (any(gla < 0, na.rm = TRUE) || anyNA(gla))
    stop("`gla` must contain non-negative proportions", call. = FALSE)
  if (length(arm_id) != n_variants || anyNA(arm_id) ||
      any(arm_id < 0 | arm_id >= nrow(gla)))
    stop("`arm_id` must have one 0-based row index into `gla` per variant ",
         "(length ", n_variants, ", values 0 to ", nrow(gla) - 1L, ")",
         call. = FALSE)
  invisible(NULL)
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
    .check_matrix(gt_genotype, values = 0:2)
    .check_matrix(ancestry)   # unrecognised codes are allowed: they contribute nothing
    .check_same_dim(gt_genotype, ancestry)
    .check_shrink_target(gla, arm_id, nrow(gt_genotype), 2L)
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
#'   Example: \code{c(AFR = 3L, EUR = 1L, NAT = 4L)}. The default is the
#'   two-population coding used by \code{\link{split_diploid}}
#'   (3 = AFR/AFR, 1 = EUR/EUR).
#' @param mixed_codes Data frame with three columns:
#'   \describe{
#'     \item{code}{Integer ancestry code for this mixed pair.}
#'     \item{pop1}{Character name of the first parent population
#'       (must be a name in \code{pure_codes}).}
#'     \item{pop2}{Character name of the second parent population
#'       (must be a name in \code{pure_codes}).}
#'   }
#'   One row per mixed-ancestry diploid type. The default is the
#'   two-population mixed code 2 = AFR/EUR, so with both defaults
#'   \code{split_diploid_multi(gt, ancestry)} matches
#'   \code{split_diploid(gt, ancestry)}.
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
                                pure_codes = c(AFR = 3L, EUR = 1L),
                                mixed_codes = data.frame(code = 2L,
                                                         pop1 = "AFR",
                                                         pop2 = "EUR"),
                                gla = NULL, arm_id = NULL) {
    .check_matrix(gt_genotype, values = 0:2)
    .check_matrix(ancestry)
    .check_same_dim(gt_genotype, ancestry)
    .check_codes(pure_codes)
    if (!is.data.frame(mixed_codes) ||
        !all(c("code", "pop1", "pop2") %in% names(mixed_codes)))
        stop("mixed_codes must be a data.frame with columns: code, pop1, pop2")
    pop_names <- names(pure_codes)
    bad <- setdiff(c(mixed_codes$pop1, mixed_codes$pop2), pop_names)
    if (length(bad))
        stop("mixed_codes references populations not in pure_codes: ",
             toString(bad))
    if (!is.null(gla) && !identical(colnames(gla), pop_names))
        stop("gla's columns must match names(pure_codes) exactly, in order: ",
             toString(pop_names))
    .check_shrink_target(gla, arm_id, nrow(gt_genotype), length(pop_names))

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
    .check_haplotype_inputs(gt_hap0, gt_hap1, anc_hap0, anc_hap1)
    .check_codes(pop_codes, min_length = 2L)
    if (length(pop_codes) != 2L)
        stop("`pop_codes` must have exactly 2 populations; use ",
             "split_haplotype_multi() for more", call. = FALSE)
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
    .check_haplotype_inputs(gt_hap0, gt_hap1, anc_hap0, anc_hap1)
    .check_codes(pop_codes, min_length = 2L)
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

.check_haplotype_inputs <- function(gt_hap0, gt_hap1, anc_hap0, anc_hap1) {
    .check_matrix(gt_hap0, values = 0:1)
    .check_matrix(gt_hap1, values = 0:1)
    .check_matrix(anc_hap0)
    .check_matrix(anc_hap1)
    .check_same_dim(gt_hap0, gt_hap1)
    .check_same_dim(gt_hap0, anc_hap0)
    .check_same_dim(gt_hap0, anc_hap1)
    invisible(NULL)
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
#'   column 2). Defaults to \code{bed} with its extension replaced by
#'   \code{.bim}, the usual PLINK file set layout.
#' @param fam Path to the \code{.fam} file (sample IDs are taken from
#'   column 2, i.e. IID). Defaults to \code{bed} with its extension replaced
#'   by \code{.fam}.
#'
#' @return Integer matrix (variants x samples), dimnamed from the
#'   \code{.bim}/\code{.fam} files. Values: 0 = no call, 1 = homozygous
#'   first allele, 2 = heterozygous, 3 = homozygous second allele.
#'
#' @seealso \code{\link{split_diploid}}, \code{\link{ancestry_split_dosage}}
#'
#' @examples
#' ext <- system.file("extdata", package = "lantern")
#' # .bim and .fam are found next to the .bed file
#' pt <- read_bed_file(file.path(ext, "toy_ancestry.bed"))
#' pt[, 1:6]   # ancestry tracts x samples, codes 1/2/3
#'
#' @export
read_bed_file <- function(bed, bim = sub("\\.bed$", ".bim", bed),
                          fam = sub("\\.bed$", ".fam", bed)) {
  .check_file(bed)
  .check_file(bim)
  .check_file(fam)
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

  pop_codes  <- .parse_msp_pop_codes(readLines(con, n = 1))
  headers    <- strsplit(readLines(con, n = 1), "\t")[[1]]
  sample_ids <- unique(sub("\\.[01]$", "", headers[7:length(headers)]))

  if (verbose) {
    message("    Found ", length(sample_ids), " samples, ", length(pop_codes), " populations")
    message("    Pop codes: ", toString(sprintf("%s = %d", names(pop_codes), pop_codes)))
  }

  data_lines <- readLines(con)
  if (length(data_lines) == 0) stop("MSP file has no data rows")
  data_mat <- do.call(rbind, strsplit(data_lines, "\t", fixed = TRUE))

  tract_df <- data.frame(
    chrom = data_mat[, 1],
    spos  = as.integer(data_mat[, 2]),
    epos  = as.integer(data_mat[, 3]),
    stringsAsFactors = FALSE
  )
  tract_order <- order(tract_df$spos)
  tract_df    <- tract_df[tract_order, ]
  rownames(tract_df) <- NULL

  # Data rows: 6 metadata fields (chm spos epos sgpos egpos "n snps"), then
  # one ancestry call per haplotype, ordered sample1.0, sample1.1, ...
  n_samples <- length(sample_ids)
  anc <- data_mat[tract_order, 6L + seq_len(2L * n_samples), drop = FALSE]
  storage.mode(anc) <- "integer"
  hap0 <- seq(1L, by = 2L, length.out = n_samples)
  anc_hap0_mat <- t(anc[, hap0, drop = FALSE])
  anc_hap1_mat <- t(anc[, hap0 + 1L, drop = FALSE])
  rownames(anc_hap0_mat) <- rownames(anc_hap1_mat) <- sample_ids

  list(pop_codes  = pop_codes,
       sample_ids = sample_ids,
       tract_df   = tract_df,
       anc_hap0   = anc_hap0_mat,
       anc_hap1   = anc_hap1_mat)
}

# Population codes from an MSP file's first line, e.g.
# "#Subpopulation order/codes: AFR=0\tEUR=1" -> c(AFR = 0L, EUR = 1L).
.parse_msp_pop_codes <- function(line1) {
  line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
  parts <- strsplit(line1, if (grepl("\t", line1)) "\t" else "\\s+")[[1]]
  m     <- regmatches(parts, regexec("([A-Za-z0-9_]+)\\s*=\\s*([0-9]+)", parts))
  m     <- m[lengths(m) == 3L]
  if (length(m) == 0) {
    warning("Could not parse pop codes from MSP; using default AFR=0, EUR=1")
    return(c(AFR = 0L, EUR = 1L))
  }
  codes <- stats::setNames(as.integer(vapply(m, `[`, "", 3L)),
                           vapply(m, `[`, "", 2L))
  codes[!duplicated(names(codes), fromLast = TRUE)]
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

  gt <- .gds_gt_strings(gds, n_samples, n_variants)
  colnames(gt) <- sample_ids
  list(chrom = vcf_chrom,
       pos   = SeqArray::seqGetData(gds, "position"),
       ref   = SeqArray::seqGetData(gds, "$ref"),
       alt   = alt,
       gt    = gt)
}

# Variants x samples matrix of GT strings, as bcftools prints them, from an
# open GDS. genotype: ploidy x sample x variant allele indices (NA =
# missing); phase: sample x variant, 1 = phased ("|"), 0 = unphased ("/").
.gds_gt_strings <- function(gds, n_samples, n_variants) {
  geno <- SeqArray::seqGetData(gds, "genotype")
  geno <- array(geno, dim = c(length(geno) / (n_samples * n_variants),
                              n_samples, n_variants))
  allele_str <- function(a) {
    s <- as.character(a)
    s[is.na(a)] <- "."
    s
  }
  gt <- if (dim(geno)[1] == 1) {
    allele_str(geno[1, , ])
  } else {
    paste0(allele_str(geno[1, , ]),
           ifelse(as.vector(SeqArray::seqGetData(gds, "phase")) == 1, "|", "/"),
           allele_str(geno[2, , ]))
  }
  t(matrix(gt, nrow = n_samples, ncol = n_variants))
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
  msp <- .restrict_tracts(.parse_msp(msp_path, verbose = verbose), chrom)
  if (verbose) message("  MSP samples: ", length(msp$sample_ids),
                       "  Tracts: ", nrow(msp$tract_df),
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
  common_samples <- intersect(vcf_samples, msp$sample_ids)
  if (length(common_samples) == 0)
    stop("No common samples between VCF and MSP")
  if (verbose) message("  Common: ", length(common_samples))

  # ---- Step 4: Parse VCF GT ----
  if (verbose) message("\nStep 4: Parsing VCF genotypes...")
  geno <- .read_vcf_haplotypes(vcf_path, chrom, reader,
                               match(common_samples, vcf_samples),
                               common_samples, verbose)

  # ---- Step 5: Map variants to tracts ----
  if (verbose) message("\nStep 5: Mapping variants to ancestry tracts...")
  geno <- .map_variants_to_tracts(geno, msp$tract_df, chrom, verbose)

  # ---- Step 6: Broadcast tract ancestry to variants ----
  # Each variant takes its tract's per-sample ancestry calls.
  if (verbose) message("\nStep 6: Broadcasting ancestry to variants...")
  anc_hap0_common <- msp$anc_hap0[common_samples, , drop = FALSE]
  anc_hap1_common <- msp$anc_hap1[common_samples, , drop = FALSE]

  list(
    gt_hap0_mat     = geno$hap0,
    gt_hap1_mat     = geno$hap1,
    anc_hap0_var    = unname(t(anc_hap0_common[, geno$tract_idx, drop = FALSE])),
    anc_hap1_var    = unname(t(anc_hap1_common[, geno$tract_idx, drop = FALSE])),
    anc_hap0_tract  = anc_hap0_common,
    anc_hap1_tract  = anc_hap1_common,
    variant_info    = data.frame(chrom = geno$chrom, pos = geno$pos,
                                 ref = geno$ref, alt = geno$alt,
                                 stringsAsFactors = FALSE),
    common_samples  = common_samples,
    pop_codes       = msp$pop_codes,
    tract_df        = msp$tract_df,
    overlap = list(
      n_vcf_samples           = length(vcf_samples),
      n_msp_samples           = length(msp$sample_ids),
      n_common                = length(common_samples),
      n_multiallelic_filtered = geno$n_multi,
      n_no_tract              = geno$n_no_tract,
      dropped_samples_vcf     = setdiff(vcf_samples, common_samples),
      dropped_samples_msp     = setdiff(msp$sample_ids, common_samples)
    )
  )
}

# Keep only `chrom`'s tracts (all tracts when chrom is NULL). Besides saving
# work, this matters for correctness: tract_df is sorted by position alone,
# so with several chromosomes findInterval() could otherwise match a variant
# to another chromosome's tract that overlaps it in position.
.restrict_tracts <- function(msp, chrom) {
  if (is.null(chrom)) return(msp)
  keep <- sub("^chr", "", msp$tract_df$chrom) == sub("^chr", "", chrom)
  if (!any(keep))
    stop("No ancestry tracts on chromosome ", chrom, " in MSP file")
  msp$tract_df <- msp$tract_df[keep, , drop = FALSE]
  rownames(msp$tract_df) <- NULL
  msp$anc_hap0 <- msp$anc_hap0[, keep, drop = FALSE]
  msp$anc_hap1 <- msp$anc_hap1[, keep, drop = FALSE]
  msp
}

# Read the VCF's biallelic variants as haplotype allele matrices (variants x
# common samples), dropping multiallelic sites.
.read_vcf_haplotypes <- function(vcf_path, chrom, reader, sample_idx,
                                 common_samples, verbose) {
  vcf <- .vcf_read_gt(vcf_path, chrom, reader, verbose = verbose)
  biallelic <- !grepl(",", vcf$alt, fixed = TRUE)
  n_multi   <- sum(!biallelic)
  if (n_multi > 0 && verbose)
    message("  Skipping ", n_multi, " multiallelic variants")

  gt <- .parse_phased_gt_matrix(vcf$gt[biallelic, sample_idx, drop = FALSE],
                                sum(biallelic), length(common_samples),
                                common_samples, verbose = verbose)
  n_missing <- sum(gt$missing)
  if (n_missing > 0 && verbose)
    message("  Missing GT: ", n_missing, " (treated as 0/0)")

  list(chrom = vcf$chrom[biallelic], pos = vcf$pos[biallelic],
       ref = vcf$ref[biallelic], alt = vcf$alt[biallelic],
       hap0 = gt$hap0, hap1 = gt$hap1, n_multi = n_multi)
}

# Subset the per-variant fields and haplotype rows of .read_vcf_haplotypes()
# output.
.subset_variants <- function(geno, keep) {
  fields <- intersect(c("chrom", "pos", "ref", "alt", "tract_idx"), names(geno))
  geno[fields] <- lapply(geno[fields], `[`, keep)
  geno$hap0 <- geno$hap0[keep, , drop = FALSE]
  geno$hap1 <- geno$hap1[keep, , drop = FALSE]
  geno
}

# Find each variant's ancestry tract (adds `tract_idx`, a row of tract_df)
# and drop variants that are off `chrom` or fall outside every tract.
.map_variants_to_tracts <- function(geno, tract_df, chrom, verbose) {
  vcf_chrom <- sub("^chr", "", geno$chrom)
  if (!is.null(chrom)) {
    keep <- vcf_chrom == sub("^chr", "", chrom)
    if (!any(keep)) stop("No variants on chromosome ", chrom)
    geno      <- .subset_variants(geno, keep)
    vcf_chrom <- vcf_chrom[keep]
    if (verbose) message("  Filtered to chrom ", chrom, ": ", sum(keep), " variants")
  }

  # Last tract starting at or before the variant; valid if the variant is
  # also before that tract's end and on the same chromosome.
  tract_idx <- findInterval(geno$pos, tract_df$spos)
  hit       <- ifelse(tract_idx > 0L, tract_idx, NA_integer_)
  in_tract  <- !is.na(hit) &
    geno$pos <= tract_df$epos[hit] &
    vcf_chrom == sub("^chr", "", tract_df$chrom)[hit]
  in_tract[is.na(in_tract)] <- FALSE

  geno$tract_idx  <- tract_idx
  geno$n_no_tract <- sum(!in_tract)
  if (geno$n_no_tract > 0) {
    if (verbose) message("  Variants with no tract: ", geno$n_no_tract, " (dropped)")
    geno <- .subset_variants(geno, in_tract)
  }
  if (length(geno$pos) == 0) stop("No variants with valid ancestry tracts")
  geno
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
  .check_file(vcf_path)
  .check_file(msp_path)
  if (is.numeric(chrom)) chrom <- as.character(chrom)   # allow chrom = 19
  .check_string(chrom, null_ok = TRUE)
  .check_flag(verbose)

  if (verbose) message("=== LANTERN Ancestry Split (Step 1, mode = ", mode, ") ===\n")

  common    <- .parse_vcf_msp_common(vcf_path, msp_path, chrom = chrom, verbose = verbose)
  pop_codes <- common$pop_codes
  if (length(pop_codes) < 2)
    stop("ancestry_split() requires at least 2 populations in the MSP")
  codes       <- .diploid_codes(names(pop_codes))
  anc_diploid <- .to_diploid_codes(common$anc_hap0_var, common$anc_hap1_var,
                                   pop_codes, codes$lookup)

  # ---- Split by ancestry ----
  # Haplotype mode is deterministic. Dosage mode shrinks each ambiguous
  # pair's split toward 1/2 (a uniform 1/K target shared by all variants;
  # see the "Shrinkage toward 1/2" section).
  if (verbose) message("\nStep 7: Splitting genotypes by ancestry (", mode, ")...")
  res <- if (mode == "haplotype") {
    split_haplotype_multi(common$gt_hap0_mat, common$gt_hap1_mat,
                          common$anc_hap0_var, common$anc_hap1_var,
                          pop_codes = pop_codes)
  } else {
    split_diploid_multi(common$gt_hap0_mat + common$gt_hap1_mat, anc_diploid,
                        codes$pure, codes$mixed,
                        gla = .uniform_shrink_target(names(pop_codes)),
                        arm_id = integer(nrow(anc_diploid)))
  }

  # Per-variant pure-ancestry sample counts (feed Step 3's gene weights)
  ancestry_counts <- .pure_ancestry_counts(anc_diploid, codes$pure)

  # ---- Filter monomorphic variants ----
  if (verbose) message("\nStep 8: Filtering monomorphic variants...")
  has_alt <- rowSums(Reduce(`+`, res)) > 0
  n_mono  <- sum(!has_alt)
  if (n_mono > 0) {
    if (verbose) message("  Removed ", n_mono, " monomorphic variants")
    res <- lapply(res, function(m) m[has_alt, , drop = FALSE])
  }
  variant_info    <- common$variant_info[has_alt, , drop = FALSE]
  ancestry_counts <- ancestry_counts[has_alt, , drop = FALSE]
  if (nrow(variant_info) == 0) stop("All variants are monomorphic after filtering")
  rownames(variant_info) <- NULL

  if (verbose) message("\n=== Step 1 Complete: ", nrow(variant_info), " variants, ",
                       length(common$common_samples), " samples, ",
                       length(pop_codes), " populations ===\n")

  invisible(c(res, list(
    variant_info    = variant_info,
    sample_ids      = common$common_samples,
    mode            = mode,
    ancestry_counts = ancestry_counts,
    tract_info      = common$tract_df,
    overlap         = c(common$overlap,
                        n_variants_kept        = nrow(variant_info),
                        n_monomorphic_filtered = n_mono)
  )))
}

# Diploid ancestry coding for K populations: population p's pure code is p
# (1..K), and the K*(K-1)/2 mixed pairs get codes K+1, K+2, ...
# `lookup[a, b]` is the diploid code of a sample whose haplotypes carry
# populations a and b.
.diploid_codes <- function(pop_names) {
  K     <- length(pop_names)
  pairs <- utils::combn(K, 2)
  mixed <- K + seq_len(ncol(pairs))
  lookup <- diag(seq_len(K), K)
  lookup[t(pairs)] <- mixed
  lookup[t(pairs[2:1, , drop = FALSE])] <- mixed
  storage.mode(lookup) <- "integer"
  list(pure   = stats::setNames(seq_len(K), pop_names),
       mixed  = data.frame(code = as.integer(mixed),
                           pop1 = pop_names[pairs[1L, ]],
                           pop2 = pop_names[pairs[2L, ]],
                           stringsAsFactors = FALSE),
       lookup = lookup)
}

# Combine two haplotype ancestry matrices (MSP population codes) into one
# diploid ancestry-code matrix.
.to_diploid_codes <- function(anc_hap0, anc_hap1, pop_codes, lookup) {
  p0 <- match(as.vector(anc_hap0), pop_codes)
  p1 <- match(as.vector(anc_hap1), pop_codes)
  matrix(lookup[cbind(p0, p1)], nrow = nrow(anc_hap0), ncol = ncol(anc_hap0))
}

# Number of samples with pure ancestry in each population, per variant.
.pure_ancestry_counts <- function(anc_diploid, pure_codes) {
  matrix(vapply(pure_codes, function(code) count_ancestry_codes(anc_diploid, code),
                numeric(nrow(anc_diploid))),
         nrow = nrow(anc_diploid), ncol = length(pure_codes),
         dimnames = list(NULL, names(pure_codes)))
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
#' @section Matching variants to regions: Each variant (row of
#'   \code{gt_matrix}) is paired with an ancestry column of \code{pt_matrix}
#'   by, in order of preference: identical names; position, when rownames
#'   look like \code{"chr:pos"} and colnames like \code{"chr:start-end"}
#'   (the first region containing the variant; several variants may share a
#'   region); or, when either matrix has no names, position in the matrix.
#'
#' @return List with elements:
#'   \item{african}{African ancestry-specific dosage matrix (variants x
#'     samples)}
#'   \item{european}{European ancestry-specific dosage matrix}
#'   \item{counts}{List (\code{african}, \code{european}, \code{mixed}) of
#'     per-variant counts of samples with ancestry code 3, 1 and 2}
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
  .check_flag(verbose)
  if (!is.null(vcf_path) || !is.null(msp_path)) {
    if (is.null(vcf_path) || is.null(msp_path))
      stop("Supply both `vcf_path` and `msp_path`", call. = FALSE)
    return(.split_dosage_from_vcf(vcf_path, msp_path, chrom, verbose))
  }
  if (is.null(gt_matrix) || is.null(pt_matrix))
    stop("Provide either (gt_matrix, pt_matrix) or (vcf_path, msp_path)")
  .check_matrix(gt_matrix, values = 0:2)
  .check_matrix(pt_matrix, values = 0:3)

  if (verbose) message("=== LANTERN Ancestry Pipeline ===\n")
  if (verbose) message("Step 1: Finding sample overlap...")
  smp <- .match_samples(gt_matrix, pt_matrix, verbose)
  if (verbose) message("\nStep 2: Finding variant/region overlap...")
  var <- .match_variants(smp$gt, smp$pt, verbose)

  # pt_matrix is samples x regions; the C splitter reads ancestry with the
  # genotype matrix's layout, so build variants x samples, one column of
  # pt_matrix per variant.
  if (verbose) message("\nStep 3: Subsetting matrices...")
  gt  <- smp$gt[var$gt_idx, smp$gt_idx, drop = FALSE]
  anc <- t(smp$pt[smp$pt_idx, var$pt_idx, drop = FALSE])

  if (verbose) message("\nStep 3b: Filtering monomorphic variants...")
  has_alt <- rowSums(gt > 0) > 0
  n_mono  <- sum(!has_alt)
  if (n_mono > 0 && verbose)
    message("  Removed ", n_mono, " monomorphic variants (no alt alleles)")
  dropped_variants <- c(var$dropped, var$gt_names[var$gt_idx][!has_alt])
  gt  <- gt[has_alt, , drop = FALSE]
  anc <- anc[has_alt, , drop = FALSE]
  if (nrow(gt) == 0) stop("All variants are monomorphic (no alt alleles)")

  # Same shrinkage toward 1/2 as ancestry_split()
  if (verbose) message("\nStep 4: Splitting genotypes by ancestry...")
  result <- split_diploid(gt, anc, gla = .uniform_shrink_target(c("AFR", "EUR")),
                          arm_id = integer(nrow(gt)))

  if (verbose) message("\nStep 5: Counting ancestries per variant...")
  counts <- list(african  = count_ancestry_codes(anc, 3),
                 european = count_ancestry_codes(anc, 1),
                 mixed    = count_ancestry_codes(anc, 2))
  if (verbose) message("\n=== Pipeline Complete ===\n")

  invisible(list(
    african  = result$african,
    european = result$european,
    counts   = counts,
    overlap  = list(
      n_samples_total        = ncol(gt_matrix) + nrow(pt_matrix) - length(smp$common),
      n_samples_kept         = length(smp$common),
      n_variants_total       = nrow(gt_matrix),
      n_variants_kept        = nrow(gt),
      n_monomorphic_filtered = n_mono,
      dropped_samples        = unique(smp$dropped),
      dropped_variants       = unique(dropped_variants)
    )
  ))
}

# ancestry_split_dosage(vcf_path =, msp_path =): run ancestry_split() and
# return its two populations in this function's african/european layout.
.split_dosage_from_vcf <- function(vcf_path, msp_path, chrom, verbose) {
  split_result <- ancestry_split(vcf_path, msp_path, mode = "dosage",
                                 chrom = chrom, verbose = verbose)
  pops <- .afr_eur_names(.split_pop_names(split_result), "ancestry_split_dosage()")
  counts <- split_result$ancestry_counts
  invisible(list(
    african  = split_result[[pops[["afr"]]]],
    european = split_result[[pops[["eur"]]]],
    counts   = list(african  = counts[, pops[["afr"]]],
                    european = counts[, pops[["eur"]]],
                    mixed    = length(split_result$sample_ids) -
                      counts[, pops[["afr"]]] - counts[, pops[["eur"]]]),
    overlap  = split_result$overlap
  ))
}

# Names of the per-population dosage matrices in an ancestry_split() result.
.split_pop_names <- function(split_result) {
  setdiff(names(split_result),
          c("variant_info", "sample_ids", "mode", "ancestry_counts",
            "tract_info", "overlap"))
}

# Pick the African and European population of a two-population split
# (AFR/EUR, or AFRICAN/EUROPEAN in any case), else the first and second.
.afr_eur_names <- function(pop_names, caller) {
  if (length(pop_names) != 2)
    stop(caller, " only supports 2-population (AFR/EUR) MSP files; use ",
         "ancestry_split() directly for K > 2 populations.", call. = FALSE)
  afr <- pop_names[toupper(pop_names) %in% c("AFR", "AFRICAN")]
  eur <- pop_names[toupper(pop_names) %in% c("EUR", "EUROPEAN")]
  if (length(afr) != 1 || length(eur) != 1) {
    afr <- pop_names[1]
    eur <- pop_names[2]
  }
  c(afr = afr, eur = eur)
}

# Give an unnamed axis positional names, keeping the first min(n_gt, n_pt)
# entries when the two lengths differ.
.name_by_position <- function(n_gt, n_pt, prefix, what) {
  n <- min(n_gt, n_pt)
  if (n_gt != n_pt)
    warning("Matrices have different numbers of ", what, " (GT: ", n_gt,
            ", PT: ", n_pt, "). Using the first ", n, ".", call. = FALSE)
  paste0(prefix, seq_len(n))
}

# Samples present in both matrices (GT columns, PT rows), matched by name.
.match_samples <- function(gt, pt, verbose) {
  if (is.null(colnames(gt)) || is.null(rownames(pt))) {
    ids <- .name_by_position(ncol(gt), nrow(pt), "sample_", "samples")
    gt  <- gt[, seq_along(ids), drop = FALSE]
    pt  <- pt[seq_along(ids), , drop = FALSE]
    colnames(gt) <- rownames(pt) <- ids
  }
  common <- intersect(colnames(gt), rownames(pt))
  if (length(common) == 0) stop("No common samples found between GT and PT matrices")
  dropped <- c(setdiff(colnames(gt), common), setdiff(rownames(pt), common))
  if (verbose) {
    message("  Samples: ", ncol(gt), " in GT, ", nrow(pt), " in PT")
    message("  Common samples: ", length(common))
    if (length(dropped) > 0)
      message("  Dropped: ", toString(utils::head(dropped, 3)),
              if (length(dropped) > 3) "...")
  }
  list(gt = gt, pt = pt, common = common, dropped = dropped,
       gt_idx = match(common, colnames(gt)), pt_idx = match(common, rownames(pt)))
}

# Pair each GT variant (row) with a PT region (column): see the "Matching
# variants to regions" section of ?ancestry_split_dosage. Returns row/column
# indices of equal length (pt_idx may repeat) and the unmatched names.
.match_variants <- function(gt, pt, verbose) {
  if (is.null(rownames(gt)) || is.null(colnames(pt))) {
    ids <- .name_by_position(nrow(gt), ncol(pt), "var_", "variants")
    rownames(gt) <- NULL
    colnames(pt) <- NULL
    gt_names <- c(ids, character(nrow(gt) - length(ids)))
    pt_names <- c(ids, character(ncol(pt) - length(ids)))
  } else {
    gt_names <- rownames(gt)
    pt_names <- colnames(pt)
  }

  pt_idx <- match(gt_names, pt_names)
  pt_idx[!nzchar(gt_names)] <- NA
  if (all(is.na(pt_idx))) {
    if (verbose) message("  No exact matches. Trying coordinate-based matching...")
    pt_idx <- .match_by_coordinates(gt_names, pt_names)
    if (all(is.na(pt_idx)))
      stop("No variant overlap between rownames(gt_matrix) and ",
           "colnames(pt_matrix), by name or by \"chr:pos\" within ",
           "\"chr:start-end\" coordinates", call. = FALSE)
    if (verbose) message("  Found ", sum(!is.na(pt_idx)), " variants in ",
                         length(unique(stats::na.omit(pt_idx))), " regions")
  }
  gt_idx  <- which(!is.na(pt_idx))
  pt_idx  <- pt_idx[gt_idx]
  dropped <- c(gt_names[-gt_idx], pt_names[-unique(pt_idx)])
  dropped <- dropped[nzchar(dropped)]
  if (verbose) {
    message("  Variants in GT: ", nrow(gt))
    message("  Regions in PT: ", ncol(pt))
    message("  Matched variants: ", length(gt_idx))
    if (length(dropped) > 0)
      message("  Dropped: ", toString(utils::head(dropped, 3)),
              if (length(dropped) > 3) "...")
  }
  list(gt_idx = gt_idx, pt_idx = pt_idx, dropped = dropped, gt_names = gt_names)
}

# For each "chr:pos" variant name, the index of the first "chr:start-end"
# region name containing it (NA if none or unparsable).
.match_by_coordinates <- function(gt_names, pt_names) {
  field2 <- function(x) ifelse(grepl(":", x, fixed = TRUE),
                               sub("^[^:]*:([^:]*).*$", "\\1", x), NA_character_)
  gt_chr <- sub(":.*$", "", gt_names)
  gt_pos <- suppressWarnings(as.numeric(field2(gt_names)))

  range_str <- field2(pt_names)
  is_range  <- grepl("^[^-]+-[^-]+$", range_str)
  pt_chr    <- sub(":.*$", "", pt_names)
  pt_start  <- suppressWarnings(as.numeric(ifelse(is_range, sub("-.*$", "", range_str), NA)))
  pt_end    <- suppressWarnings(as.numeric(ifelse(is_range, sub("^.*-", "", range_str), NA)))

  vapply(seq_along(gt_names), function(i) {
    if (is.na(gt_pos[i])) return(NA_integer_)
    hit <- which(pt_chr == gt_chr[i] & pt_start <= gt_pos[i] & pt_end >= gt_pos[i])
    if (length(hit)) hit[1] else NA_integer_
  }, integer(1))
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
#' @param out_path Output directory for the GDS files. Only needed (and
#'   then required) when \code{write_vcf = TRUE}.
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
ancestry_split_phased <- function(vcf_path, msp_path, out_path = NULL,
                                chrom = NULL, write_vcf = TRUE,
                                verbose = TRUE) {
  .check_flag(write_vcf)
  .check_flag(verbose)
  if (write_vcf) .check_string(out_path)

  if (verbose) message("=== LANTERN Phased Ancestry Pipeline ===\n")
  split_result <- ancestry_split(vcf_path, msp_path, mode = "haplotype",
                                 chrom = chrom, verbose = verbose)
  pops <- .afr_eur_names(.split_pop_names(split_result), "ancestry_split_phased()")

  vcf_paths <- NULL
  if (write_vcf) {
    gds_paths <- write_ancestry_gds(split_result, out_path, verbose = verbose)
    vcf_paths <- list(african_gds  = unname(gds_paths[[pops[["afr"]]]]),
                      european_gds = unname(gds_paths[[pops[["eur"]]]]))
  }
  if (verbose) message("\n=== Pipeline Complete ===\n")

  invisible(list(african      = split_result[[pops[["afr"]]]],
                 european     = split_result[[pops[["eur"]]]],
                 variant_info = split_result$variant_info,
                 sample_ids   = split_result$sample_ids,
                 tract_info   = split_result$tract_info,
                 overlap      = split_result$overlap,
                 vcf_paths    = vcf_paths))
}
