# ============================================================================
# ancestry_smmat (Step 3): SMMAT association testing + Cauchy combination
# ============================================================================

#' Per-gene ancestry weights for the Cauchy combination (Step 3 helper)
#'
#' Generalises the retired \code{src/step3_weight_finding.R} script's
#' median per-gene ancestry count logic to K populations (that
#' functionality now lives here, fed by \code{ancestry_counts}/
#' \code{variant_info} from Step 1): for each gene in \code{gene_group_file},
#' takes the median (across the gene's variants) of the per-variant pure-
#' ancestry sample counts in \code{ancestry_counts}.
#'
#' @param ancestry_counts Numeric matrix (variants x populations), as
#'   returned by \code{\link{ancestry_split}}.
#' @param variant_info data.frame with columns chrom/pos, row-aligned to
#'   \code{ancestry_counts}, as returned by \code{\link{ancestry_split}}.
#' @param gene_group_file Path to the gene_group file (no header; columns
#'   gene, chr, pos, ref, alt, weight) -- the same file passed to
#'   \code{GMMAT::SMMAT()}.
#' @param pop_names Character vector of population names to produce weight
#'   columns for (typically \code{names(gds_paths)} from \code{\link{ancestry_smmat}}).
#'
#' @return Numeric matrix (genes x \code{pop_names}), or \code{NULL} if
#'   none of \code{pop_names} match \code{colnames(ancestry_counts)}.
#'
#' @keywords internal
.compute_gene_weights <- function(ancestry_counts, variant_info,
                                   gene_group_file, pop_names) {
  gg <- data.table::fread(gene_group_file, header = FALSE,
                          col.names = c("gene", "chr", "pos", "ref", "alt", "weight"))

  avail_pops <- intersect(pop_names, colnames(ancestry_counts))
  if (length(avail_pops) == 0) {
    warning("None of the population names (", toString(pop_names),
            ") match ancestry_counts columns (",
            toString(colnames(ancestry_counts)),
            "); falling back to equal weights.")
    return(NULL)
  }

  vi_key  <- paste0(sub("^chr", "", as.character(variant_info$chrom)), ":",
                    variant_info$pos)
  gg_key  <- paste0(sub("^chr", "", as.character(gg$chr)), ":", gg$pos)
  row_idx <- match(gg_key, vi_key)
  genes   <- factor(gg$gene, levels = unique(gg$gene))

  # Weight 1 for populations without counts and for genes none of whose
  # variants are in variant_info; otherwise the median count over the
  # gene's variants.
  weights_mat <- matrix(1, nrow = nlevels(genes), ncol = length(pop_names),
                        dimnames = list(levels(genes), pop_names))
  hit <- !is.na(row_idx)
  med <- vapply(avail_pops, function(pop)
    as.numeric(tapply(ancestry_counts[row_idx[hit], pop], genes[hit], stats::median)),
    numeric(nlevels(genes)))
  med <- matrix(med, nrow = nlevels(genes), dimnames = list(levels(genes), avail_pops))
  weights_mat[, avail_pops] <- ifelse(is.na(med), 1, med)
  weights_mat
}

#' Run ancestry-stratified SMMAT + Cauchy combination (Step 3)
#'
#' Wraps \code{GMMAT::SMMAT()} across all ancestry-split GDS files (as
#' produced by \code{\link{write_ancestry_gds}}), computes per-gene
#' ancestry weights from Step 1's \code{ancestry_counts} (this replaces the
#' retired \code{src/step3_weight_finding.R} script), and combines the per-gene p-values
#' across populations with \code{\link{cauchy_combine}}. This is the
#' package's Step 3 entry point; it does not run Step 1/2 for you (SMMAT
#' runs can take hours, so this is always an explicit, separate call).
#'
#' @param gds_paths Named list/vector of GDS file paths, one per ancestry
#'   population, e.g. \code{c(AFR = "afr.gds", EUR = "eur.gds")} -- as
#'   returned by \code{\link{write_ancestry_gds}}. Any K >= 1 populations
#'   are supported; add an extra named entry (e.g. \code{OBSERVED = ...})
#'   to include an additional (e.g. non-ancestry-split) GDS in the same
#'   SMMAT + Cauchy combination run.
#' @param pheno data.frame with an \code{id_col} column plus the response
#'   and covariates referenced in \code{formula}.
#' @param formula Model formula for the null model, e.g. \code{y ~ age + sex}.
#' @param kinship Numeric matrix (or object coercible via \code{as.matrix})
#'   with row/column names equal to sample IDs.
#' @param gene_group_file Either a path to the gene_group file passed to
#'   \code{GMMAT::SMMAT()} (no header; columns gene, chr, pos, ref, alt,
#'   weight), or a data.frame/matrix with those same six columns. Either way,
#'   it is read, has any leading \code{"chr"} stripped from the chr column
#'   (matching \code{SeqArray::seqVCF2GDS()}'s convention -- see the "chr
#'   column" note below), and written to a fresh temporary file that is
#'   deleted with \code{unlink()} when \code{ancestry_smmat()} returns; the
#'   caller's own path/data.frame is never read-modified or deleted, and can
#'   safely be reused across repeated calls (e.g. one call per phenotype).
#'   \strong{chr column}: pass either \code{"19"} or \code{"chr19"} --
#'   normalized automatically. Without this, a mismatch between the group
#'   file's chr convention and the GDS's (always bare, no \code{"chr"}
#'   prefix) would silently drop every affected variant, and the whole gene
#'   if all its variants are affected, since \code{GMMAT::SMMAT()} matches
#'   variants by exact \code{chr:pos:ref:alt} string equality with no
#'   fuzzy/positional alignment.
#' @param ancestry_counts,variant_info Optional: pass
#'   \code{ancestry_split()$ancestry_counts} / \code{$variant_info}
#'   straight through to compute per-gene ancestry weights for the Cauchy
#'   combination. If omitted, equal weights are used.
#' @param id_col Name of the sample ID column in \code{pheno} and
#'   \code{formula}'s \code{id} argument to \code{glmmkin} (default \code{"id"}).
#' @param family A family object for \code{GMMAT::glmmkin()} (default
#'   \code{gaussian()}).
#' @param gene_col Name of the gene/group ID column in \code{SMMAT()}'s
#'   output (default \code{"group"}, matching \code{GMMAT::SMMAT()}).
#' @param p_col Name of the p-value column in \code{SMMAT()}'s output to
#'   combine (default \code{"E.pval"}, matching \code{GMMAT::SMMAT()}'s
#'   default \code{tests = "E"}; use \code{"O.pval"}/\code{"S.pval"}/
#'   \code{"B.pval"} if you pass a different \code{tests} via \code{...}).
#' @param MAF.range,miss.cutoff,method,ncores Passed straight through to
#'   \code{GMMAT::SMMAT()} -- same names, same meaning as in a direct
#'   \code{SMMAT()} call.
#' @param verbose Print this function's own step-by-step progress messages
#'   (Steps 1-6 below). Distinct from \code{GMMAT::SMMAT()}'s own internal
#'   \code{verbose} (its per-call progress bar), which is not exposed
#'   separately and stays at SMMAT's default (\code{FALSE}); every other
#'   \code{SMMAT()} argument is available via \code{...}.
#' @param ... Any other \code{GMMAT::SMMAT()} argument, forwarded as-is --
#'   for example \code{tests} (\code{"B"}, \code{"S"}, \code{"O"}, or the
#'   default \code{"E"}), \code{rho}, \code{MAF.weights.beta},
#'   \code{missing.method}, \code{use.minor.allele}, \code{auto.flip},
#'   \code{Garbage.Collection}, \code{group.file.sep}, or
#'   \code{meta.file.prefix}. \code{null.obj}, \code{geno.file},
#'   \code{group.file}, and \code{is.dosage} are always supplied internally
#'   by \code{ancestry_smmat()} and cannot be overridden this way. See
#'   \code{?GMMAT::SMMAT} for the full list and defaults.
#'
#' @return A list with two elements:
#'   \item{results}{data.frame with one row per gene: \code{gene}, one
#'     \code{p_<POP>} column per population, one \code{w_<POP>} weight
#'     column per population (only if \code{ancestry_counts}/
#'     \code{variant_info} were supplied), and \code{p_cauchy}.}
#'   \item{smmat_results}{Named list (one entry per population, same names
#'     as \code{gds_paths}) of the raw \code{GMMAT::SMMAT()} result
#'     data.frames.}
#'
#' @seealso \code{\link{ancestry_split}}, \code{\link{write_ancestry_gds}},
#'   \code{\link{cauchy_combine}}
#'
#' @examples
#' if (requireNamespace("GMMAT", quietly = TRUE) &&
#'     requireNamespace("SeqVarTools", quietly = TRUE)) {
#'   ext <- system.file("extdata", package = "lantern")
#'   split <- ancestry_split(file.path(ext, "toy.vcf.gz"),
#'                           file.path(ext, "toy.msp.tsv"),
#'                           mode = "dosage", chrom = "chr19", verbose = FALSE)
#'   out_dir <- tempfile()
#'   gds <- write_ancestry_gds(split, out_dir, verbose = FALSE)
#'
#'   pheno <- read.delim(file.path(ext, "toy_pheno.tsv"))
#'   kinship <- diag(nrow(pheno))   # unrelated samples
#'   dimnames(kinship) <- list(pheno$id, pheno$id)
#'
#'   # gene_group_file: a path, or a data.frame with the same six columns
#'   out <- ancestry_smmat(gds, pheno, y ~ age + sex, kinship,
#'                         file.path(ext, "toy_genes.tsv"),
#'                         ancestry_counts = split$ancestry_counts,
#'                         variant_info    = split$variant_info,
#'                         verbose = FALSE)
#'   out$results              # per-gene p-values + Cauchy combination
#'   head(out$smmat_results$AFR)   # raw GMMAT::SMMAT() output for AFR
#'   unlink(out_dir, recursive = TRUE)
#' }
#'
#' @export
ancestry_smmat <- function(gds_paths, pheno, formula, kinship, gene_group_file,
                            ancestry_counts = NULL, variant_info = NULL,
                            id_col = "id", family = stats::gaussian(),
                            gene_col = "group", p_col = "E.pval",
                            MAF.range = c(0, 0.5), miss.cutoff = 1,
                            method = "davies", ncores = 1, verbose = TRUE, ...) {
  if (!requireNamespace("GMMAT", quietly = TRUE))
    stop("Package 'GMMAT' is required for ancestry_smmat(). Install with: ",
         "install.packages('GMMAT') (not available via conda/Bioconductor).")
  # GMMAT::SMMAT() needs SeqVarTools to read GDS input, but GMMAT only
  # Suggests it.
  if (!requireNamespace("SeqVarTools", quietly = TRUE))
    stop("Package 'SeqVarTools' is required for ancestry_smmat(). Install ",
         "with: BiocManager::install('SeqVarTools')")
  .check_smmat_args(gds_paths, pheno, formula, kinship, ancestry_counts,
                    variant_info, id_col, family, gene_col, p_col, ncores,
                    verbose)

  gene_group_file <- .normalized_gene_groups(gene_group_file)
  on.exit(unlink(gene_group_file), add = TRUE)
  gds_paths <- as.list(gds_paths)
  pop_names <- names(gds_paths)

  if (verbose) message("=== LANTERN Association Testing (Step 3) ===\n")

  # ---- Step 1: Sample overlap across pheno, kinship, and all GDS files ----
  if (verbose) message("Step 1: Intersecting sample IDs...")
  ids <- .smmat_common_samples(gds_paths, pheno, kinship, id_col)
  if (verbose) message("  Common samples: ", length(ids))
  kin_idx <- match(ids, rownames(kinship))

  # ---- Step 2: Fit the null model once ----
  if (verbose) message("\nStep 2: Fitting null model (glmmkin)...")
  model0 <- GMMAT::glmmkin(formula, data = pheno[match(ids, pheno[[id_col]]), , drop = FALSE],
                           kins = as.matrix(kinship[kin_idx, kin_idx, drop = FALSE]),
                           id = id_col, family = family)

  # ---- Step 3: Run SMMAT for each population ----
  smmat_results <- lapply(stats::setNames(seq_along(pop_names), pop_names), function(i) {
    if (verbose) message("\nStep 3: Running SMMAT for ", pop_names[i],
                         " (", i, "/", length(pop_names), ")...")
    GMMAT::SMMAT(model0, gds_paths[[i]], gene_group_file,
                 MAF.range = MAF.range, miss.cutoff = miss.cutoff,
                 method = method, is.dosage = TRUE, ncores = ncores, ...)
  })

  # ---- Step 4: Merge per-gene p-values ----
  if (verbose) message("\nStep 4: Merging per-gene p-values...")
  p_mat <- .merge_smmat_pvalues(smmat_results, gene_col, p_col)

  # ---- Step 5: Per-gene ancestry weights ----
  weights_mat <- NULL
  if (!is.null(ancestry_counts)) {
    if (verbose) message("\nStep 5: Computing per-gene ancestry weights...")
    weights_mat <- .compute_gene_weights(ancestry_counts, variant_info,
                                         gene_group_file, pop_names)
  } else if (verbose) {
    message("\nStep 5: ancestry_counts/variant_info not supplied; ",
            "Cauchy combination will use equal weights.")
  }

  # ---- Step 6: Cauchy-combine p-values per gene ----
  if (verbose) message("\nStep 6: Cauchy-combining p-values per gene...")
  results <- .combine_gene_pvalues(p_mat, weights_mat)
  if (verbose) message("\n=== Step 3 Complete: ", nrow(results), " genes ===\n")

  list(results = results, smmat_results = smmat_results)
}

.check_smmat_args <- function(gds_paths, pheno, formula, kinship,
                              ancestry_counts, variant_info, id_col, family,
                              gene_col, p_col, ncores, verbose) {
  pop_names <- names(gds_paths)
  if (is.null(pop_names) || any(pop_names == "") || anyDuplicated(pop_names))
    stop("gds_paths must be a named list/vector, e.g. c(AFR = ..., EUR = ...)")
  lapply(gds_paths, .check_file, arg = "gds_paths")
  if (!is.data.frame(pheno))
    stop("`pheno` must be a data.frame", call. = FALSE)
  .check_string(id_col)
  if (!id_col %in% names(pheno))
    stop("pheno must contain an id column named '", id_col, "'")
  if (!inherits(formula, "formula"))
    stop("`formula` must be a model formula, e.g. y ~ age + sex", call. = FALSE)
  if (is.null(rownames(kinship)) || is.null(colnames(kinship)))
    stop("kinship must have row and column names corresponding to sample IDs")
  if (nrow(kinship) != ncol(kinship))
    stop("`kinship` must be a square matrix", call. = FALSE)
  if (!inherits(family, "family"))
    stop("`family` must be a family object, e.g. gaussian() or binomial()",
         call. = FALSE)
  if (is.null(ancestry_counts) != is.null(variant_info))
    stop("Supply `ancestry_counts` and `variant_info` together (both from ",
         "ancestry_split()), or neither for equal weights", call. = FALSE)
  if (!is.null(variant_info) && nrow(variant_info) != nrow(ancestry_counts))
    stop("`ancestry_counts` and `variant_info` must have one row per variant",
         call. = FALSE)
  .check_string(gene_col)
  .check_string(p_col)
  .check_count(ncores)
  .check_flag(verbose)
  invisible(NULL)
}

# Gene-group file (path or data.frame) copied to a fresh temp file with any
# leading "chr" stripped from the chr column, so it matches the GDS files:
# SeqArray::seqVCF2GDS() strips "chr", and GMMAT::SMMAT() matches variants
# by exact chr:pos:ref:alt strings, so a "chr19" group file would silently
# drop every variant. The caller's file or data.frame is never modified.
# Returns the temp file's path; the caller deletes it.
.normalized_gene_groups <- function(gene_group_file) {
  if (is.character(gene_group_file)) {
    .check_file(gene_group_file)
    gg_df <- data.table::fread(gene_group_file, header = FALSE)
  } else if (is.data.frame(gene_group_file) || is.matrix(gene_group_file)) {
    gg_df <- as.data.frame(gene_group_file)
  } else {
    stop("gene_group_file must be a file path (character) or a ",
         "data.frame/matrix with columns gene, chr, pos, ref, alt, weight")
  }
  if (ncol(gg_df) != 6)
    stop("gene_group_file must have exactly 6 columns: gene, chr, pos, ref, alt, weight")
  gg_df[[2]] <- sub("^chr", "", as.character(gg_df[[2]]))
  gg_tmp <- tempfile(fileext = ".tsv")
  utils::write.table(gg_df, gg_tmp, sep = "\t", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)
  gg_tmp
}

# Sample IDs present in the phenotype table, the kinship matrix and every
# GDS file.
.smmat_common_samples <- function(gds_paths, pheno, kinship, id_col) {
  gds_ids <- lapply(gds_paths, function(p) {
    g <- SeqArray::seqOpen(p)
    on.exit(SeqArray::seqClose(g), add = TRUE)
    SeqArray::seqGetData(g, "sample.id")
  })
  ids <- Reduce(intersect, c(list(pheno[[id_col]], rownames(kinship)), gds_ids))
  if (length(ids) == 0)
    stop("No overlapping sample IDs found among phenotype, kinship, and GDS files")
  ids
}

# Genes x populations matrix of SMMAT p-values (NA where a population's
# SMMAT run did not test the gene).
.merge_smmat_pvalues <- function(smmat_results, gene_col, p_col) {
  lapply(names(smmat_results), function(pop) {
    cols <- names(smmat_results[[pop]])
    if (!gene_col %in% cols)
      stop("Column '", gene_col, "' not found in SMMAT() output for ", pop)
    if (!p_col %in% cols)
      stop("Column '", p_col, "' not found in SMMAT() output for ", pop,
           "; available columns: ", toString(cols))
  })
  gene_ids <- Reduce(union, lapply(smmat_results, function(r) as.character(r[[gene_col]])))
  p_mat <- vapply(smmat_results, function(r)
    as.numeric(r[[p_col]])[match(gene_ids, as.character(r[[gene_col]]))],
    numeric(length(gene_ids)))
  matrix(p_mat, nrow = length(gene_ids),
         dimnames = list(gene_ids, names(smmat_results)))
}

# One row per gene: p_<POP> columns, w_<POP> weight columns (if weights
# were computed), and the Cauchy-combined p_cauchy (NA when fewer than two
# populations have a p-value for the gene).
.combine_gene_pvalues <- function(p_mat, weights_mat) {
  gene_ids  <- rownames(p_mat)
  pop_names <- colnames(p_mat)
  w_rows    <- if (!is.null(weights_mat)) match(gene_ids, rownames(weights_mat))
  p_cauchy  <- vapply(seq_along(gene_ids), function(i) {
    keep <- !is.na(p_mat[i, ])
    if (sum(keep) < 2) return(NA_real_)
    w <- if (!is.null(weights_mat) && !is.na(w_rows[i]))
      weights_mat[w_rows[i], pop_names][keep]
    cauchy_combine(p_mat[i, keep], weights = w)
  }, numeric(1))

  results <- data.frame(gene = gene_ids, stringsAsFactors = FALSE)
  results[paste0("p_", pop_names)] <- as.data.frame(unname(p_mat))
  if (!is.null(weights_mat))
    results[paste0("w_", pop_names)] <-
      as.data.frame(unname(weights_mat[w_rows, pop_names, drop = FALSE]))
  results$p_cauchy <- p_cauchy
  rownames(results) <- NULL
  results
}
