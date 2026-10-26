# ============================================================================
# Argument checks of the exported functions
# ============================================================================

test_that("file and flag arguments are validated", {
  ext <- system.file("extdata", package = "lantern")
  msp <- file.path(ext, "toy.msp.tsv")
  expect_error(ancestry_split("no-such.vcf", msp), "`vcf_path` does not exist")
  expect_error(ancestry_split(c("a", "b"), msp), "`vcf_path` must be a single")
  expect_error(ancestry_split(file.path(ext, "toy.vcf.gz"), msp, verbose = NA),
               "`verbose` must be TRUE or FALSE")
  expect_error(read_bed_file("no-such.bed"), "`bed` does not exist")
  expect_error(ancestry_split_phased(file.path(ext, "toy.vcf.gz"), msp,
                                     write_vcf = TRUE, verbose = FALSE),
               "`out_path` must be a single")
})

test_that("split primitives validate matrix values and dimensions", {
  gt <- matrix(c(0, 1, 2, 1), 2, 2)
  an <- matrix(c(3, 2, 1, 3), 2, 2)
  expect_error(split_diploid(gt + 1, an), "`gt_genotype` must only contain")
  expect_error(split_diploid(gt, an[, 1, drop = FALSE]), "same dimensions")
  expect_error(split_diploid(gt, an, gla = matrix(.5, 1, 2)),
               "must be supplied together")
  expect_error(split_diploid(gt, an, gla = matrix(.5, 1, 2), arm_id = 0L),
               "one 0-based row index into `gla` per variant")
  expect_error(split_diploid_multi(gt, an, pure_codes = c(3L, 1L)),
               "`pure_codes` must be a named vector")
  h <- matrix(c(0L, 1L, 1L, 0L), 2, 2)
  expect_error(split_haplotype(h + 1L, h, h, h), "`gt_hap0` must only contain")
  expect_error(split_haplotype(h, h, h, h, pop_codes = c(A = 0L, B = 1L, C = 2L)),
               "exactly 2 populations")
})

test_that("split_diploid_multi defaults reproduce split_diploid", {
  set.seed(3)
  gt <- matrix(sample(0:2, 40, TRUE, c(.6, .3, .1)), 8, 5)
  an <- matrix(sample(1:3, 40, TRUE), 8, 5)
  two  <- split_diploid(gt, an)
  mult <- split_diploid_multi(gt, an)
  expect_equal(mult$AFR, two$african)
  expect_equal(mult$EUR, two$european)
})

test_that("read_bed_file finds .bim/.fam next to the .bed by default", {
  ext <- system.file("extdata", package = "lantern")
  bed <- file.path(ext, "toy_ancestry.bed")
  expect_identical(read_bed_file(bed),
                   read_bed_file(bed, sub("bed$", "bim", bed), sub("bed$", "fam", bed)))
})

test_that("write_dosage_gds validates its inputs", {
  m  <- matrix(1, 2, 2)
  vi <- data.frame(chrom = "19", pos = 1:2, ref = "A", alt = "T")
  expect_error(write_dosage_gds(m, vi[, 1:2], c("a", "b"), tempfile()),
               "`variant_info` must be a data.frame with columns")
  expect_error(write_dosage_gds(m, vi, 1:2, tempfile()),
               "`sample_ids` must be a character vector")
  expect_error(write_dosage_gds(m, vi, c("a", "b"), file.path(tempfile(), "x.gds")),
               "Directory for `gds_path` does not exist")
  expect_error(write_ancestry_gds(list(AFR = m), tempfile()),
               "`split_result` must be a list as returned by ancestry_split")
})

test_that("ancestry_smmat validates its inputs before fitting anything", {
  skip_if_not_installed("GMMAT")
  skip_if_not_installed("SeqVarTools")
  ext <- system.file("extdata", package = "lantern")
  ph  <- read.delim(file.path(ext, "toy_pheno.tsv"))
  K   <- diag(nrow(ph)); dimnames(K) <- list(ph$id, ph$id)
  gg  <- file.path(ext, "toy_genes.tsv")
  gds <- c(AFR = file.path(ext, "toy.vcf.gz"))   # exists; never opened
  expect_error(ancestry_smmat(c(AFR = "no.gds"), ph, y ~ age, K, gg),
               "`gds_paths` does not exist")
  expect_error(ancestry_smmat(gds, as.list(ph), y ~ age, K, gg),
               "`pheno` must be a data.frame")
  expect_error(ancestry_smmat(gds, ph, "y ~ age", K, gg), "`formula` must be")
  expect_error(ancestry_smmat(gds, ph, y ~ age, K, gg, family = "gaussian"),
               "`family` must be a family object")
  expect_error(ancestry_smmat(gds, ph, y ~ age, K, gg, ancestry_counts = matrix(1)),
               "together")
  expect_error(ancestry_smmat(gds, ph, y ~ age, K, gg, ncores = 0),
               "`ncores` must be a single positive whole number")
})

test_that("cauchy_combine rejects non-numeric input", {
  expect_error(cauchy_combine(c("0.1", "0.2")), "p_values must be numeric")
  expect_error(cauchy_combine(c(.1, .2), weights = c("1", "2")),
               "weights must be numeric")
})
