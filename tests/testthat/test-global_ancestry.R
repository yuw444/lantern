# ============================================================================
# GLA (global local ancestry) shrinkage: .assign_arm, .compute_arm_gla, and
# the w-weighted blend in split_diploid()/split_diploid_multi()
# ============================================================================

test_that(".assign_arm classifies p/q arms using the hg38 centromere table", {
  # chr19 centromere (inst/scripts/centromeres_hg38.R): p_end=24498980, q_start=27190874
  expect_equal(.assign_arm("chr19", 100), 0L)
  expect_equal(.assign_arm("19", 24498979), 0L)
  expect_equal(.assign_arm("chr19", 24498980), 1L)   # inside gap -> q tie-break
  expect_equal(.assign_arm("chr19", 30000000), 1L)
  expect_equal(.assign_arm(c("chr19", "chr19"), c(100, 30000000)), c(0L, 1L))
})

test_that(".assign_arm errors on non-autosomal chromosomes", {
  expect_error(.assign_arm("chrX", 100), "autosomal")
  expect_error(.assign_arm("chr23", 100), "autosomal")
})

test_that(".compute_arm_gla computes length-weighted per-arm proportions", {
  # 2 samples x 2 tracts, both on chr19's p-arm (well below its ~24.5Mb centromere)
  # diploid codes: 1=AFR, 2=EUR, 3=mixed(AFR/EUR)
  tract_df <- data.frame(chrom = c("chr19", "chr19"),
                          spos  = c(0, 1000000),
                          epos  = c(1000000, 5000000))
  anc_diploid_tract <- matrix(c(1, 3,
                                 3, 2),
                               nrow = 2, ncol = 2, byrow = TRUE)
  pure_codes_named <- c(AFR = 1L, EUR = 2L)
  mixed_codes_df <- data.frame(code = 3L, pop1 = "AFR", pop2 = "EUR",
                                stringsAsFactors = FALSE)

  gla <- .compute_arm_gla(tract_df, anc_diploid_tract, pure_codes_named, mixed_codes_df)

  expect_named(gla, "chr19")
  mat <- gla$chr19
  expect_equal(rownames(mat), c("p", "q"))
  expect_equal(colnames(mat), c("AFR", "EUR"))
  # T_AFR=1e6, T_EUR=4e6, T_mixed=5e6, total=1e7
  # GLA_AFR = (2*1e6+5e6)/2e7 = 0.35 ; GLA_EUR = (2*4e6+5e6)/2e7 = 0.65
  expect_equal(unname(mat["p", "AFR"]), 0.35, tolerance = 1e-8)
  expect_equal(unname(mat["p", "EUR"]), 0.65, tolerance = 1e-8)
  # No q-arm tracts at all -> falls back to whole-chromosome GLA, which here
  # equals the p-arm-only computation exactly (all tracts are on the p-arm)
  expect_equal(mat["q", ], mat["p", ])
})

test_that(".compute_arm_gla computes p/q arms independently when tracts span both", {
  tract_df <- data.frame(chrom = c("chr19", "chr19"),
                          spos  = c(0, 26000000),          # p-arm, q-arm
                          epos  = c(1000000, 27000000))
  anc_diploid_tract <- matrix(c(1, 2),   # sample1: p-arm=AFR, q-arm=EUR
                               nrow = 1, ncol = 2, byrow = TRUE)
  pure_codes_named <- c(AFR = 1L, EUR = 2L)
  mixed_codes_df <- data.frame(code = 3L, pop1 = "AFR", pop2 = "EUR",
                                stringsAsFactors = FALSE)

  gla <- .compute_arm_gla(tract_df, anc_diploid_tract, pure_codes_named, mixed_codes_df)
  mat <- gla$chr19
  expect_equal(unname(mat["p", "AFR"]), 1.0)
  expect_equal(unname(mat["p", "EUR"]), 0.0)
  expect_equal(unname(mat["q", "AFR"]), 0.0)
  expect_equal(unname(mat["q", "EUR"]), 1.0)
})

# ============================================================================
# split_diploid() (2-population) GLA blending
# ============================================================================

test_that("split_diploid with gla=NULL reproduces the pre-shrinkage estimator exactly", {
  gt <- matrix(c(1, 1, 1), nrow = 1, ncol = 3)
  an <- matrix(c(3, 2, 1), nrow = 1, ncol = 3)
  result <- split_diploid(gt, an)   # no gla/arm_id
  expect_equal(result$african[1, 2], 0.5)
  expect_equal(result$european[1, 2], 0.5)
})

test_that("split_diploid singleton case blends fully to GLA (w = 1)", {
  gt <- matrix(1, nrow = 1, ncol = 1)
  an <- matrix(2, nrow = 1, ncol = 1)   # single mixed het, no other evidence
  gla <- matrix(c(0.7, 0.3), nrow = 1, ncol = 2, dimnames = list("p", c("AFR", "EUR")))
  result <- split_diploid(gt, an, gla = gla, arm_id = 0L)
  expect_equal(result$african[1, 1], 0.7, tolerance = 1e-8)
  expect_equal(result$european[1, 1], 0.3, tolerance = 1e-8)
})

test_that("split_diploid with abundant unambiguous evidence ignores GLA (w = 0)", {
  # No mixed hets (N5=0) at all -> w=0, raw ratio unaffected by gla
  gt <- matrix(c(1, 1, 1), nrow = 1, ncol = 3)
  an <- matrix(c(3, 3, 1), nrow = 1, ncol = 3)   # 2 pure-AFR het, 1 pure-EUR het
  gla <- matrix(c(0.1, 0.9), nrow = 1, ncol = 2, dimnames = list("p", c("AFR", "EUR")))
  result_shrunk   <- split_diploid(gt, an, gla = gla, arm_id = c(0L, 0L, 0L))
  result_unshrunk <- split_diploid(gt, an)
  expect_equal(result_shrunk$african, result_unshrunk$african)
  expect_equal(result_shrunk$european, result_unshrunk$european)
})

test_that("split_diploid partial-w blends raw ratio and GLA proportionally", {
  # variant: 1 pure-AFR het (N2=1), 1 mixed het (N5=1); no hom-alt carriers
  # here, so allele-weighted total_alt happens to equal the unweighted
  # carrier count (2). w = N5/total_alt = 1/2
  # raw: denom = total_alt - N5 = (1+1) - 1 = 1 ; p1_raw = (0+1+0)/1 = 1, p2_raw = 0
  gt <- matrix(c(1, 1), nrow = 1, ncol = 2)
  an <- matrix(c(3, 2), nrow = 1, ncol = 2)
  gla <- matrix(c(0.2, 0.8), nrow = 1, ncol = 2, dimnames = list("p", c("AFR", "EUR")))
  result <- split_diploid(gt, an, gla = gla, arm_id = c(0L, 0L))
  w <- 0.5
  expect_equal(result$african[1, 2], (1 - w) * 1.0 + w * 0.2, tolerance = 1e-8)
  expect_equal(result$european[1, 2], (1 - w) * 0.0 + w * 0.8, tolerance = 1e-8)
})

test_that("split_diploid partial-w gives hom-alt carriers 2x the confidence weight of het carriers", {
  # variant: 1 pure-AFR HOM-ALT (N1=1, worth 2 alt alleles of confidence),
  # 1 mixed het (N5=1, ambiguous, worth 1 allele of confidence)
  # total_alt = 2*N1 + N5 = 2 + 1 = 3 ; w = N5/total_alt = 1/3
  # (a naive unweighted w = N5/carriers would give 1/2, which is wrong: the
  # hom-alt carrier's evidence should outweigh the ambiguous het's 2-to-1)
  gt <- matrix(c(2, 1), nrow = 1, ncol = 2)
  an <- matrix(c(3, 2), nrow = 1, ncol = 2)
  gla <- matrix(c(0.2, 0.8), nrow = 1, ncol = 2, dimnames = list("p", c("AFR", "EUR")))
  result <- split_diploid(gt, an, gla = gla, arm_id = c(0L, 0L))
  w <- 1 / 3
  expect_equal(result$african[1, 2], (1 - w) * 1.0 + w * 0.2, tolerance = 1e-8)
  expect_equal(result$european[1, 2], (1 - w) * 0.0 + w * 0.8, tolerance = 1e-8)
})

# ============================================================================
# split_diploid_multi() (K-population) GLA blending
# ============================================================================

test_that("split_diploid_multi with gla=NULL reproduces original behavior (K=3)", {
  gt  <- matrix(c(1, 1, 1, 1), nrow = 1, ncol = 4)
  anc <- matrix(c(1, 2, 4, 3), nrow = 1, ncol = 4)  # AFR pure, AFR/EUR mixed, NAT pure, AFR/NAT mixed
  pure  <- c(AFR = 1L, EUR = 2L, NAT = 3L)
  mixed <- data.frame(code = c(2L, 4L, 5L),
                      pop1 = c("AFR", "AFR", "EUR"),
                      pop2 = c("EUR", "NAT", "NAT"))
  out <- split_diploid_multi(gt, anc, pure, mixed)
  expect_named(out, c("AFR", "EUR", "NAT"))
})

test_that("split_diploid_multi singleton mixed pair blends fully to conditional GLA (w=1)", {
  # Single AFR/EUR mixed het, nothing else -> pk unavailable for this pair,
  # w=1, result should be gla's conditional (AFR/EUR)-pair fraction
  gt  <- matrix(1, nrow = 1, ncol = 1)
  anc <- matrix(4, nrow = 1, ncol = 1)   # mixed code 4 = AFR/EUR pair below
  pure  <- c(AFR = 1L, EUR = 2L, NAT = 3L)
  mixed <- data.frame(code = c(4L, 5L, 6L),
                      pop1 = c("AFR", "AFR", "EUR"),
                      pop2 = c("EUR", "NAT", "NAT"))
  gla <- matrix(c(0.5, 0.3, 0.2), nrow = 1, ncol = 3,
                dimnames = list("p", c("AFR", "EUR", "NAT")))
  out <- split_diploid_multi(gt, anc, pure, mixed, gla = gla, arm_id = 0L)
  # conditional fraction within the AFR/EUR pair: 0.5/(0.5+0.3) = 0.625
  expect_equal(out$AFR[1, 1], 0.625, tolerance = 1e-8)
  expect_equal(out$EUR[1, 1], 0.375, tolerance = 1e-8)
  expect_equal(out$NAT[1, 1], 0.0)
})

test_that("split_diploid_multi validates gla column names against pure_codes", {
  gt  <- matrix(1, nrow = 1, ncol = 1)
  anc <- matrix(4, nrow = 1, ncol = 1)
  pure  <- c(AFR = 1L, EUR = 2L, NAT = 3L)
  mixed <- data.frame(code = c(4L, 5L, 6L),
                      pop1 = c("AFR", "AFR", "EUR"),
                      pop2 = c("EUR", "NAT", "NAT"))
  bad_gla <- matrix(c(0.5, 0.3, 0.2), nrow = 1, ncol = 3,
                     dimnames = list("p", c("EUR", "AFR", "NAT")))  # wrong order
  expect_error(split_diploid_multi(gt, anc, pure, mixed, gla = bad_gla, arm_id = 0L),
               "gla's columns must match")
})

# ============================================================================
# End-to-end: ancestry_split() shrinks toward 1/2 (not GLA) in dosage mode
# ============================================================================

test_that("ancestry_split dosage mode shrinks toward 1/2 proportionally to w", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # 3 samples, chr19 p-arm tract: S1 pure AFR, S2 mixed AFR/EUR, S3 pure AFR
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1\tS3.0\tS3.1"
  ), con)
  writeLines(c(
    "chr19\t100\t500\t0.1\t0.2\t3\t0\t0\t0\t1\t0\t0"
  ), con)
  close(con)

  # Variant: S1 het (pure AFR carrier), S2 het (mixed, ambiguous), S3 hom-ref
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr19\t200\t.\tG\tC\t.\tPASS\t.\tGT\t0|1\t0|1\t0|0"
  ), vcf)

  result <- ancestry_split(vcf, msp, mode = "dosage", chrom = "chr19", verbose = FALSE)

  # Shrinkage target is 1/2 regardless of the tract's ancestry mix (this
  # tract's GLA would be 5/6 AFR -- the old target -- and must NOT be used).
  # w = N5/total_alt: N2=1 (S1 pure AFR het), N5=1 (S2 mixed het), no hom-alt
  # carriers here so total_alt equals the unweighted carrier count (2) -> w=1/2
  # raw p1 = (0+1+0)/((1+1)-1) = 1
  w <- 0.5
  expected_afr <- (1 - w) * 1.0 + w * 0.5      # = 0.75 (GLA would give 11/12)
  s2_idx <- match("S2", result$sample_ids)
  expect_equal(result$AFR[1, s2_idx], expected_afr, tolerance = 1e-8)
  expect_equal(result$AFR[1, s2_idx] + result$EUR[1, s2_idx], 1, tolerance = 1e-8)
})
