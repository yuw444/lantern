# Generates the small simulated example dataset shipped in inst/extdata/,
# used by the package examples and vignette. Run from the package root:
#
#   Rscript inst/scripts/make_toy_data.R
#
# 60 admixed samples, two populations (AFR=0, EUR=1), 8 RFMix ancestry
# windows on chr19 (4 per arm, hg38 centromere ~24.5-27.2 Mb), 25 rare
# variants per window. Each window is treated as one gene. The phenotype
# carries an AFR-specific burden effect in GENE2.
#
# Outputs (all in inst/extdata/):
#   toy.vcf.gz              phased VCF (GT only)
#   toy.msp.tsv             RFMix MSP local-ancestry file
#   toy_genes.tsv           SMMAT gene-group file (no header)
#   toy_pheno.tsv           phenotype + covariates (id, y, age, sex)
#   toy_ancestry.{bed,bim,fam}  per-window diploid ancestry codes as PLINK
#                           genotypes (1=EUR/EUR, 2=AFR/EUR, 3=AFR/AFR),
#                           for read_bed_file()

set.seed(2026)
out_dir <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_samples <- 60
ids       <- sprintf("S%03d", seq_len(n_samples))
spos      <- as.integer(c(1, 5, 10, 15, 30, 35, 40, 45) * 1e6)
epos      <- spos + 4000000L - 1L
n_win     <- length(spos)
n_per_win <- 25

# ---- Local ancestry: each haplotype AFR (0) w.p. 0.8 per window ----
anc_h0 <- matrix(rbinom(n_samples * n_win, 1, 0.2), n_samples, n_win)
anc_h1 <- matrix(rbinom(n_samples * n_win, 1, 0.2), n_samples, n_win)

hap_cols <- as.vector(rbind(paste0(ids, ".0"), paste0(ids, ".1")))
msp_rows <- vapply(seq_len(n_win), function(w) {
  calls <- as.vector(rbind(anc_h0[, w], anc_h1[, w]))
  paste(c("chr19", spos[w], epos[w], sprintf("%.2f", spos[w] / 1e6),
          sprintf("%.2f", epos[w] / 1e6), n_per_win, calls), collapse = "\t")
}, character(1))
writeLines(c("#Subpopulation order/codes: AFR=0\tEUR=1",
             paste(c("#chm", "spos", "epos", "sgpos", "egpos", "n snps",
                     hap_cols), collapse = "\t"),
             msp_rows),
           file.path(out_dir, "toy.msp.tsv"))

# ---- Variants: rare, typically more common on AFR haplotypes ----
win_of <- rep(seq_len(n_win), each = n_per_win)
pos    <- unlist(lapply(seq_len(n_win), function(w)
  sort(sample(spos[w]:epos[w], n_per_win))))
bases  <- c("A", "C", "G", "T")
ref    <- sample(bases, length(pos), replace = TRUE)
alt    <- vapply(ref, function(r) sample(setdiff(bases, r), 1), character(1))
f_afr  <- runif(length(pos), 0.005, 0.06)
f_eur  <- runif(length(pos), 0, 0.02) * rbinom(length(pos), 1, 0.6)

draw_allele <- function(anc, v) rbinom(length(anc), 1,
                                       ifelse(anc == 0, f_afr[v], f_eur[v]))
h0 <- t(vapply(seq_along(pos), function(v) draw_allele(anc_h0[, win_of[v]], v),
               numeric(n_samples)))
h1 <- t(vapply(seq_along(pos), function(v) draw_allele(anc_h1[, win_of[v]], v),
               numeric(n_samples)))

con <- gzfile(file.path(out_dir, "toy.vcf.gz"), "wt")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr19,length=58617616>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "FORMAT", ids), collapse = "\t"),
  vapply(seq_along(pos), function(v)
    paste(c("chr19", pos[v], ".", ref[v], alt[v], ".", "PASS", ".", "GT",
            paste0(h0[v, ], "|", h1[v, ])), collapse = "\t"),
    character(1))
), con)
close(con)

# ---- Gene groups: one gene per ancestry window ----
genes <- data.frame(gene = paste0("GENE", win_of), chr = "19", pos = pos,
                    ref = ref, alt = alt, weight = 1)
write.table(genes, file.path(out_dir, "toy_genes.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- Phenotype: AFR-specific burden effect in GENE2 ----
in_gene2   <- win_of == 2
afr_burden <- colSums(h0[in_gene2, ] * (anc_h0[, 2] == 0)[col(h0[in_gene2, ])] +
                      h1[in_gene2, ] * (anc_h1[, 2] == 0)[col(h1[in_gene2, ])])
age <- round(rnorm(n_samples, 50, 10))
sex <- rbinom(n_samples, 1, 0.5)
y   <- 0.02 * age + 0.3 * sex + 1.5 * afr_burden + rnorm(n_samples)
write.table(data.frame(id = ids, y = round(y, 4), age = age, sex = sex),
            file.path(out_dir, "toy_pheno.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# ---- PLINK BED of diploid ancestry codes (windows as "variants") ----
# Ancestry code (1=EUR/EUR, 2=AFR/EUR, 3=AFR/AFR) -> PLINK 2-bit code in
# snpStats numeric order: 1 = hom first (00), 2 = het (10), 3 = hom second (11).
anc_code <- 3L - ((anc_h0 == 1) + (anc_h1 == 1))   # samples x windows
bits     <- c(`1` = 0L, `2` = 2L, `3` = 3L)
bytes_per_var <- ceiling(n_samples / 4)
bed_body <- unlist(lapply(seq_len(n_win), function(w) {
  b <- bits[as.character(anc_code[, w])]
  b <- c(b, rep(0L, bytes_per_var * 4 - n_samples))
  vapply(seq_len(bytes_per_var), function(k) {
    q <- b[(4 * k - 3):(4 * k)]
    as.integer(q[1] + bitwShiftL(q[2], 2) + bitwShiftL(q[3], 4) +
               bitwShiftL(q[4], 6))
  }, integer(1))
}))
con <- file(file.path(out_dir, "toy_ancestry.bed"), "wb")
writeBin(as.raw(c(0x6c, 0x1b, 0x01, bed_body)), con)
close(con)
writeLines(paste("19", sprintf("chr19:%d-%d", spos, epos), 0, spos, "A", "G",
                 sep = "\t"),
           file.path(out_dir, "toy_ancestry.bim"))
writeLines(paste("FAM", ids, 0, 0, 0, -9, sep = "\t"),
           file.path(out_dir, "toy_ancestry.fam"))
