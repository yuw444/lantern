# lantern 0.99.0

* Initial Bioconductor submission.
* `ancestry_split()` splits a phased VCF into per-ancestry dosages using
  RFMix local ancestry, in `"haplotype"` (phased) or `"dosage"` (unphased)
  mode, for any number of populations. Dosage mode shrinks ambiguous
  mixed-ancestry splits toward 1/2, in proportion to the share of the
  variant's evidence that is ambiguous. (An earlier development version
  shrank toward each chromosome arm's global local ancestry; that target and
  the `use_gla` argument were removed.)
* `write_ancestry_gds()` / `write_dosage_gds()` write SeqArray GDS files.
* `ancestry_smmat()` runs ancestry-stratified GMMAT SMMAT gene tests and
  Cauchy-combines the per-ancestry p-values (`cauchy_combine()`).
* VCFs are read with `bcftools` when it is on `PATH`, otherwise with
  SeqArray. Choose one explicitly with `options(lantern.vcf_reader = ...)`.
* Example dataset in `inst/extdata` (simulated; see `inst/scripts/make_toy_data.R`).
