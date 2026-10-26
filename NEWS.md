# lantern 0.99.1

Changes in response to the Bioconductor review.

* Bug fixes in the matrix input of `ancestry_split_dosage(gt_matrix,
  pt_matrix)`. The VCF/MSP input and `ancestry_split()` were not affected.
  * Genotypes were split with the wrong samples' ancestry calls, because the
    ancestry matrix was passed to the C code without transposing it.
  * Variants and regions matched by name are now paired by name; they were
    paired by position, so a reordered `pt_matrix` gave wrong results.
  * Coordinate matching (`"chr:pos"` within `"chr:start-end"`) now allows
    several variants per region, which previously gave an error.
  * `counts` are now per variant, as documented (they were per sample).
  * Unnamed matrices of different sizes no longer fail, and inputs with no
    variant overlap give a clear error.
* Exported functions now check their arguments and give informative errors.
* New defaults: `split_diploid_multi()` defaults to the two-population codes
  of `split_diploid()`; `read_bed_file()` finds the `.bim` and `.fam` files
  next to the `.bed` file; `ancestry_split_phased()` needs `out_path` only
  when `write_vcf = TRUE`.
* `ancestry_smmat()` now requires `ancestry_counts` and `variant_info`
  together (previously one without the other was silently ignored).
* The C code checks for user interrupts.
* Long functions were split into smaller helpers and loops vectorized;
  results are unchanged.
* Removed the unused global-ancestry helpers and the internal
  `centromeres_hg38` table.
* Vignettes rewritten for readability; the real-data vignette's source is
  no longer shipped in `vignettes/`. The Getting started vignette now
  explains the method's motivation with figures from the LANTERN preprint.
* Added `citation("lantern")` for the preprint (Wang et al. 2026, medRxiv).

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
