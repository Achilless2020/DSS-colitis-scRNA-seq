# Verification record

Verification date: 2026-10-20 (original-NicheNet replacement). Environment: R 4.5.1, Windows x64, with the installed packages and source commits recorded under `environment/`.

## Checks actually completed

- All R files (entry point, configuration, modules, tests, exporter and environment recorder) parsed successfully. The replacement adds `test_original_nichenet.R` to the preceding suite.
- `tests/test_contracts.R` passed: known SHA-256 value, missing/overlapping/foreign barcode rejection, complete subtype projection, preservation of coarse labels, unmapped cluster rejection, state-restricted root selection, sparse identity-specific prevalence, current manuscript/figure subset definitions, and exclusion of stale skipped DE files from enrichment.
- `tests/test_original_nichenet.R` passed: all ten manuscript comparisons, all six original A-F definitions, fixed references, exactly three specified senders, two-sided DE targets, receptor union, retired-design rejection and historical label aliases. The machine-readable panel CSV agrees with the executable definitions.
- `tests/smoke_pipeline.R` completed with exit code 0. It generated synthetic data and ran Seurat normalization/PCA/neighbors/clustering/UMAP, checkpoint import through the CLI, identity-marker summaries, two independent Monocle3 root orderings, all four original Figure 5 and six retained-original Figure 6 NicheNet contrasts, and 19 sample-level edgeR analyses (9 global identities, 4 fibroblast states and 6 neutrophil states).
- For all ten NicheNet contrasts, the actual original `nichenet_seuratobj_aggregate` / `nichenet_seuratobj_cluster_de` call also completed. Its target set, expression background and prioritized ligands matched the sparse audit calculation; AUROC, AUPR, corrected AUPR, Pearson and rank values agreed within 1e-10. Native combined activity/target plots, ligand/receptor plots and sender dotplots were exported. This uses the weighted-network resource, rather than hashing it without running the original wrapper.
- The same smoke test wrote six compressed matrix/feature/barcode input triplets and ran raw import/QC through the CLI. All six deliberately introduced high-mitochondrial cells were removed, and all retained cells met mitochondrial/hemoglobin thresholds. Two other synthetic cells also exceeded the hemoglobin cutoff and were correctly removed (1,020 imported, 1,012 retained). DoubletFinder was explicitly disabled for this synthetic integration run.
- The archive test passed: the ZIP includes `run.R` and its configuration and excludes private configuration, expression/network binaries and runtime outputs. The builder also checks ZIP integrity, per-file SHA-256 values and personal absolute-path strings.
- The preceding package's extracted entry point passed its dependency check. The replacement archive is independently integrity-checked and its own extracted entry point is checked before delivery; the local delivery manifest records the result.
- Independent code review identified and prompted fixes for an omitted archive entry point and potential reuse of stale DE results. Both have regression checks.

## Checks against the study's designated metadata

The author's already-extracted, designated saved-object metadata were verified separately. The global object contains 94,205 cells, the fibroblast-enriched subobject 1,031 cells, and the neutrophil subobject 12,562 cells. Sample-to-group assignments matched the six-sample manifest. Neutrophil labels included all six states, including 210 IS_I_Inflam labels. Barcode-based fine-label projection contained no missing labels, and detected the known 263 cells labeled coarse Fibroblasts but epithelial/glial in the subtype object. The original source SHA-256 matched the author-designated manifest and remained unchanged.

This check used extracted metadata, not the complete expression workspace, and does not independently verify the biological annotations. The actual local metadata verification JSON and complete synthetic execution log are retained outside the public code-only archive.

## What these checks do not establish

- There was no full study-data rerun from the original count matrices, and no claim of numerical reproduction of submitted figures. The raw route includes documented changes to the source workflow.
- The per-sample DoubletFinder sweep/classification was not exercised in the synthetic test. Its full-scale run, expected-rate justification and effects on downstream cells remain to be validated on the actual libraries.
- Raw-route biological cluster assignments still require marker review. Synthetic Seurat execution does not validate transfer of the original numeric cluster mappings.
- Optional `graph_test` and the GO/GSEA enrichment stage were parsed but not executed in this verification. The sample-level edgeR stage was executed on synthetic data only.
- The test used fresh `Rscript --vanilla` sessions with an existing installed package library. It was not a fresh dependency installation, container restore or an `renv.lock` restoration.
- No trajectory direction, cell identity, NicheNet mechanism, positive feedback, clinical generalization or biological significance is established by passing a software test.
- Public deposition, a code license, permanent version/DOI, shareable study-input access, final annotation reconciliation and a full study-data rerun remain author release steps. No external upload occurred.

The Monocle3 test emitted upstream ggplot2 deprecation warnings, and NicheNet emitted an upstream import warning. They did not prevent the tested stages from completing; the versions and the warning context are preserved in the local execution record.

The wrapper's sender DotPlots also warned about scaling a small number of synthetic groups. This describes the artificial test fixture, not a new study-data analysis. Source-script provenance is recorded in `ORIGINAL_NICHENET.md`; neither source R script, original RData, historical figure nor manuscript Word file was modified.
