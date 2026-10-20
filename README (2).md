# DSS colitis single-cell analysis: R4-F2 review package

This package reorganizes the author-designated `scRNA_analysis_codenew.R` into an ordered, configurable analysis. It contains raw matrix import, explicit QC and per-sample DoubletFinder, global Seurat clustering, fibroblast and neutrophil subclustering, marker exports, Monocle3 trajectories under two alternative roots, and NicheNet analyses. Optional sample-level pseudobulk DE and enrichment are revision analyses. No study expression matrices or network binaries are bundled.

**Status: local code prepared for review; not publicly deposited.** Synthetic integration tests assess code execution, not reproduction of the study's figures. See `docs/PROVENANCE.md` and `docs/VALIDATION.md` before interpreting a rerun.

## Inputs and environment

Use R 4.5.1 and the versions/source commits in `environment/`. Copy `config.R` to `config.local.R`, then edit private paths and analysis settings. Paths are resolved relative to this package, so the commands work from another working directory as well.

`metadata/samples.csv` defines S1–S3 as Chronic and S4–S6 as Acute. These are six independent biological samples, not six technical replicates. No healthy control is inferred. Metadata must contain exactly one group per sample; the loader rejects inconsistent checkpoint labels.

For each raw sample, supply `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz`. Inspect the feature table before setting `feature_column`: the source uses column 1, but the public loader rejects predominantly Ensembl IDs because this workflow requires mouse symbols. It also rejects duplicate symbols rather than silently renaming or aggregating them. Document any required gene-ID conversion or duplicate aggregation upstream.

For NicheNet, provide the three mouse RDS resources listed in `metadata/network_hashes.csv` under the configured network directory. SHA-256 verification is on by default. Obtain files through the NicheNet authors' resource distribution. For optional GSEA, supply the named mouse MSigDB GMT under its terms; it is not redistributed here.

```sh
Rscript --vanilla run.R --config config.local.R --stages check
```

## Route A: start from raw matrices

```sh
Rscript --vanilla run.R --config config.local.R --stages import,qc,global
```

Global clustering saves unannotated objects, complete returned cluster-marker tables, cell-to-cluster mappings and UMAPs before applying labels. The distributed mapping CSVs record the original numeric cluster assignments. **They are not universal annotations.** Review the marker evidence from this rerun, update the relevant map if needed, and then set `config$global$annotation_reviewed <- TRUE` in the local config. Rerun `global` to save its annotated object.

```sh
Rscript --vanilla run.R --config config.local.R --stages global
Rscript --vanilla run.R --config config.local.R --stages subcluster
```

Subclustering follows the same review step, separately for fibroblasts and neutrophils (`annotation_reviewed` in each setting). When a population is awaiting review, its unannotated object, marker table and plots are retained. The stage stops with an explanatory message. After both mappings are reviewed, run:

```sh
Rscript --vanilla run.R --config config.local.R --stages subcluster,summary,monocle,nichenet
```

QC defaults use >300 and <5000 detected genes, genes detected in >3 cells, mitochondrial percentage <5, and hemoglobin percentage <1. Ribosomal percentage is recorded but not filtered, matching the executed source operations. The public mitochondrial pattern is `^mt-`, case-insensitive, replacing the source's broad case-sensitive `^Mt`; inspect the exported gene list against the actual reference. **This change and the repaired QC-to-DoubletFinder chain make a raw rerun a revised analysis, not a guaranteed reconstruction of saved study results.**

DoubletFinder runs on the QC-retained cells of each sample, with PCs 1:15, pN=0.25, a recorded pK sweep, and singlet retention. The default `source_heuristic` uses the notebook's `0.008 * n_cells/1000` expected rate. This is an unverified source assumption, not a platform-specific validated rate. To use measured library-specific rates, select `sample_manifest` and fill `expected_doublet_rate` (fraction, not percent). No homotypic adjustment is invented. Every sample's settings and classifications are exported. Disabling DoubletFinder is explicit and recorded.

Global neighbor PCs are 1:15 and UMAP PCs 1:30; resolution 0.5 is selected. Fibroblast neighbor/UMAP PCs are 1:15, resolution 0.5, with all-gene scaling. Neutrophil neighbor/UMAP PCs are 1:15, resolution 0.3. Candidate resolutions and PCA dimensions are in `config.R`. Neutrophil secondary QC uses >max(200, 0.25×median features), <2×median features, <3×median counts, mitochondrial percentage <10, and hemoglobin percentage <5. All boundaries are strict; the actual computed thresholds are exported. Both subpopulations are normalized once in their own analysis. No Harmony/CellChat analysis is claimed merely because those packages appeared in the source's imports.

## Route B: continue from saved annotated Seurat objects

Set `config$input_mode <- "checkpoints"` and provide global, fibroblast, and neutrophil RDS files. The expected source objects are `sce`, `sce.F`, and `sce.neu`. If these are held in a private RData:

```sh
Rscript --vanilla tools/export_workspace.R /private/author_workspace.RData /private/checkpoints
```

Set the three corresponding paths in `config$objects`, then run:

```sh
Rscript --vanilla run.R --config config.local.R --stages import,summary,monocle,nichenet
```

Required fields are `sample`, `group`, global/fibroblast `celltype`, and neutrophil `Neu_subtype`. RNA must contain counts and normalized data; saved UMAP embeddings are required. This route preserves the supplied cells, clustering and annotations; it does not claim that raw QC, doublet removal or subclustering have been rerun. SHA-256 input hashes identify the actual files.

## Monocle3

The neutrophil raw counts initialize a cell data set. Preprocessing is followed by Seurat UMAP transfer through the supported reduced-dimension accessor, Monocle3 clustering and graph learning (`use_partition=FALSE`, `close_loop=FALSE`, `minimal_branch_len=5`). No internal UMAP slots are edited.

The same learned graph is ordered independently from an `Immature_Neu` root and an `Immuno_silent_Neu` root. The selected cell is nearest its state's UMAP centroid **among cells of that state**. Both objects, both pseudotime columns, root candidate barcodes and the number of unreachable cells are saved. The old script searched all cells for the nearest root and overwrote the first ordering; those behaviors are explicitly changed. Pseudotime is conditional on the chosen root and does not independently establish developmental direction, disease chronology or lineage transitions.

The gene curves use prespecified genes; they are not discovered trajectory-associated genes. `config$monocle$graph_test=TRUE` optionally exports the full Monocle3 graph-association test. The pipeline does not call the Monocle2 `differentialGeneTest` function.

## NicheNet: choose the comparison set explicitly

`config$nichenet$contrast_set <- "recorded"` runs the four global comparisons actually specified in the source (Acute vs Chronic and Chronic vs Acute for each receiver lineage) and the two explicitly written subtype comparisons (IFN_response_Neu vs IS_I_like; Slpi+ iFs vs vFs). The global analyses use upregulated targets and log2FC ≥0.4, top30. The two original subtype comparisons preserve the installed wrapper's two-sided DE default, absolute log2FC ≥0.25, top15. This reconstructs explicit settings; it does not recover missing historical Figure 6 panels.

`config$nichenet$contrast_set <- "revision_figure6"` runs six **new revision comparisons**, matching the separately documented R4-M4 design: each of the three selected neutrophil receiver states vs the other two, and each of three selected fibroblast states vs the other two. Sender pools are all four specified fibroblast states or all six specified neutrophil states, respectively. Stages are pooled. Targets are upregulated (log2FC ≥0.25); receptors must be expressed in the affected receiver alone. The original six source-code comparisons and the new six-panel design are saved in distinct output folders. To produce both, edit the local setting and rerun only `nichenet`.

All analyses use RNA data, ≥10% nonzero expression within each identity, and target `p_val_adj ≤0.05`. Seurat `p_val_adj` is Bonferroni-adjusted over assay genes; it must not be described as BH. The explicitly recorded FindMarkers logFC prefilter is 0.1. Receiver-expression backgrounds are intersected with model targets. For subtype analyses the background is the union over compared receiver states; global expression filtering pools both conditions. Targets outside this expression background remain identified in the exported evaluation union, as in the recorded NicheNet implementation. The corrected AUPR uses the target prevalence over that union.

Ligands are ranked by corrected AUPR, with ties retained at top-N. **There is no ligand P/FDR significance cutoff.** The 0.33 setting is the 33rd percentile of target-link weights for visualization. Each comparison exports all returned DE results, targets, background, evaluation universe, candidate ligand–receptor edges, all ligand scores, selected ligand links, all modeled target weights for every candidate (RDS), cell/sample coverage, parameters and heatmaps. The weighted-network resource is hashed for provenance; its optional receptor/ligand heatmaps from the old wrapper are not generated here. Predictions do not demonstrate active signaling, secretion, disease causation or positive feedback.

## Optional sample-level DE and enrichment

```sh
Rscript --vanilla run.R --config config.local.R --stages de
Rscript --vanilla run.R --config config.local.R --stages enrichment
```

The `de` stage aggregates raw counts per biological sample within an identity, fits edgeR quasi-likelihood models with a Chronic-minus-Acute contrast, and exports all tested genes, BH FDR, sample coverage, gene filters, count matrices and design matrices. It requires the configured cell count and replicate coverage; unsupported identities are explicitly skipped. Fibroblast DE uses the four named fibroblast states and excludes the epithelial/glial labels in that enriched object. Global DE continues to use the supplied coarse global labels. This is a revision extension and not the source's pooled cell-level disease test. The descriptive per-sample proportions have no P values; proportion inference belongs to the separate R4-M6 analysis.

Enrichment uses the sample-level DE results, the provided mouse GMT for GSEA, and org.Mm.eg.db GO overrepresentation with the tested-gene universe. It requires the current DE stage to have completed, and reads only its completed identity entries, so stale results from identities skipped on a rerun cannot enter enrichment. Up/down ORA inputs require DE FDR <0.05. Outputs retain all returned pathway tests and their adjusted P values. GSEA ranks all tested genes by log2FC, as specified explicitly in code. These optional outputs have not been verified on the study data.

## Execution records and verification

Outputs include effective configuration, per-stage `config_*.R` records, input/network/GMT hashes, running/completed/failed stage status files, start/end sessions, QC step counts and `last_failure.txt` on failure. Use a new `output_dir` for each scientific analysis. Rerunning a stage replaces produced files; older files for skipped identities can remain on disk, but enrichment excludes them by consulting current DE status. A successful later stage is not proof that an earlier analysis used the same inputs/configuration: retain a complete run folder with its manifests.

```sh
Rscript --vanilla tests/test_contracts.R
Rscript --vanilla tests/smoke_pipeline.R
```

The smoke test creates synthetic matrices in a temporary directory and runs the command-line workflow. Its DoubletFinder step is deliberately disabled; raw-data doublet sweeps and the study-scale analysis need separate validation. See `docs/VALIDATION.md` for the actual verification status and unresolved annotation matters.

## Public deposition

`tools/build_archive.py` builds an allowlisted code-only ZIP and SHA-256 inventory. It excludes local configuration, inputs, runtime outputs and study data. The author still needs to choose a code license, approve the shareable input scope, deposit an immutable version in a public repository/archive, and insert its real URL/DOI in the manuscript. No repository URL, DOI, data accession or public availability statement is fabricated by this package.

Official API references: [Seurat v5 commands](https://satijalab.org/seurat/articles/seurat5_essential_commands), [Monocle3 trajectories](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/), and [the recorded NicheNet source commit](https://github.com/saeyslab/nichenetr/blob/2d5c1ab5e135ebe6634d56c5ee13b8207ac178a1/R/application_prediction.R).
