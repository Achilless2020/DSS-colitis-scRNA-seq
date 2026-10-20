# Source provenance and deliberate changes

The author-designated source is `scRNA_analysis_codenew.R`, SHA-256 `8fb51c68eab1a7e1b5dbfe3716045956d5e95a96a1dc6d9a3e668cc30c86636c`. It is retained unchanged in the author's working directory. Original desktop paths are not redistributed in the code archive. The current author-designated workspace SHA-256 recorded in the local input audit is `fae48fceafe818c84bd3f9a422cc2beb2df2c2302169881af825afa2eff019fc`; the workspace itself is not included.

| Source section | Public implementation | Interpretation |
|---|---|---|
| Matrix import, QC, DoubletFinder | `R/01_input_qc.R` | Repaired QC object chain; sample prefixing; explicit feature-column and gene checks; per-step counts and doublet parameters |
| Global Seurat preprocessing/annotation | `R/02_seurat.R`, global mapping CSV | Explicit final resolution 0.5; source neighbor/UMAP PCs retained; review required for numeric mapping after reclustering |
| Fibroblast subclustering and annotation | `R/02_seurat.R`, fibro mapping CSV | Resolution 0.5 selected explicitly; original celltype names preserved, including source spelling lY6C1; contaminants not silently deleted |
| Repeated neutrophil subclustering sections | `R/02_seurat.R`, neu mapping CSV | Consolidated into one execution; actual formula multipliers retained, rather than the conflicting prose comments; all six levels preserved |
| Neutrophil Monocle3 | `R/03_monocle.R` | Public accessors; independent root-specific outputs; restricted root-cell search; no unsupported direction claims |
| Four global and two explicit subtype NicheNet calls | `R/contrasts.R`, recorded set | Explicit original settings; saved-data reconstruction only if provided checkpoints are used |
| Separately documented R4-M4 revision design | `R/contrasts.R`, revision_figure6 set | Six new receiver/reference definitions and full sender pools; affected-only receptors; upregulated targets; stages pooled |
| Source pooled disease DE/GSEA | `R/05_de_enrichment.R` | Optional sample-level edgeR revision extension; not asserted to reproduce old cell-level DE or pathway results |

Changes that can alter scientific results are explicit:

1. The original script created `sceilt` for QC but split the unassigned `sce` for DoubletFinder. The raw route now sends the QC-retained cells into sample-specific doublet removal. It will therefore not necessarily reproduce saved cell counts, clusters, labels or plots.
2. The source's case-sensitive mitochondrial `^Mt` is replaced by case-insensitive `^mt-`, the common mouse mitochondrial symbol prefix. This avoids treating every `Mt...` symbol as mitochondrial. The actual reference still needs inspection; all matching genes are exported. The old pattern is documented rather than silently described as the new one.
3. Global and subtype annotation is obtained by a complete mapping and checked for missing values. The public code never truncates six neutrophil labels to the old five-state `progression_levels`. Fibroblast projection uses the real `sce.F$celltype` field, not the absent `Fib_subtype` field.
4. Fibroblast normalization is performed explicitly before variable-feature selection/PCA. In the source it initially inherited parent normalization and was normalized again after clustering. This public raw route is a consistent revised workflow, not a numerical identity claim.
5. Monocle preprocessing occurs before Seurat UMAP import; overwritten manual PCA assignments and direct internal-slot edits are removed. The root candidate pool is restricted to the named root state. Both orderings are preserved. These are deliberate implementation changes.
6. The two subtype calls in the source use `geneset="DE"` by the installed NicheNet wrapper default. The new six-panel revision uses `geneset="up"` and a different receptor-expression scope. These are separate named comparison sets, not interchangeable results.
7. Marker tests use Seurat Wilcoxon and Seurat's Bonferroni `p_val_adj`. Optional disease DE uses biological-sample pseudobulk and edgeR BH FDR. NicheNet scores themselves do not carry ligand P/FDR values.
8. Original `rm(list=ls())`, `setwd`, repeated library blocks, internal assay casts and unrelated unused imports are removed. Counts/data layers are joined through the Seurat v5 API. The original source files and saved objects are never overwritten.

## Identity questions intentionally left visible

The current saved global object labels 1,031 cells as Fibroblasts, while the revised fibroblast-enriched subobject contains 255 epithelial and 8 enteric glial labels (263 cells in total), with four named fibroblast states comprising 768 cells. The public code preserves the coarse labels, projects fine labels in a separate vector, and exports a coarse-to-fine audit. It does not claim to have reconciled final global proportions or validated all four fibroblast identities.

The Slpi clusters' immune/mesothelial signals, the former mFs assignment, neutrophil cluster 4's digestive-enzyme signals, and hemoglobin/mitochondrial-rich states remain biological annotation questions. Restoring a label or making a script run is not validation of the cell identity. The package must be rerun after any scientific change to cell membership or annotation.

No animal-level proportion significance test is inferred from a stacked bar chart. No NicheNet axis is presented as functionally demonstrated. No trajectory is presented as independently validated development or an Acute-to-Chronic transition.
