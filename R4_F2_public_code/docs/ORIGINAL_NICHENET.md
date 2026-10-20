# Original NicheNet source and retained manuscript definitions

Updated 2026-10-20. Source filename: `nichenet简化版.R`; SHA-256 `49c17519e6a187d6573b083760e6e8c6e4259ff31eb96cb7ed148773624fc5be`. The original private script is not modified or copied with personal paths into this public package.

The supplied script contains four global calls and two executed subtype examples (IFN-response versus ISG-I-like; Slpi+ iFs versus vFs), with alternative receiver names in comments. The original six panel definitions and explicit author confirmation establish the complete Figure 6 mapping in `metadata/figure6_contrasts.csv`. Expanding the same call pattern to all six is a documented reconstruction of those definitions, not recovery of previously missing numerical results.

The default `manuscript` set executes four Figure 5 plus six Figure 6 comparisons. `figure5` and `figure6` are the only subset selections. The earlier `revision_figure6` redesign and incomplete `recorded` set have been retired; old local configurations using those names fail with an actionable message.

| Parameter | Global Figure 5 | Retained Figure 6 |
|---|---|---|
| Function | nichenet_seuratobj_aggregate | nichenet_seuratobj_cluster_de |
| Targets | upregulated, log2FC >=0.4 | both directions, abs(log2FC) >=0.25; geneset=DE |
| Receiver-gene adjusted P | <=0.05, Seurat Bonferroni | same |
| Expression prevalence | >=0.10 in each identity | same |
| Background/receptor expression | receiver pooled across stages | union of affected and reference receiver expression sets |
| Sender pool | other eight major identities | exactly three specified states in each direction |
| Ligand ranking | aupr_corrected; top30 and ties | aupr_corrected; top15 and ties |
| Prior target extraction | top200 per ligand intersect receiver targets | same |
| Visualization cutoff | weight quantile 0.33 | same |
| Ligand P/FDR cutoff | none | none |

The execution definitions preserve the fixed IS-I-like reference for neutrophils and fixed vFs reference for iFs. Figure 6A's receiver is IS-I inflammatory, as confirmed by the author. The old `ISG_*` labels are aliases only; source labels are exported alongside canonical analysis labels. No biological identities are reassigned by this naming operation. Barcode projection uses the saved fibroblast `celltype` column rather than the unassigned `Fib_subtype` field in the private script.

`R/04_nichenet.R` contains both a sparse audit calculation and `original_nichenet_wrapper`, which directly calls the same NicheNet wrappers as the supplied script. With the default crosscheck enabled, target/background sets, top ligands and AUROC/AUPR/corrected-AUPR/Pearson/rank values must agree (numeric tolerance 1e-10). The actual native combined activity/target heatmap is saved explicitly; saving never depends on whichever plot happened to be drawn last. Native plot files use the source's 60 x 20 cm export dimensions. Receptor and sender-expression plots and the full wrapper result are retained separately.

A run on current saved counts is a fresh execution of the documented original definitions. It is not proof of identical historical Figure 6 rankings or values. Original figures are retained in the manuscript workflow; no source plot, Word file or original RData is overwritten by this package.
