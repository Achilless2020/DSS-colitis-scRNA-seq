# Execution and output map

Paths below are relative to `config$output_dir`. Source figure comments are navigation hints, not proof that a generated panel is numerically identical to the submitted figure.

| Stage / source navigation | Inputs | Main outputs |
|---|---|---|
| import: original lines 69–86 | Six matrix/feature/barcode triplets, or three annotated RDS files | objects/raw.rds or objects/global.rds, fibro.rds, neu.rds; input_manifest.csv; input_scope.txt for checkpoint route |
| qc: original lines 89–177 | raw.rds | objects/qc.rds; qc/cell_counts_by_step.csv; gene lists; per-sample pK sweeps and doublet classifications/settings |
| global: original lines 182–244; source Fig1 labels | qc.rds, global mapping | global_clusters_unannotated.rds; seurat/global_cluster_markers_full.csv; annotated global.rds after review |
| subcluster: original fibro lines 496–519/1133 onward, neutrophil lines 1242–1516; source Fig2/3 labels | global.rds, subtype maps | fibro/neu unannotated and annotated objects; cluster markers/UMAPs; actual secondary QC thresholds and per-step counts |
| summary | Three annotated objects | full returned identity markers; metadata; per-sample counts/proportions including zeros; identity/group UMAPs; coarse_to_fine audit |
| monocle: original lines 1611–2015 | annotated neu.rds and its UMAP | graph_before_ordering.rds; separate cds_root_* objects; pseudotime_both_roots.csv; root candidates, selected roots and coverage; prespecified gene curves; optional graph_test_full.csv |
| nichenet manuscript: author-supplied simplified script plus original panel definitions | global/fibro/neu checkpoints, three network resources | nichenet/manuscript/*: four Figure 5 + six retained-original Figure 6 comparisons; native-wrapper agreement records and plots |
| nichenet figure5 / figure6 | Same inputs; explicit figure subset in config | nichenet/figure5/* or nichenet/figure6/*: the selected four or six comparisons; original references, sender pools and target directions |
| de: optional revision extension | raw counts and sample/group metadata | de/*/DE_full.csv; pseudobulk_counts.rds; design_matrix.csv; sample_coverage.csv; gene_filter.csv; analysis_status.csv |
| enrichment: optional revision extension | sample-level DE, supplied GMT, mouse OrgDb | enrichment/*/GSEA_full.csv; GO_ORA_up/down_full.csv; ranked genes, input gene lists and GMT hash |

For each NicheNet comparison the public script produces:

- `receiver_DE_full.csv`: all rows returned by the specified Seurat test/prefilter, including P, adjusted P, log2FC and detection proportions.
- `target_genes.csv`, `background_genes.csv`, `evaluation_universe.csv`: the exact gene sets and target/background membership.
- `candidate_ligands.csv`, `candidate_ligand_receptor_edges.csv`: expressed candidate ligands and receptor edges.
- `ligand_activities_full.csv`, `prioritized_ligands.csv`: all candidate activity scores and explicitly selected top-N/ties.
- `all_candidate_target_weights.rds`: modeled target-by-candidate weights for every candidate, without the visualization cutoff.
- `top_target_links_before_visualization.csv`, `ligand_target_visualization_matrix.csv`, `ligand_target_heatmap.pdf`: selected-link and visualization outputs, when visualization succeeds.
- `cells_used.csv`, `cells_by_sample.csv`, `parameters.json`, `complete_numeric_result.rds`: sample coverage, assumptions and the complete numeric comparison object.
- `original_wrapper_equivalence.csv`, `original_wrapper_result.rds`, `original_ligand_activity_target_heatmap.pdf`, `original_ligand_receptor_heatmap.pdf`, `original_ligand_expression_dotplot.pdf`: original-wrapper agreement and outputs when the default wrapper crosscheck is enabled.
- `neutrophil_label_alias_audit.csv` at the comparison-set level: original historical labels and the corresponding canonical analysis labels, without modifying the input objects.

Visualization errors are retained in `analysis_summary.csv` and the numeric result rather than causing valid scores to disappear. Empty or degenerate gene/candidate sets produce explicit failures. `last_failure.txt` and the failing process exit code identify a stopped stage. Per-stage completion files only certify that stage ran, not that all manuscript analyses have been validated.
