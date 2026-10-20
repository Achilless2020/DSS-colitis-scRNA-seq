# Paths are relative to this package; use config.local.R for private paths.
config <- list(
  seed = 123L,
  library_paths = character(),
  output_dir = "outputs/analysis",
  input_mode = "raw", # "raw" or "checkpoints"; checkpoints preserve saved clustering/QC
  samples = "metadata/samples.csv",
  objects = c(global = "inputs/global.rds", fibro = "inputs/fibro.rds", neu = "inputs/neu.rds"),
  feature_column = 1L, # Original used V1. Set 2 ONLY after inspecting features.tsv.gz.
  min_cells_per_gene = 3L, # keep genes detected in >3 cells, as in the source
  qc = list(min_features = 300, max_features = 5000, max_mito = 5, max_hb = 1,
            mito_pattern = "^mt-", ribo_pattern = "^Rp[sl]", hb_pattern = "^Hb[^(p)]"),
  doublets = list(enabled = TRUE, rate_mode = "source_heuristic", # or "sample_manifest"
                  pN = 0.25, dims = 1:15),
  global = list(npcs = 30L, neighbor_dims = 1:15, umap_dims = 1:30,
                resolutions = c(0.2, 0.5), final_resolution = 0.5,
                map = "metadata/global_clusters.csv", annotation_reviewed = FALSE),
  fibro = list(npcs = 50L, neighbor_dims = 1:15, umap_dims = 1:15,
               resolutions = c(0.2, 0.5, 1), final_resolution = 0.5,
               map = "metadata/fibro_clusters.csv", annotation_reviewed = FALSE),
  neu = list(npcs = 30L, neighbor_dims = 1:15, umap_dims = 1:15,
             resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1), final_resolution = 0.3,
             map = "metadata/neu_clusters.csv", annotation_reviewed = FALSE),
  markers = list(min_pct = 0.25, logfc = 0.25),
  monocle = list(num_dim = 50L, use_partition = FALSE, close_loop = FALSE,
                 minimal_branch_len = 5L, cores = 1L, graph_test = FALSE,
                 roots = c("Immature_Neu", "Immuno_silent_Neu")),
  nichenet = list(network_dir = "inputs/networks", verify_hashes = TRUE,
                 # Retained manuscript: four Figure 5 + six original Figure 6 comparisons.
                 # Use "figure5" or "figure6" for one figure only.
                 contrast_set = "manuscript", expression_pct = 0.1,
                 # Run the supplied script's wrappers too, verify agreement and export their plots.
                 wrapper_crosscheck = TRUE,
                 target_padj = 0.05, marker_prefilter_lfc = 0.1,
                 top_targets = 200L, visualization_quantile = 0.33),
  de = list(min_cells_per_sample = 10L, min_replicates_per_group = 2L),
  enrichment = list(gmt = "inputs/m5.go.v2024.1.Mm.symbols.gmt", min_size = 10L, max_size = 500L)
)
