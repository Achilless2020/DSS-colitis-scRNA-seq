stage_monocle <- function(cfg, root, out, samples) {
  require_packages(c("monocle3", "SingleCellExperiment", "SummarizedExperiment"))
  x <- join_rna(read_object(out, "neu"))
  validate_object(x, samples, "Neu_subtype")
  counts <- SeuratObject::GetAssayData(x, assay = "RNA", layer = "counts")
  md <- x[[]][colnames(counts), c("sample", "group", "Neu_subtype"), drop = FALSE]
  genes <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- monocle3::new_cell_data_set(counts, cell_metadata = md, gene_metadata = genes)
  nd <- min(cfg$monocle$num_dim, nrow(cds) - 1L, ncol(cds) - 1L)
  cds <- monocle3::preprocess_cds(cds, num_dim = nd)
  umap <- SeuratObject::Embeddings(x, reduction = "umap")[colnames(cds), 1:2, drop = FALSE]
  if (anyNA(umap) || any(!is.finite(umap))) stop("Invalid Seurat UMAP")
  SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- umap
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP", random_seed = cfg$seed)
  cds <- monocle3::learn_graph(cds, use_partition = cfg$monocle$use_partition,
            close_loop = cfg$monocle$close_loop,
            learn_graph_control = list(minimal_branch_len = cfg$monocle$minimal_branch_len))
  dest <- file.path(out, "monocle")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  saveRDS(cds, file.path(dest, "graph_before_ordering.rds"))
  save_plot(monocle3::plot_cells(cds, color_cells_by = "Neu_subtype", label_cell_groups = FALSE), file.path(dest, "trajectory_by_state.pdf"))
  save_plot(monocle3::plot_cells(cds, color_cells_by = "group", label_cell_groups = FALSE), file.path(dest, "trajectory_by_group.pdf"))
  # Each root uses the SAME graph and its OWN output. Neither is validated chronology.
  pseudotimes <- data.frame(cell_id = colnames(cds), md, check.names = FALSE)
  roots <- list()
  for (state in cfg$monocle$roots) {
    candidates <- colnames(cds)[as.character(md$Neu_subtype) == state]
    id <- root_near_centroid(umap, candidates)
    ordered <- monocle3::order_cells(cds, root_cells = id)
    pt <- monocle3::pseudotime(ordered)[colnames(cds)]
    if (all(!is.finite(pt))) stop("No finite pseudotimes for root ", state)
    pseudotimes[[paste0("pseudotime_root_", state)]] <- unname(pt)
    roots[[state]] <- data.frame(root_state = state, cell_id = id, candidates = length(candidates),
                                 finite_cells = sum(is.finite(pt)), unreachable_cells = sum(!is.finite(pt)))
    saveRDS(ordered, file.path(dest, paste0("cds_root_", state, ".rds")))
    save_plot(monocle3::plot_cells(ordered, color_cells_by = "pseudotime", label_cell_groups = FALSE,
               label_leaves = FALSE, label_branch_points = FALSE), file.path(dest, paste0("pseudotime_root_", state, ".pdf")))
    write_table(data.frame(cell_id = candidates), file.path(dest, paste0("root_candidates_", state, ".csv")))
    selected_genes <- intersect(c("Camp", "Ltf", "Ngp", "Il1b", "Tnf", "Cxcl2", "Ifit1", "Ifit3", "Isg15", "Rsad2", "Cd274"), rownames(ordered))
    if (length(selected_genes)) {
      finite <- colnames(ordered)[is.finite(pt)]
      plot <- monocle3::plot_genes_in_pseudotime(ordered[selected_genes, finite], color_cells_by = "Neu_subtype", min_expr = 0.5)
      save_plot(plot, file.path(dest, paste0("prespecified_genes_root_", state, ".pdf")), height = max(6, length(selected_genes) * 1.4))
    }
  }
  write_table(pseudotimes, file.path(dest, "pseudotime_both_roots.csv"))
  write_table(do.call(rbind, roots), file.path(dest, "root_selection_and_coverage.csv"))
  if (isTRUE(cfg$monocle$graph_test)) {
    gt <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = cfg$monocle$cores)
    write_table(data.frame(gene = rownames(gt), gt), file.path(dest, "graph_test_full.csv"))
  }
  jsonlite::write_json(cfg$monocle, file.path(dest, "parameters.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(c("Pseudotime is conditional on an assumed root and is not elapsed time or lineage tracing.",
               "use_partition=FALSE follows the source but does not validate connections between disconnected biological populations.",
               "Gene curves are prespecified visualizations. graph_test, if enabled, tests graph association, not causation or direction."), file.path(dest, "interpretation.txt"))
}
