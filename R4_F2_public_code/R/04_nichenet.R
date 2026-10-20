original_nichenet_wrapper <- function(x, ct, ltm, lr, weighted, settings) {
  args <- list(seurat_obj = x, sender = ct$sender, ligand_target_matrix = ltm, lr_network = lr,
       weighted_networks = weighted, assay_oi = "RNA", expression_pct = settings$expression_pct,
       lfc_cutoff = ct$lfc, geneset = ct$geneset, filter_top_ligands = TRUE,
       top_n_ligands = ct$top_ligands, top_n_targets = settings$top_targets,
       cutoff_visualization = settings$visualization_quantile, verbose = FALSE)
  if (ct$mode == "global") {
    SeuratObject::Idents(x) <- "celltype"; args$seurat_obj <- x
    args$receiver <- ct$receiver; args$condition_colname <- "group"
    args$condition_oi <- ct$condition; args$condition_reference <- ct$reference
    do.call(nichenetr::nichenet_seuratobj_aggregate, args)
  } else {
    SeuratObject::Idents(x) <- "celltype_fine"; args$seurat_obj <- x
    args$receiver_affected <- ct$receiver; args$receiver_reference <- ct$reference
    do.call(nichenetr::nichenet_seuratobj_cluster_de, args)
  }
}
stage_nichenet <- function(cfg, root, out, samples) {
  require_packages(c("nichenetr", "dplyr", "tidyr"))
  files <- c("ligand_target_matrix_nsga2r_final_mouse.rds", "lr_network_mouse_21122021.rds", "weighted_networks_nsga2r_final_mouse.rds")
  paths <- file.path(path_from(root, cfg$nichenet$network_dir), files)
  dest <- file.path(out, "nichenet", cfg$nichenet$contrast_set)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  hashes <- hash_files(paths)
  expected <- utils::read.csv(file.path(root, "metadata", "network_hashes.csv"))
  if (isTRUE(cfg$nichenet$verify_hashes) && !identical(hashes$sha256, expected$sha256[match(hashes$file, expected$file)])) stop("NicheNet network hash mismatch; do not relabel a different network as the recorded model")
  write_table(hashes, file.path(dest, "network_manifest.csv"))
  ltm <- readRDS(paths[1]); lr <- readRDS(paths[2])
  if (!is.matrix(ltm) || !all(c("from", "to") %in% names(lr))) stop("Invalid network format")
  weighted <- if (isTRUE(cfg$nichenet$wrapper_crosscheck)) readRDS(paths[3]) else NULL
  x <- join_rna(read_object(out, "global"))
  validate_object(x, samples, "celltype")
  fib <- read_object(out, "fibro"); neu <- read_object(out, "neu")
  fine <- project_labels(x[[]], fib[[]], neu[[]])
  x$celltype_fine <- unname(fine[colnames(x)])
  write_table(data.frame(cell_id = colnames(neu), source_label = as.character(neu$Neu_subtype),
             analysis_label = canonical_neutrophil_labels(neu$Neu_subtype)), file.path(dest, "neutrophil_label_alias_audit.csv"))
  data <- SeuratObject::GetAssayData(x, assay = "RNA", layer = "data")
  if (!nrow(data) || !identical(colnames(data), colnames(x))) stop("Normalized RNA data missing/misaligned")
  global_ids <- checked_labels(x$celltype, "global")
  fine_ids <- checked_labels(x$celltype_fine, "fine")
  global_ex <- expressed_by_identity(data, global_ids, unique(global_ids), cfg$nichenet$expression_pct)
  fine_ex <- expressed_by_identity(data, fine_ids, unique(fine_ids), cfg$nichenet$expression_pct)
  rm(data); gc(verbose = FALSE)
  summary <- list()
  for (ct in make_contrasts(cfg$nichenet$contrast_set)) {
    message("NicheNet: ", ct$id)
    folder <- file.path(dest, ct$id)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    is_global <- ct$mode == "global"
    ids <- if (is_global) global_ids else fine_ids
    ex <- if (is_global) global_ex else fine_ex
    receiver_types <- if (is_global) ct$receiver else c(ct$receiver, ct$reference)
    required <- unique(c(receiver_types, ct$sender))
    if (!all(required %in% names(ex))) stop("Missing NicheNet identities: ", paste(setdiff(required, names(ex)), collapse = ", "))
    background <- intersect(unique(unlist(ex[receiver_types], use.names = FALSE)), rownames(ltm))
    receptor_types <- if (ct$affected_receptors_only) ct$receiver else receiver_types
    receptors <- unique(unlist(ex[receptor_types], use.names = FALSE))
    sender_genes <- unique(unlist(ex[ct$sender], use.names = FALSE))
    edges <- lr[lr$from %in% sender_genes & lr$to %in% receptors, , drop = FALSE]
    candidates <- unique(edges$from)
    if (!all(candidates %in% colnames(ltm))) stop("LR ligands absent from ligand-target model; network versions are inconsistent")
    edges <- edges[edges$from %in% candidates, , drop = FALSE]
    if (!length(candidates)) stop("No candidate ligands: ", ct$id)
    z <- subset(x, cells = colnames(x)[ids %in% receiver_types])
    SeuratObject::Idents(z) <- if (is_global) as.character(z$group) else as.character(z$celltype_fine)
    de <- Seurat::FindMarkers(z, ident.1 = if (is_global) ct$condition else ct$receiver,
           ident.2 = ct$reference, assay = "RNA", test.use = "wilcox", min.pct = cfg$nichenet$expression_pct,
           logfc.threshold = cfg$nichenet$marker_prefilter_lfc, only.pos = FALSE, verbose = FALSE)
    de$gene <- rownames(de)
    keep <- !is.na(de$p_val_adj) & de$p_val_adj <= cfg$nichenet$target_padj &
      if (ct$geneset == "up") de$avg_log2FC >= ct$lfc else abs(de$avg_log2FC) >= ct$lfc
    targets <- intersect(de$gene[keep], rownames(ltm))
    if (!length(targets) || !length(setdiff(background, targets))) stop("Empty target or negative background set: ", ct$id)
    ct$expression_pct <- cfg$nichenet$expression_pct
    ct$target_padj <- cfg$nichenet$target_padj
    ct$target_padj_method <- "Seurat p_val_adj: Bonferroni over all assay genes; not BH"
    ct$marker_prefilter_lfc <- cfg$nichenet$marker_prefilter_lfc
    ct$activity_metric <- "aupr_corrected"; ct$ligand_significance_threshold <- NULL
    ct$top_targets <- cfg$nichenet$top_targets; ct$visualization_quantile <- cfg$nichenet$visualization_quantile
    ct$expression_scope <- "Identity-specific prevalence, Acute/Chronic pooled"
    ct$background_scope <- "Union of compared receiver identities expressed genes intersect model targets"
    ct$top_ligand_ties <- "All ligands tied at the requested top-N score are retained"
    activities <- nichenetr::predict_ligand_activities(geneset = targets, background_expressed_genes = background,
                       ligand_target_matrix = ltm, potential_ligands = candidates)
    activities <- as.data.frame(activities)
    activities <- activities[order(-activities$aupr_corrected, activities$test_ligand), ]
    if (!nrow(activities) || any(!is.finite(activities$aupr_corrected))) stop("Invalid ligand activity scores: ", ct$id)
    prevalence <- length(targets) / length(union(background, targets))
    if (max(abs(activities$aupr - prevalence - activities$aupr_corrected)) > 1e-10) stop("Unexpected corrected AUPR definition")
    activities$rank <- rank(-activities$aupr_corrected)
    cutoff <- activities$aupr_corrected[min(ct$top_ligands, nrow(activities))]
    best <- activities$test_ligand[activities$aupr_corrected >= cutoff]
    links <- dplyr::bind_rows(lapply(best, nichenetr::get_weighted_ligand_target_links,
                geneset = targets, ligand_target_matrix = ltm, n = cfg$nichenet$top_targets))
    links <- tidyr::drop_na(links)
    ct$wrapper_crosscheck <- isTRUE(cfg$nichenet$wrapper_crosscheck)
    if (ct$wrapper_crosscheck) {
      native <- original_nichenet_wrapper(x, ct, ltm, lr, weighted, cfg$nichenet)
      if (!setequal(native$geneset_oi, targets) || !setequal(native$background_expressed_genes, background) ||
          !setequal(native$top_ligands, best) || !setequal(native$ligand_activities$test_ligand, activities$test_ligand)) stop("Original-wrapper input/selection mismatch: ", ct$id)
      scores <- c("auroc", "aupr", "aupr_corrected", "pearson", "rank")
      old <- as.matrix(native$ligand_activities[match(activities$test_ligand, native$ligand_activities$test_ligand), scores])
      current <- as.matrix(activities[, scores])
      if (!identical(unname(is.na(old)), unname(is.na(current)))) stop("Original-wrapper missing-score mismatch: ", ct$id)
      delta <- max(abs(old - current), na.rm = TRUE)
      if (!is.finite(delta) || delta > 1e-10) stop("Original-wrapper score mismatch: ", ct$id)
      write_table(data.frame(targets_match = TRUE, background_match = TRUE, prioritized_ligands_match = TRUE,
                  max_absolute_score_difference = delta), file.path(folder, "original_wrapper_equivalence.csv"))
      saveRDS(native, file.path(folder, "original_wrapper_result.rds"))
      for (field in c("ligand_activity_target_heatmap", "ligand_receptor_heatmap", "ligand_expression_dotplot")) {
        if (!is.null(native[[field]])) save_plot(native[[field]], file.path(folder, paste0("original_", field, ".pdf")), width = 60 / 2.54, height = 20 / 2.54)
      }
      rm(native); gc(verbose = FALSE)
    }
    # Export every modeled target/ligand weight for every candidate, not only topN.
    # Restrict rows first to avoid duplicating the large full prior matrix in memory.
    saveRDS(ltm[targets, candidates, drop = FALSE], file.path(folder, "all_candidate_target_weights.rds"))
    vis <- NULL; vis_error <- NULL
    if (nrow(links)) vis <- tryCatch(nichenetr::prepare_ligand_target_visualization(links, ltm, cutoff = cfg$nichenet$visualization_quantile),
                    error = function(e) {vis_error <<- conditionMessage(e); NULL})
    write_table(de, file.path(folder, "receiver_DE_full.csv"))
    write_table(data.frame(gene = background), file.path(folder, "background_genes.csv"))
    write_table(data.frame(gene = targets), file.path(folder, "target_genes.csv"))
    universe <- union(background, targets)
    write_table(data.frame(gene = universe, target = universe %in% targets, expression_background = universe %in% background), file.path(folder, "evaluation_universe.csv"))
    write_table(edges, file.path(folder, "candidate_ligand_receptor_edges.csv"))
    write_table(data.frame(ligand = candidates), file.path(folder, "candidate_ligands.csv"))
    write_table(activities, file.path(folder, "ligand_activities_full.csv"))
    write_table(activities[activities$test_ligand %in% best, ], file.path(folder, "prioritized_ligands.csv"))
    write_table(links, file.path(folder, "top_target_links_before_visualization.csv"))
    if (!is.null(vis) && length(vis)) {
      utils::write.csv(vis, file.path(folder, "ligand_target_visualization_matrix.csv"))
      h <- as.data.frame(as.table(vis), stringsAsFactors = FALSE)
      names(h) <- c("target", "ligand", "weight")
      save_plot(ggplot2::ggplot(h, ggplot2::aes(x = ligand, y = target, fill = weight)) + ggplot2::geom_tile() +
                  ggplot2::scale_fill_viridis_c() + ggplot2::theme_bw() +
                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)),
                file.path(folder, "ligand_target_heatmap.pdf"), width = max(8, ncol(vis) * 0.3), height = max(6, nrow(vis) * 0.12))
    }
    cells <- data.frame(cell_id = colnames(x), sample = x$sample, group = x$group, identity = ids,
                        receiver = ids %in% receiver_types, sender = ids %in% ct$sender)
    cells <- cells[cells$receiver | cells$sender, ]
    write_table(cells, file.path(folder, "cells_used.csv"))
    write_table(as.data.frame(with(cells, table(sample = factor(sample, levels = samples$sample), group, identity))), file.path(folder, "cells_by_sample.csv"))
    jsonlite::write_json(ct, file.path(folder, "parameters.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
    saveRDS(list(parameters = ct, activities = activities, targets = targets, background = background,
                  candidates = candidates, prioritized_ligands = best, links = links, visualization = vis, visualization_error = vis_error), file.path(folder, "complete_numeric_result.rds"))
    summary[[ct$id]] <- data.frame(analysis = ct$id, targets = length(targets), background = length(background),
        targets_outside_background = length(setdiff(targets, background)), candidates = length(candidates),
        prioritized = length(best), visualization_error = if (is.null(vis_error)) "" else vis_error)
    write_table(do.call(rbind, summary), file.path(dest, "analysis_summary.csv"))
    rm(z); gc(verbose = FALSE)
  }
}
