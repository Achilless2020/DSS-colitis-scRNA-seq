cluster_object <- function(x, settings, seed, all_genes_scale = FALSE) {
  x <- join_rna(x)
  x <- Seurat::NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x <- Seurat::ScaleData(x, features = if (all_genes_scale) rownames(x) else SeuratObject::VariableFeatures(x), verbose = FALSE)
  npcs <- min(settings$npcs, ncol(x) - 1L, length(SeuratObject::VariableFeatures(x)) - 1L)
  if (npcs < max(c(settings$neighbor_dims, settings$umap_dims))) stop("Insufficient cells/variable genes for configured PCA dimensions")
  x <- Seurat::RunPCA(x, npcs = npcs, seed.use = seed, verbose = FALSE)
  x <- Seurat::FindNeighbors(x, dims = settings$neighbor_dims, verbose = FALSE)
  x <- Seurat::FindClusters(x, resolution = settings$resolutions, algorithm = 1, random.seed = seed, verbose = FALSE)
  x <- Seurat::RunUMAP(x, reduction = "pca", dims = settings$umap_dims, n.neighbors = min(30L, ncol(x) - 1L), seed.use = seed, verbose = FALSE)
  col <- paste0("RNA_snn_res.", settings$final_resolution)
  if (!col %in% names(x[[]])) stop("Final resolution was not computed")
  SeuratObject::Idents(x) <- col
  x
}
annotate_object <- function(x, settings, label, name, cfg, root, out) {
  col <- paste0("RNA_snn_res.", settings$final_resolution)
  save_object(x, out, paste0(name, "_clusters_unannotated"))
  write_table(data.frame(cell_id = colnames(x), sample = x$sample, cluster = as.character(x[[]][[col]])), file.path(out, "seurat", paste0(name, "_clusters.csv")))
  markers <- Seurat::FindAllMarkers(x, only.pos = TRUE, test.use = "wilcox",
            min.pct = cfg$markers$min_pct, logfc.threshold = cfg$markers$logfc, verbose = FALSE)
  write_table(markers, file.path(out, "seurat", paste0(name, "_cluster_markers_full.csv")))
  save_plot(Seurat::DimPlot(x, group.by = col, label = TRUE), file.path(out, "seurat", paste0(name, "_clusters.pdf")))
  save_plot(Seurat::DimPlot(x, group.by = "sample"), file.path(out, "seurat", paste0(name, "_samples.pdf")))
  if (!isTRUE(settings$annotation_reviewed)) stop(name, ": cluster outputs saved. Review markers and metadata mapping; set this population's annotation_reviewed=TRUE only for the reviewed clustering. Original numeric cluster IDs are not transferable automatically.")
  mapping <- utils::read.csv(path_from(root, settings$map), colClasses = "character")
  x[[label]] <- map_clusters(x[[]][[col]], mapping)
  # Derive levels from all mapped labels, retaining cluster 4 IS_I_Inflam.
  x[[label]] <- factor(x[[]][[label]], levels = unique(mapping$label))
  SeuratObject::Idents(x) <- label
  save_object(x, out, name)
  x
}
stage_global <- function(cfg, root, out, samples) {
  if (cfg$input_mode != "raw") stop("global reclustering requires raw mode")
  x <- cluster_object(read_object(out, "qc"), cfg$global, cfg$seed)
  x <- annotate_object(x, cfg$global, "celltype", "global", cfg, root, out)
  audit_counts(x, "global_annotated", out, samples)
}
stage_subcluster <- function(cfg, root, out, samples) {
  if (cfg$input_mode != "raw") stop("subcluster rerun requires raw mode; checkpoint mode preserves saved subclusters")
  x <- read_object(out, "global")
  validate_object(x, samples, "celltype")
  for (name in c("fibro", "neu")) {
    parent <- if (name == "fibro") "Fibroblasts" else "Neutrophils"
    cells <- colnames(x)[as.character(x$celltype) == parent]
    if (!length(cells)) stop("Parent population absent: ", parent)
    z <- subset(x, cells = cells)
    audit_counts(z, paste0(name, "_before_secondary_qc"), out, samples)
    if (name == "neu") {
      low <- max(200, stats::median(z$nFeature_RNA) * 0.25)
      high <- stats::median(z$nFeature_RNA) * 2
      count_high <- stats::median(z$nCount_RNA) * 3
      keep <- colnames(z)[z$nFeature_RNA > low & z$nFeature_RNA < high & z$nCount_RNA < count_high & z$percent_mito < 10 & z$percent_hb < 5]
      if (!length(keep)) stop("No neutrophils survive secondary QC")
      write_table(data.frame(min_features_exclusive = low, max_features_exclusive = high, max_counts_exclusive = count_high,
                     max_mito_exclusive = 10, max_hb_exclusive = 5), file.path(out, "qc", "neu_secondary_parameters.csv"))
      z <- subset(z, cells = keep)
    }
    audit_counts(z, paste0(name, "_after_secondary_qc"), out, samples)
    z <- cluster_object(z, cfg[[name]], cfg$seed, all_genes_scale = name == "fibro")
    annotate_object(z, cfg[[name]], if (name == "neu") "Neu_subtype" else "celltype", name, cfg, root, out)
  }
}
stage_summary <- function(cfg, root, out, samples) {
  objs <- lapply(c("global", "fibro", "neu"), function(nm) read_object(out, nm))
  names(objs) <- c("global", "fibro", "neu")
  project_labels(objs$global[[]], objs$fibro[[]], objs$neu[[]])
  for (nm in names(objs)) {
    x <- join_rna(objs[[nm]])
    label <- if (nm == "neu") "Neu_subtype" else "celltype"
    validate_object(x, samples, label)
    ids <- checked_labels(x[[]][[label]], label)
    tab <- as.data.frame(table(sample = factor(x$sample, levels = samples$sample), identity = factor(ids, levels = sort(unique(ids)))))
    tab$group <- samples$group[match(tab$sample, samples$sample)]
    denom <- ave(tab$Freq, tab$sample, FUN = sum)
    tab$denominator <- denom
    tab$proportion <- ifelse(denom > 0, tab$Freq / denom, NA_real_)
    tab$denominator_definition <- paste0("All cells retained in the ", nm, " object; no identity cleanup implied")
    write_table(tab, file.path(out, "summary", paste0(nm, "_per_sample_proportions.csv")))
    write_table(data.frame(cell_id = colnames(x), x[[]], check.names = FALSE), file.path(out, "summary", paste0(nm, "_cell_metadata.csv")))
    save_plot(Seurat::DimPlot(x, group.by = label, label = TRUE), file.path(out, "summary", paste0(nm, "_identities.pdf")))
    save_plot(Seurat::DimPlot(x, group.by = "group"), file.path(out, "summary", paste0(nm, "_groups.pdf")))
    SeuratObject::Idents(x) <- label
    markers <- Seurat::FindAllMarkers(x, assay = "RNA", test.use = "wilcox", only.pos = TRUE,
                 min.pct = cfg$markers$min_pct, logfc.threshold = cfg$markers$logfc, verbose = FALSE)
    write_table(markers, file.path(out, "summary", paste0(nm, "_identity_markers_full.csv")))
  }
  # Preserve coarse labels; show every refinement so the 263-cell issue remains visible.
  fine <- project_labels(objs$global[[]], objs$fibro[[]], objs$neu[[]])
  md <- objs$global[[]]
  write_table(data.frame(cell_id = names(fine), sample = md$sample, coarse = md$celltype, fine = unname(fine)), file.path(out, "summary", "coarse_to_fine_audit.csv"))
  write_table(as.data.frame(table(coarse = md$celltype, fine = fine)), file.path(out, "summary", "coarse_to_fine_counts.csv"))
}
