require_packages <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing R packages: ", paste(absent, collapse = ", "),
                           ". See environment/README.md.")
}
path_from <- function(root, path) {
  if (grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)) path else file.path(root, path)
}
write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
checked_labels <- function(x, what) {
  if (is.null(x) || anyNA(x) || any(!nzchar(as.character(x)))) stop(what, ": missing labels")
  as.character(x)
}
map_clusters <- function(clusters, mapping) {
  stopifnot(all(c("cluster", "label") %in% names(mapping)))
  if (anyDuplicated(mapping$cluster)) stop("duplicate cluster mapping")
  result <- as.character(mapping$label[match(as.character(clusters), as.character(mapping$cluster))])
  if (anyNA(result) || any(!nzchar(result))) stop("unmapped clusters: ", paste(setdiff(clusters, mapping$cluster), collapse = ", "))
  result
}
canonical_neutrophil_labels <- function(values) {
  values <- as.character(values)
  aliases <- c(ISG_I_Inflam = "IS_I_Inflam", ISG_I_Neu = "IS_I_Neu", ISG_I_like = "IS_I_like")
  matched <- values %in% names(aliases)
  values[matched] <- unname(aliases[values[matched]])
  values
}
project_labels <- function(global_md, fib_md, neu_md) {
  labels <- setNames(checked_labels(global_md$celltype, "global"), rownames(global_md))
  if (length(intersect(rownames(fib_md), rownames(neu_md)))) stop("subtype barcode overlap")
  for (item in list(list(md = fib_md, col = "celltype"), list(md = neu_md, col = "Neu_subtype"))) {
    md <- item$md
    if (anyDuplicated(rownames(md)) || !all(rownames(md) %in% names(labels))) stop("invalid subset barcode")
    values <- checked_labels(md[[item$col]], item$col)
    if (item$col == "Neu_subtype") values <- canonical_neutrophil_labels(values)
    if (length(values) != nrow(md)) stop("missing subtype labels")
    labels[rownames(md)] <- values
  }
  labels
}
root_near_centroid <- function(umap, candidates) {
  if (!length(candidates) || !all(candidates %in% rownames(umap))) stop("invalid root candidates")
  # Search ONLY within the requested root state; stable lexical tie break.
  candidates <- sort(candidates)
  u <- umap[candidates, , drop = FALSE]
  distances <- rowSums(sweep(u, 2, colMeans(u), "-")^2)
  candidates[which.min(distances)]
}
expressed_by_identity <- function(data, identities, celltypes, pct = 0.1) {
  stopifnot(ncol(data) == length(identities), !anyNA(identities), pct > 0, pct <= 1)
  result <- lapply(celltypes, function(z) {
    cells <- which(identities == z)
    if (!length(cells)) stop("missing identity: ", z)
    rownames(data)[Matrix::rowSums(data[, cells, drop = FALSE] > 0) / length(cells) >= pct]
  })
  stats::setNames(result, celltypes)
}
join_rna <- function(x) {
  SeuratObject::DefaultAssay(x) <- "RNA"
  if (inherits(x[["RNA"]], "Assay5")) x <- SeuratObject::JoinLayers(x, assay = "RNA")
  x
}
validate_object <- function(x, samples, label = NULL) {
  if (!inherits(x, "Seurat")) stop("Expected a Seurat object")
  if (!"RNA" %in% names(x@assays)) stop("RNA assay missing")
  md <- x[[]]
  if (anyDuplicated(colnames(x)) || !identical(rownames(md), colnames(x))) stop("barcode order mismatch")
  if (!all(c("sample", "group") %in% names(md))) stop("sample/group metadata missing")
  if (!all(md$sample %in% samples$sample)) stop("unknown sample")
  expected <- samples$group[match(md$sample, samples$sample)]
  if (anyNA(md$group) || any(as.character(md$group) != expected)) stop("sample/group mapping mismatch")
  counts <- SeuratObject::GetAssayData(x, assay = "RNA", layer = "counts")
  if (!identical(colnames(counts), colnames(x))) stop("RNA count barcode order mismatch")
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) stop("invalid gene names")
  if (!is.null(label)) checked_labels(md[[label]], label)
  invisible(TRUE)
}
audit_counts <- function(x, stage, out, samples) {
  tab <- table(factor(x$sample, levels = samples$sample))
  z <- data.frame(stage = stage, sample = samples$sample, group = samples$group,
                  cells = as.integer(tab), genes = nrow(x))
  path <- file.path(out, "qc", "cell_counts_by_step.csv")
  if (file.exists(path)) {
    prev <- utils::read.csv(path)
    z <- rbind(prev[prev$stage != stage, ], z)
  }
  write_table(z, path)
}
save_plot <- function(p, path, width = 9, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = p, width = width, height = height, limitsize = FALSE)
}
read_object <- function(out, name) readRDS(file.path(out, "objects", paste0(name, ".rds")))
save_object <- function(x, out, name) {
  dir.create(file.path(out, "objects"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file.path(out, "objects", paste0(name, ".rds")))
}
sha256_file <- function(path) {
  # Stream an actual read connection. Windows sandbox ACL metadata can make
  # file.access() return -1 for files that are successfully readable.
  con <- file(path, open = "rb")
  on.exit(close(con))
  unclass(as.character(openssl::sha256(con)))
}
hash_files <- function(paths, root = NULL) {
  if (!all(file.exists(paths))) stop("input file missing: ", paste(paths[!file.exists(paths)], collapse = ", "))
  data.frame(file = if (is.null(root)) basename(paths) else substring(normalizePath(paths, winslash = "/"), nchar(normalizePath(root, winslash = "/")) + 2L),
             bytes = file.info(paths)$size,
             sha256 = vapply(paths, sha256_file, character(1)))
}
