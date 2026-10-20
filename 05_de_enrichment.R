# Optional revision analyses: these are NOT recovered outputs from the old notebook.
stage_de <- function(cfg, root, out, samples) {
  require_packages("edgeR")
  status <- list()
  write_table(data.frame(analysis = character(), status = character()), file.path(out, "de", "analysis_status.csv"))
  for (nm in c("global", "fibro", "neu")) {
    x <- join_rna(read_object(out, nm))
    label <- if (nm == "neu") "Neu_subtype" else "celltype"
    validate_object(x, samples, label)
    identities <- unique(as.character(x[[]][[label]]))
    if (nm == "fibro") identities <- intersect(identities, c("vFs", "Cxcl12+ iFs", "lY6C1+ iFs", "Slpi+ iFs"))
    for (id in identities) {
      key <- paste(nm, gsub("[^A-Za-z0-9]+", "_", id), sep = "__")
      folder <- file.path(out, "de", key)
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      cells <- colnames(x)[as.character(x[[]][[label]]) == id]
      z <- subset(x, cells = cells)
      n <- as.integer(table(factor(z$sample, levels = samples$sample)))
      use <- n >= cfg$de$min_cells_per_sample
      coverage <- data.frame(sample = samples$sample, group = samples$group, cells = n, included = use)
      write_table(coverage, file.path(folder, "sample_coverage.csv"))
      reps <- table(factor(samples$group[use], levels = c("Acute", "Chronic")))
      if (any(reps < cfg$de$min_replicates_per_group)) {
        status[[key]] <- data.frame(analysis = key, status = "skipped: insufficient biological replicates")
        next
      }
      counts <- SeuratObject::GetAssayData(z, assay = "RNA", layer = "counts")
      design_cells <- Matrix::sparseMatrix(i = seq_len(ncol(z)), j = match(z$sample, samples$sample),
                       x = 1, dims = c(ncol(z), nrow(samples)), dimnames = list(colnames(z), samples$sample))
      pb <- counts %*% design_cells
      # Dense conversion is only genes x biological samples (six columns), never genes x cells.
      pb <- as.matrix(pb[, samples$sample[use], drop = FALSE])
      saveRDS(pb, file.path(folder, "pseudobulk_counts.rds"))
      group <- factor(samples$group[use], levels = c("Acute", "Chronic"))
      design <- stats::model.matrix(~group)
      rownames(design) <- colnames(pb)
      write_table(data.frame(sample = rownames(design), design), file.path(folder, "design_matrix.csv"))
      y <- edgeR::DGEList(counts = pb, group = group)
      keep <- edgeR::filterByExpr(y, design = design)
      write_table(data.frame(gene = rownames(pb), tested = keep), file.path(folder, "gene_filter.csv"))
      if (sum(keep) < 2) stop("Insufficient tested genes for ", key)
      y <- y[keep, , keep.lib.sizes = FALSE]
      y <- edgeR::calcNormFactors(y)
      y <- edgeR::estimateDisp(y, design, robust = TRUE)
      fit <- edgeR::glmQLFit(y, design, robust = TRUE)
      test <- edgeR::glmQLFTest(fit, coef = "groupChronic")
      result <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table
      result$gene <- rownames(result)
      result$contrast <- "Chronic minus Acute"
      result$identity <- id
      write_table(result, file.path(folder, "DE_full.csv"))
      write_table(data.frame(sample = rownames(y$samples), y$samples), file.path(folder, "library_sizes_and_norm_factors.csv"))
      status[[key]] <- data.frame(analysis = key, status = "completed")
    }
  }
  write_table(do.call(rbind, status), file.path(out, "de", "analysis_status.csv"))
}
current_de_files <- function(out) {
  state_file <- file.path(out, "status_de.csv")
  if (!file.exists(state_file) || !identical(utils::read.csv(state_file)$status, "completed")) stop("The current DE stage has not completed successfully")
  status_file <- file.path(out, "de", "analysis_status.csv")
  if (!file.exists(status_file)) stop("DE analysis status missing")
  status <- utils::read.csv(status_file)
  files <- file.path(out, "de", status$analysis[status$status == "completed"], "DE_full.csv")
  if (!length(files)) stop("No completed sample-level DE outputs")
  if (!all(file.exists(files))) stop("Completed DE result missing")
  files
}
stage_enrichment <- function(cfg, root, out, samples) {
  require_packages(c("clusterProfiler", "org.Mm.eg.db"))
  gmt <- path_from(root, cfg$enrichment$gmt)
  if (!file.exists(gmt)) stop("Provide the documented mouse GO GMT file; it is not redistributed in this package")
  terms <- clusterProfiler::read.gmt(gmt)
  files <- current_de_files(out)
  write_table(hash_files(gmt), file.path(out, "enrichment", "gmt_manifest.csv"))
  for (file in files) {
    de <- utils::read.csv(file)
    folder <- file.path(out, "enrichment", basename(dirname(file)))
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    ranked <- de$logFC; names(ranked) <- de$gene
    ranked <- sort(ranked[is.finite(ranked)], decreasing = TRUE)
    if (anyDuplicated(names(ranked))) stop("Duplicate symbols in GSEA ranking")
    gsea <- clusterProfiler::GSEA(geneList = ranked, TERM2GENE = terms, pvalueCutoff = 1,
                minGSSize = cfg$enrichment$min_size, maxGSSize = cfg$enrichment$max_size,
                pAdjustMethod = "BH", seed = TRUE, verbose = FALSE)
    write_table(as.data.frame(gsea), file.path(folder, "GSEA_full.csv"))
    write_table(data.frame(gene = names(ranked), logFC = unname(ranked)), file.path(folder, "GSEA_ranked_genes.csv"))
    for (direction in c("up", "down")) {
      genes <- de$gene[!is.na(de$FDR) & de$FDR < 0.05 & if (direction == "up") de$logFC > 0 else de$logFC < 0]
      write_table(data.frame(gene = genes), file.path(folder, paste0("ORA_", direction, "_input.csv")))
      if (!length(genes)) {
        writeLines("No genes pass the specified FDR/direction rule.", file.path(folder, paste0("ORA_", direction, "_empty.txt")))
        next
      }
      ego <- clusterProfiler::enrichGO(gene = genes, universe = de$gene, OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                   keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                   minGSSize = cfg$enrichment$min_size, maxGSSize = cfg$enrichment$max_size, readable = TRUE)
      write_table(as.data.frame(ego), file.path(folder, paste0("GO_ORA_", direction, "_full.csv")))
    }
  }
}
