stage_import <- function(cfg, root, out, samples) {
  if (cfg$input_mode == "checkpoints") {
    paths <- vapply(cfg$objects[c("global", "fibro", "neu")], function(p) path_from(root, p), character(1))
    write_table(hash_files(paths), file.path(out, "input_manifest.csv"))
    objects <- lapply(paths, readRDS)
    for (nm in names(objects)) {
      x <- join_rna(objects[[nm]])
      validate_object(x, samples, if (nm == "neu") "Neu_subtype" else "celltype")
      if (nm != "global" && !all(colnames(x) %in% colnames(objects$global))) stop("subset cells absent from global checkpoint")
      save_object(x, out, nm)
      audit_counts(x, paste0("checkpoint_", nm), out, samples)
    }
    project_labels(objects$global[[]], objects$fibro[[]], objects$neu[[]])
    writeLines("Saved annotated objects imported. Raw QC/doublet removal/clustering were NOT rerun.", file.path(out, "input_scope.txt"))
    return(invisible(NULL))
  }
  if (cfg$input_mode != "raw") stop("input_mode must be raw or checkpoints")
  all_files <- character()
  objects <- lapply(seq_len(nrow(samples)), function(i) {
    dir <- path_from(root, samples$matrix_dir[i])
    files <- file.path(dir, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
    if (!all(file.exists(files))) stop("Missing count files for ", samples$sample[i])
    all_files <<- c(all_files, files)
    m <- Matrix::readMM(gzfile(files[1]))
    f <- utils::read.delim(gzfile(files[2]), header = FALSE, stringsAsFactors = FALSE)
    b <- utils::read.delim(gzfile(files[3]), header = FALSE, stringsAsFactors = FALSE)
    if (cfg$feature_column > ncol(f)) stop("feature_column not present")
    genes <- as.character(f[[cfg$feature_column]])
    if (length(genes) != nrow(m) || nrow(b) != ncol(m)) stop("matrix/feature/barcode dimension mismatch")
    if (anyNA(genes) || anyDuplicated(genes) || any(!nzchar(genes))) stop("Gene symbols must be nonempty and unique; supply a documented duplicate-gene aggregation upstream")
    if (mean(grepl("^ENSMUSG", genes)) > 0.5) stop("Ensembl IDs detected. Select a verified gene-symbol column or supply a documented symbol mapping before QC/NicheNet.")
    if (anyDuplicated(b[[1]])) stop("duplicate input barcodes")
    rownames(m) <- genes
    colnames(m) <- paste(samples$sample[i], b[[1]], sep = "_")
    x <- SeuratObject::CreateSeuratObject(counts = as(m, "CsparseMatrix"), project = samples$sample[i])
    x$sample <- samples$sample[i]; x$group <- samples$group[i]
    x
  })
  names(objects) <- samples$sample
  x <- if (length(objects) == 1) objects[[1]] else merge(objects[[1]], y = objects[-1])
  x <- join_rna(x)
  validate_object(x, samples)
  save_object(x, out, "raw")
  # Keep sample-relative filenames to distinguish the six matrix triplets.
  manifest <- hash_files(all_files)
  manifest$file <- paste(rep(samples$sample, each = 3), basename(all_files), sep = "/")
  write_table(manifest, file.path(out, "input_manifest.csv"))
  audit_counts(x, "raw", out, samples)
}
stage_qc <- function(cfg, root, out, samples) {
  if (cfg$input_mode != "raw") stop("qc requires raw mode; checkpoint imports must not silently reapply QC")
  x <- join_rna(read_object(out, "raw"))
  for (key in c("mito", "ribo", "hb")) {
    genes <- grep(cfg$qc[[paste0(key, "_pattern")]], rownames(x), value = TRUE, ignore.case = TRUE)
    if (!length(genes)) stop("No ", key, " genes found; inspect feature identifiers and regex")
    write_table(data.frame(gene = genes), file.path(out, "qc", paste0(key, "_genes.csv")))
    x[[paste0("percent_", key)]] <- Seurat::PercentageFeatureSet(x, features = genes)
  }
  write_table(data.frame(cell_id = colnames(x), x[[]], check.names = FALSE), file.path(out, "qc", "metadata_before_filters.csv"))
  gene_keep <- rownames(x)[Matrix::rowSums(SeuratObject::GetAssayData(x, assay = "RNA", layer = "counts") > 0) > cfg$min_cells_per_gene]
  keep <- colnames(x)[x$nFeature_RNA > cfg$qc$min_features & x$nFeature_RNA < cfg$qc$max_features]
  if (!length(keep) || !length(gene_keep)) stop("No cells/genes survive feature QC")
  x <- subset(x, cells = keep, features = gene_keep)
  audit_counts(x, "feature_filter", out, samples)
  keep <- colnames(x)[x$percent_mito < cfg$qc$max_mito]
  if (!length(keep)) stop("No cells survive mitochondrial QC")
  x <- subset(x, cells = keep); audit_counts(x, "mito_filter", out, samples)
  keep <- colnames(x)[x$percent_hb < cfg$qc$max_hb]
  if (!length(keep)) stop("No cells survive hemoglobin QC")
  x <- subset(x, cells = keep); audit_counts(x, "hb_filter", out, samples)
  # Ribo percentage is recorded, but NOT filtered: original selected_ribo was unused.
  if (cfg$doublets$enabled) {
    require_packages("DoubletFinder")
    objects <- Seurat::SplitObject(x, split.by = "sample")
    records <- list()
    objects <- lapply(names(objects), function(id) {
      z <- objects[[id]]
      z <- Seurat::NormalizeData(z, verbose = FALSE)
      z <- Seurat::FindVariableFeatures(z, nfeatures = 2000, verbose = FALSE)
      z <- Seurat::ScaleData(z, verbose = FALSE)
      z <- Seurat::RunPCA(z, npcs = 30, seed.use = cfg$seed, verbose = FALSE)
      sweep <- DoubletFinder::paramSweep(z, PCs = cfg$doublets$dims, sct = FALSE)
      bcmvn <- DoubletFinder::find.pK(DoubletFinder::summarizeSweep(sweep, GT = FALSE))
      pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      n <- ncol(z)
      rate <- switch(cfg$doublets$rate_mode,
        source_heuristic = 0.008 * (n / 1000),
        sample_manifest = samples$expected_doublet_rate[match(id, samples$sample)],
        stop("Unknown doublet rate_mode"))
      if (length(rate) != 1 || is.na(rate) || rate <= 0 || rate >= 1) stop("Invalid expected doublet rate for ", id)
      n_exp <- round(rate * n)
      if (n_exp < 1 || n_exp >= n) stop("Invalid expected doublet count for ", id)
      z <- DoubletFinder::doubletFinder(z, PCs = cfg$doublets$dims, pN = cfg$doublets$pN,
              pK = pk, nExp = n_exp, reuse.pANN = FALSE, sct = FALSE)
      col <- grep("^DF.classifications", names(z[[]]), value = TRUE)
      if (length(col) != 1) stop("Ambiguous DoubletFinder classification columns")
      z$Doublet_Class <- as.character(z[[]][[col]])
      records[[id]] <<- data.frame(sample = id, cells = n, pK = pk, pN = cfg$doublets$pN, rate = rate,
                    nExp = n_exp, rate_mode = cfg$doublets$rate_mode, homotypic_adjustment = FALSE)
      write_table(bcmvn, file.path(out, "qc", paste0(id, "_pK_sweep.csv")))
      write_table(data.frame(cell_id = colnames(z), sample = id, classification = z$Doublet_Class), file.path(out, "qc", paste0(id, "_doublets.csv")))
      subset(z, cells = colnames(z)[z$Doublet_Class == "Singlet"])
    })
    x <- if (length(objects) == 1) objects[[1]] else merge(objects[[1]], y = objects[-1])
    x <- join_rna(x)
    write_table(do.call(rbind, records), file.path(out, "qc", "doublet_parameters.csv"))
    audit_counts(x, "singlets", out, samples)
  } else {
    writeLines("DoubletFinder explicitly disabled in effective configuration.", file.path(out, "qc", "doublets_disabled.txt"))
  }
  validate_object(x, samples)
  save_object(x, out, "qc")
}
