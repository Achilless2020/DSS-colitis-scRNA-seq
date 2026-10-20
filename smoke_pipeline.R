# Synthetic integration test. This never loads study data or validates a biological conclusion.
# Run from package root: Rscript --vanilla tests/smoke_pipeline.R
source("R/helpers.R")
source("R/02_seurat.R")
require_packages(c("Seurat", "Matrix", "jsonlite", "digest", "monocle3", "nichenetr", "edgeR"))
source("config.R")
set.seed(123)
test_parent <- Sys.getenv("R4F2_TEST_DIR", unset = tempdir())
dir.create(test_parent, recursive = TRUE, showWarnings = FALSE)
work <- tempfile("r4f2_smoke_", tmpdir = test_parent)
dir.create(work)
ids <- c("IS_I_Neu", "IS_I_like", "IFN_response_Neu", "Immuno_silent_Neu", "IS_I_Inflam", "Immature_Neu",
         "vFs", "Cxcl12+ iFs", "lY6C1+ iFs", "Slpi+ iFs", "Epithelial cells", "Endothelial cells", "Myeloid cells", "T cells", "B cells", "Plasma cell", "ILC")
fine <- rep(ids, each = 60)
sample <- rep(rep(paste0("S", 1:6), each = 10), length(ids))
genes <- c(paste0("g", 1:260), "lig1", "lig2", "lig3", "lig4", "rec1", "rec2",
           "Camp", "Ltf", "Ngp", "Il1b", "Tnf", "Cxcl2", "Ifit1", "Ifit3", "Isg15", "Rsad2", "Cd274", "mt-Nd1", "Rpl3", "Hba-a1")
counts <- matrix(stats::rpois(length(genes) * length(fine), 1), nrow = length(genes), dimnames = list(genes, paste0("cell", seq_along(fine))))
for (i in seq_along(ids)) counts[(i - 1) * 10 + 1:10, fine == ids[i]] <- stats::rpois(10 * 60, 12)
chronic <- sample %in% c("S1", "S2", "S3")
counts[201:210, chronic] <- stats::rpois(10 * sum(chronic), 10)
counts[211:220, !chronic] <- stats::rpois(10 * sum(!chronic), 10)
x <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
x$sample <- sample; x$group <- ifelse(chronic, "Chronic", "Acute")
x$celltype <- ifelse(fine %in% ids[1:6], "Neutrophils", ifelse(fine %in% ids[7:10], "Fibroblasts", fine))
settings <- config$neu; settings$npcs <- 20L
x <- cluster_object(x, settings, 123L)
fib <- subset(x, cells = colnames(x)[fine %in% ids[7:10]])
fib$celltype <- fine[match(colnames(fib), colnames(x))]
fib <- cluster_object(fib, settings, 123L)
neu <- subset(x, cells = colnames(x)[fine %in% ids[1:6]])
neu$Neu_subtype <- fine[match(colnames(neu), colnames(x))]
neu <- cluster_object(neu, settings, 123L)
for (nm in c("global", "fibro", "neu")) saveRDS(switch(nm, global = x, fibro = fib, neu = neu), file.path(work, paste0(nm, ".rds")))
config$input_mode <- "checkpoints"
config$objects <- stats::setNames(file.path(work, paste0(c("global", "fibro", "neu"), ".rds")), c("global", "fibro", "neu"))
config$output_dir <- file.path(work, "outputs")
config$library_paths <- .libPaths()
config$de$min_cells_per_sample <- 5L
config$monocle$num_dim <- 20L
config$monocle$minimal_branch_len <- 2L
config$nichenet$network_dir <- file.path(work, "networks")
config$nichenet$verify_hashes <- FALSE # Synthetic matrices cannot match study-model hashes.
dir.create(config$nichenet$network_dir)
ltm <- matrix(stats::runif(length(genes) * 4), nrow = length(genes), dimnames = list(genes, paste0("lig", 1:4)))
saveRDS(ltm, file.path(config$nichenet$network_dir, "ligand_target_matrix_nsga2r_final_mouse.rds"))
saveRDS(data.frame(from = paste0("lig", 1:4), to = rep(c("rec1", "rec2"), 2)), file.path(config$nichenet$network_dir, "lr_network_mouse_21122021.rds"))
saveRDS(list(synthetic = TRUE), file.path(config$nichenet$network_dir, "weighted_networks_nsga2r_final_mouse.rds"))
config_path <- file.path(work, "config.R")
write_cfg <- function() writeLines(c("config <-", capture.output(dput(config))), config_path)
run <- function(stages) {
  write_cfg()
  code <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "run.R", "--config", shQuote(config_path), "--stages", stages))
  if (code != 0L) stop("Smoke pipeline failed at ", stages, "; fixture retained at ", work)
}
run("import,summary,monocle,nichenet,de")
config$nichenet$contrast_set <- "recorded"
run("nichenet")
both <- utils::read.csv(file.path(config$output_dir, "monocle/pseudotime_both_roots.csv"))
roots <- utils::read.csv(file.path(config$output_dir, "monocle/root_selection_and_coverage.csv"))
stopifnot(nrow(both) == ncol(neu), nrow(roots) == 2L, all(roots$finite_cells > 0))
for (set in c("revision_figure6", "recorded")) {
  summary <- utils::read.csv(file.path(config$output_dir, "nichenet", set, "analysis_summary.csv"))
  stopifnot(nrow(summary) == 6L, all(summary$targets > 0), all(summary$candidates == 4L))
}
stopifnot(!anyNA(read_object(config$output_dir, "neu")$Neu_subtype))
stopifnot(length(list.files(file.path(config$output_dir, "de"), "^DE_full.csv$", recursive = TRUE)) == 19L)
# Raw import and QC exercise the same retained-cell chain, with explicit test thresholds.
for (i in 1:6) {
  folder <- file.path(work, "counts", paste0("S", i)); dir.create(folder, recursive = TRUE)
  cells <- which(sample == paste0("S", i))
  m <- counts[, cells, drop = FALSE]
  m["mt-Nd1", 1] <- 10000 # must be removed by mito filter
  p <- tempfile(tmpdir = folder)
  Matrix::writeMM(Matrix::Matrix(m, sparse = TRUE), p)
  con <- gzfile(file.path(folder, "matrix.mtx.gz"), "wb"); writeBin(readBin(p, "raw", n = file.info(p)$size), con); close(con)
  con <- gzfile(file.path(folder, "features.tsv.gz"), "wt"); writeLines(genes, con); close(con)
  con <- gzfile(file.path(folder, "barcodes.tsv.gz"), "wt"); writeLines(colnames(m), con); close(con)
}
manifest <- utils::read.csv("metadata/samples.csv")
manifest$matrix_dir <- file.path(work, "counts", manifest$sample)
config$samples <- file.path(work, "samples.csv")
write_table(manifest, config$samples)
config$input_mode <- "raw"; config$output_dir <- file.path(work, "raw_output")
config$qc$min_features <- 30L; config$doublets$enabled <- FALSE
run("import,qc")
raw <- read_object(config$output_dir, "raw"); qc <- read_object(config$output_dir, "qc")
forced_mito_ids <- paste(paste0("S", 1:6), colnames(counts)[vapply(1:6, function(i) which(sample == paste0("S", i))[1], integer(1))], sep = "_")
stopifnot(ncol(qc) <= ncol(raw) - 6L, !any(forced_mito_ids %in% colnames(qc)),
          all(qc$percent_mito < config$qc$max_mito), all(qc$percent_hb < config$qc$max_hb))
cat("PASS: synthetic Seurat preprocessing; raw import/QC; checkpoint CLI; summaries; both Monocle3 roots; 12 NicheNet contrasts; 19 sample-level DE analyses.\n")
cat("Synthetic fixture and outputs:", work, "\n")
