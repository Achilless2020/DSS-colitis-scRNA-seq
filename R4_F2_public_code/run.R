#!/usr/bin/env Rscript
# Usage: Rscript --vanilla run.R --config config.local.R --stages import,summary,monocle,nichenet
args <- commandArgs(trailingOnly = TRUE)
arg <- function(key, default) {
  at <- match(key, args)
  if (is.na(at)) return(default)
  if (at == length(args)) stop("Missing value for ", key)
  args[at + 1L]
}
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
source(file.path(root, "R/helpers.R"))
cfg_env <- new.env(parent = baseenv())
sys.source(path_from(root, arg("--config", "config.R")), envir = cfg_env)
cfg <- cfg_env$config
if (length(cfg$library_paths)) .libPaths(c(cfg$library_paths, .libPaths()))
default_stages <- if (cfg$input_mode == "raw") "import,qc,global,subcluster,summary,monocle,nichenet" else "import,summary,monocle,nichenet"
stages <- strsplit(arg("--stages", default_stages), ",", fixed = TRUE)[[1]]
valid <- c("check", "import", "qc", "global", "subcluster", "summary", "monocle", "nichenet", "de", "enrichment")
if (!all(stages %in% valid)) stop("Unknown stage; use: ", paste(valid, collapse = ", "))
require_packages(c("Seurat", "Matrix", "ggplot2", "openssl", "jsonlite", "future"))
future::plan("sequential")
options(warn = 1)
set.seed(cfg$seed)
out <- path_from(root, cfg$output_dir)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
samples <- utils::read.csv(path_from(root, cfg$samples), stringsAsFactors = FALSE)
stopifnot(all(c("sample", "group", "matrix_dir", "expected_doublet_rate") %in% names(samples)),
          !anyDuplicated(samples$sample), !anyNA(samples$sample), !anyNA(samples$group),
          all(samples$group %in% c("Acute", "Chronic")))
for (file in c("01_input_qc.R", "02_seurat.R", "03_monocle.R", "contrasts.R", "04_nichenet.R", "05_de_enrichment.R")) source(file.path(root, "R", file))
writeLines(capture.output(dput(cfg)), file.path(out, "effective_config.R"))
writeLines(capture.output(sessionInfo()), file.path(out, "session_start.txt"))
for (stage in stages) {
  message("START ", stage)
  status <- "running"
  writeLines(capture.output(dput(cfg)), file.path(out, paste0("config_", stage, ".R")))
  write_table(data.frame(stage = stage, status = status, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(out, paste0("status_", stage, ".csv")))
  tryCatch({
    switch(stage,
      check = {require_packages(c("DoubletFinder", "monocle3", "nichenetr", "SingleCellExperiment", "SummarizedExperiment")); message("Core dependencies available; no analysis executed.")},
      import = stage_import(cfg, root, out, samples),
      qc = stage_qc(cfg, root, out, samples),
      global = stage_global(cfg, root, out, samples),
      subcluster = stage_subcluster(cfg, root, out, samples),
      summary = stage_summary(cfg, root, out, samples),
      monocle = stage_monocle(cfg, root, out, samples),
      nichenet = stage_nichenet(cfg, root, out, samples),
      de = stage_de(cfg, root, out, samples),
      enrichment = stage_enrichment(cfg, root, out, samples))
    status <- "completed"
  }, error = function(e) {
    write_table(data.frame(stage = stage, status = "failed", time = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(out, paste0("status_", stage, ".csv")))
    writeLines(c(paste("FAILED stage:", stage), conditionMessage(e), capture.output(sessionInfo())), file.path(out, "last_failure.txt"))
    stop(conditionMessage(e), call. = FALSE)
  })
  write_table(data.frame(stage = stage, status = status, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(out, paste0("status_", stage, ".csv")))
  message("DONE ", stage)
}
writeLines(capture.output(sessionInfo()), file.path(out, "session_end.txt"))
