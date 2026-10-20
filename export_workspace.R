#!/usr/bin/env Rscript
# Explicit local conversion; no source file is modified, no data are published.
# Rscript --vanilla tools/export_workspace.R /private/author_workspace.RData /private/checkpoints
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: export_workspace.R WORKSPACE.RData OUTPUT_DIRECTORY")
if (!file.exists(args[1])) stop("Workspace does not exist")
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
targets <- file.path(args[2], paste0(c("global", "fibro", "neu"), ".rds"))
if (any(file.exists(targets))) stop("Refusing to overwrite existing checkpoints; choose a new output directory")
e <- new.env(parent = baseenv())
loaded <- load(args[1], envir = e)
object_names <- c("sce", "sce.F", "sce.neu")
if (!all(object_names %in% loaded)) stop("Expected sce, sce.F, sce.neu in workspace")
for (i in seq_along(object_names)) saveRDS(e[[object_names[i]]], targets[i])
writeLines(c("Converted saved Seurat objects only; historical upstream steps not reexecuted.",
             capture.output(sessionInfo())), file.path(args[2], "export_session.txt"))
cat("Exported three local RDS files. Keep all expression data out of the public code archive.\n")
