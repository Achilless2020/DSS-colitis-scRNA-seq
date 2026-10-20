# Run only after verifying the analysis environment. Records observations, not a fabricated lockfile.
packages <- c("Seurat", "SeuratObject", "Matrix", "DoubletFinder", "monocle3", "nichenetr", "SingleCellExperiment",
              "SummarizedExperiment", "ggplot2", "dplyr", "tidyr", "future", "openssl", "jsonlite",
              "edgeR", "clusterProfiler", "org.Mm.eg.db")
db <- installed.packages()
deps <- tools::package_dependencies(packages, db = db, recursive = TRUE)
all_packages <- sort(unique(c(packages, unlist(deps))))
fields <- c("Package", "Version", "Repository", "RemoteHost", "RemoteUsername", "RemoteRepo", "RemoteRef", "RemoteSha", "Built")
entries <- lapply(all_packages, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(NULL)
  d <- read.dcf(file.path(find.package(p), "DESCRIPTION"))
  vals <- stats::setNames(lapply(fields, function(f) if (f %in% colnames(d)) unname(d[1, f]) else ""), fields)
  as.data.frame(vals, stringsAsFactors = FALSE)
})
utils::write.csv(do.call(rbind, entries), "environment/observed_packages.csv", row.names = FALSE)
writeLines(c(R.version.string, "Observed installed packages and dependency origins; no clean-environment restoration claim."), "environment/runtime.txt")
# sessionInfo may include user-specific library paths; the public record uses package metadata instead.
cat("Recorded", sum(vapply(entries, Negate(is.null), logical(1))), "installed package versions.\n")
