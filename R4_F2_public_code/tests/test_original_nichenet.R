# Regression for the retained manuscript Figure 6, not the withdrawn redesign.
source("R/helpers.R")
source("R/contrasts.R")
source("config.R")
stopifnot(identical(config$nichenet$contrast_set, "manuscript"))
cts <- make_contrasts("manuscript")
stopifnot(length(cts) == 10L)
fig6 <- make_contrasts("figure6")
stopifnot(length(fig6) == 6L)
stopifnot(identical(vapply(fig6, `[[`, character(1), "receiver"),
                   c("IS_I_Inflam", "IFN_response_Neu", "IS_I_Neu", "Cxcl12+ iFs", "Slpi+ iFs", "lY6C1+ iFs")))
stopifnot(identical(vapply(fig6, `[[`, character(1), "reference"), c(rep("IS_I_like", 3), rep("vFs", 3))))
for (i in seq_along(fig6)) {
  a <- fig6[[i]]
  stopifnot(a$geneset == "DE", a$lfc == 0.25, a$top_ligands == 15L, !a$affected_receptors_only)
  stopifnot(identical(a$sender, if (i <= 3) c("Cxcl12+ iFs", "lY6C1+ iFs", "Slpi+ iFs") else c("IS_I_Inflam", "IS_I_Neu", "IFN_response_Neu")))
}
fig5 <- make_contrasts("figure5")
stopifnot(length(fig5) == 4L, all(vapply(fig5, function(a) a$geneset == "up" && a$lfc == 0.4 && a$top_ligands == 30L && length(a$sender) == 8L, logical(1))))
for (old in c("revision_figure6", "recorded")) {
  err <- tryCatch({make_contrasts(old); NULL}, error = identity)
  stopifnot(inherits(err, "error"), grepl("retired", conditionMessage(err)))
}
legacy <- c("ISG_I_Inflam", "ISG_I_Neu", "ISG_I_like", "IFN_response_Neu")
stopifnot(identical(canonical_neutrophil_labels(legacy), c("IS_I_Inflam", "IS_I_Neu", "IS_I_like", "IFN_response_Neu")))
stopifnot(identical(legacy, c("ISG_I_Inflam", "ISG_I_Neu", "ISG_I_like", "IFN_response_Neu")))
neu_md <- data.frame(Neu_subtype = legacy, row.names = paste0("legacy", 1:4))
global_md <- data.frame(celltype = rep("Neutrophils", 4), row.names = rownames(neu_md))
fib_md <- data.frame(celltype = character())
stopifnot(identical(unname(project_labels(global_md, fib_md, neu_md)), canonical_neutrophil_labels(legacy)),
          identical(neu_md$Neu_subtype, legacy))
table <- utils::read.csv("metadata/figure6_contrasts.csv", check.names = FALSE)
stopifnot(identical(table$receiver, vapply(fig6, `[[`, character(1), "receiver")),
          identical(table$reference, vapply(fig6, `[[`, character(1), "reference")))
cat("PASS: retained Figure 6 A-F, all ten manuscript comparisons, three-sender pools, DE defaults, retired-design rejection and legacy aliases\n")
