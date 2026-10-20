make_contrasts <- function(set) {
  all_global <- c("Epithelial cells", "Endothelial cells", "Myeloid cells", "T cells", "B cells",
                  "Neutrophils", "Fibroblasts", "Plasma cell", "ILC")
  neu <- c("IS_I_like", "IFN_response_Neu", "IS_I_Neu")
  fib <- c("Cxcl12+ iFs", "Slpi+ iFs", "lY6C1+ iFs")
  if (set == "revision_figure6") {
    receivers <- c(neu, fib)
    ids <- c("Fig6A_IS_I_like", "Fig6B_IFN_response_Neu", "Fig6C_IS_I_Neu", "Fig6D_Cxcl12_iFs", "Fig6E_Slpi_iFs", "Fig6F_Ly6c1_iFs")
    return(lapply(seq_along(receivers), function(i) list(id = ids[i], mode = "subtype",
      receiver = receivers[i], reference = setdiff(if (i <= 3) neu else fib, receivers[i]),
      sender = if (i <= 3) c(fib, "vFs") else c(neu, "IS_I_Inflam", "Immuno_silent_Neu", "Immature_Neu"),
      geneset = "up", lfc = 0.25, top_ligands = 15L, affected_receptors_only = TRUE,
      provenance = "New revision contrast; stages pooled; not historical reconstruction")))
  }
  if (set != "recorded") stop("Unknown NicheNet contrast_set")
  globals <- lapply(seq_len(4), function(i) {
    rec <- if (i <= 2) "Neutrophils" else "Fibroblasts"
    condition <- if (i %% 2 == 1) "Acute" else "Chronic"
    list(id = paste0("Fig5", LETTERS[i], "_", rec, "_", condition), mode = "global", receiver = rec,
         condition = condition, reference = setdiff(c("Acute", "Chronic"), condition),
         sender = setdiff(all_global, rec), geneset = "up", lfc = 0.4, top_ligands = 30L,
         affected_receptors_only = FALSE, provenance = "Comparison explicitly recorded in original code")
  })
  c(globals, list(
    list(id = "recorded_IFN_vs_IS_I_like", mode = "subtype", receiver = "IFN_response_Neu", reference = "IS_I_like",
         sender = fib, geneset = "DE", lfc = 0.25, top_ligands = 15L, affected_receptors_only = FALSE,
         provenance = "Original wrapper default: two-sided DE targets; stages pooled"),
    list(id = "recorded_Slpi_vs_vFs", mode = "subtype", receiver = "Slpi+ iFs", reference = "vFs",
         sender = c("IS_I_Inflam", "IS_I_Neu", "IFN_response_Neu"), geneset = "DE", lfc = 0.25,
         top_ligands = 15L, affected_receptors_only = FALSE,
         provenance = "Original wrapper default: two-sided DE targets; stages pooled")))
}
