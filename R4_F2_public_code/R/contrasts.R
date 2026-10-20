make_contrasts <- function(set = "manuscript") {
  if (set %in% c("revision_figure6", "recorded")) stop("This comparison set is retired. Use manuscript, figure5 or figure6 for the retained manuscript figures.")
  if (!set %in% c("manuscript", "figure5", "figure6")) stop("Unknown NicheNet contrast_set")
  all_global <- c("Epithelial cells", "Endothelial cells", "Myeloid cells", "T cells", "B cells",
                  "Neutrophils", "Fibroblasts", "Plasma cell", "ILC")
  globals <- lapply(seq_len(4), function(i) {
    rec <- if (i <= 2) "Neutrophils" else "Fibroblasts"
    condition <- if (i %% 2 == 1) "Acute" else "Chronic"
    list(id = paste0("Fig5", LETTERS[i], "_", rec, "_", condition), mode = "global", receiver = rec,
         condition = condition, reference = setdiff(c("Acute", "Chronic"), condition),
         sender = setdiff(all_global, rec), geneset = "up", lfc = 0.4, top_ligands = 30L,
         affected_receptors_only = FALSE, provenance = "Original script's four global calls; ordered by manuscript panels")
  })
  receivers <- c("IS_I_Inflam", "IFN_response_Neu", "IS_I_Neu", "Cxcl12+ iFs", "Slpi+ iFs", "lY6C1+ iFs")
  ids <- c("Fig6A_IS_I_Inflam", "Fig6B_IFN_response_Neu", "Fig6C_IS_I_Neu", "Fig6D_Cxcl12_iFs", "Fig6E_Slpi_iFs", "Fig6F_Ly6c1_iFs")
  subtypes <- lapply(seq_along(receivers), function(i) list(
    id = ids[i], mode = "subtype", receiver = receivers[i],
    reference = if (i <= 3) "IS_I_like" else "vFs",
    sender = if (i <= 3) c("Cxcl12+ iFs", "lY6C1+ iFs", "Slpi+ iFs") else c("IS_I_Inflam", "IS_I_Neu", "IFN_response_Neu"),
    geneset = "DE", lfc = 0.25, top_ligands = 15L, affected_receptors_only = FALSE,
    provenance = "Retained original Figure 6: script call patterns plus original panel definitions; A confirmed as IS-I inflammatory; stages pooled"))
  switch(set, manuscript = c(globals, subtypes), figure5 = globals, figure6 = subtypes)
}
