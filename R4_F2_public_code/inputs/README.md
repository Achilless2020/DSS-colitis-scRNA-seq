No expression data, animal identifiers beyond the six analysis sample aliases, or network binaries are bundled.

Raw route: put each sample's matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz in counts/S1 through counts/S6, or configure private paths.

Checkpoint route: provide global.rds (sce), fibro.rds (sce.F), neu.rds (sce.neu), with raw RNA counts, normalized RNA data, sample/group metadata and annotated subtypes. UMAP is required by summary/Monocle3. Run tools/export_workspace.R explicitly to convert a local RData if needed.

Networks: place the three mouse RDS files listed in metadata/network_hashes.csv in networks/. Obtain them from the NicheNet authors' resource distribution; respect their data terms. Their hashes identify the model used in the revision. Binary files are deliberately excluded from the code archive.

Optional GSEA: supply m5.go.v2024.1.Mm.symbols.gmt under the MSigDB terms. The package does not redistribute it or silently substitute a newer gene-set collection.
