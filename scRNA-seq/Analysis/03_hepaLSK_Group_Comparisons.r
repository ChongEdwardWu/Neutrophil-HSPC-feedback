### ---------------- Section 0: Preparation -----------------------------------
rm(list = ls()); gc()

# Load personal R configuration (if applicable)
suppressMessages(source("path_to_.radian_profile"))
# .libPaths()

# Set working directory
workdir <- "path_to_data"
setwd(workdir)

# Create output folders if missing
if (!file.exists(file.path(workdir, "figures"))) dir.create(file.path(workdir, "figures"))
if (!file.exists(file.path(workdir, "results"))) dir.create(file.path(workdir, "results"))

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(openxlsx2)
  library(future)
  library(msigdbr)
  library(clusterProfiler)
})

# Optional: interactive graphics device
# httpgd::hgd()

set.seed(123)
nworkers <- 20
plan("sequential")
# plan("multicore", workers = nworkers)

# Species and group settings
species <- "mm"   # "mm" for mouse, "hs" for human
groups_order <- c("BM", "SP", "SP_NAC")  # Adjust if your groups differ

# Load annotated Seurat object
seu <- readRDS("results/02_Annotation_LSK.rds")

# Exclude unwanted genes (e.g., TCR, BCR, MT, HB, pseudogenes)
rna_genes <- rownames(seu[["RNA"]])

# Build a blacklist of low-informative or artifactual features
hist_genes   <- grep("^Hist", rna_genes, ignore.case = TRUE, value = TRUE)
hb_genes     <- grep("^Hb[ab]-|^HB(?!P)", rna_genes, perl = TRUE, value = TRUE)
mt_genes     <- grep("^mt-|^MT-", rna_genes, ignore.case = TRUE, value = TRUE)
rps_genes    <- grep("^Rp[sl]|^RP[SL]", rna_genes, ignore.case = TRUE, value = TRUE)
rik_genes    <- grep("Rik$", rna_genes, value = TRUE)
alu_genes    <- if (species == "hs") grep("^AL[0-9]+", rna_genes, value = TRUE) else character(0)
pseudo_genes <- grep("-(ps|ps[0-9]+)$", rna_genes, ignore.case = TRUE, value = TRUE)
mir_genes    <- grep("^Mir", rna_genes, ignore.case = TRUE, value = TRUE)
gencode_genes <- grep("^Gm[0-9]+", rna_genes, value = TRUE)

bad_features <- unique(c(rik_genes, gencode_genes, pseudo_genes, alu_genes))
features <- setdiff(rna_genes, bad_features)


### ---------------- Section 1: Differential regulon analysis -----------------
# Subset Seurat object for MyePro cells only
mye <- seu[, seu$CellType_l1 == "MyePro"]

# Define identity by sample group
Idents(mye) <- "group"

### ---------------- Section 2: Regulon activity (Bin assay) ------------------
# Compare regulon binary activity between groups

# (1) SP vs BM — regulons upregulated in SP
markers_MyePro_SPvsBM_bin <- FindMarkers(
  mye, ident.1 = "SP", ident.2 = "BM",
  assay = "Bin", test.use = "wilcox"
) %>%
  as.data.frame() %>%
  rownames_to_column("Regulon") %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC))

# (2) SP_NAC vs SP — regulons downregulated in SP_NAC
markers_MyePro_SPNACvsSP_bin <- FindMarkers(
  mye, ident.1 = "SP_NAC", ident.2 = "SP",
  assay = "Bin", test.use = "wilcox"
) %>%
  as.data.frame() %>%
  rownames_to_column("Regulon") %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(avg_log2FC)


### ---------------- Section 3: Identify shared regulons ----------------------
# Identify regulons that are upregulated in SP vs BM
# and downregulated in SP_NAC vs SP (i.e., shared response)

common_MyePro_bin <- intersect(
  markers_MyePro_SPvsBM_bin$Regulon,
  markers_MyePro_SPNACvsSP_bin$Regulon
)

# Combine information from both comparisons
common_MyePro_regulonBin_df <- inner_join(
  markers_MyePro_SPvsBM_bin,
  markers_MyePro_SPNACvsSP_bin,
  by = "Regulon",
  suffix = c("_SPvsBM", "_SPNACvsSP")
)


### ---------------- Section 4: Export results --------------------------------
# Export all results to an Excel file (multi-sheet)
out_file <- "results/Step3_Group_Comparisons.xlsx"
wb <- wb_workbook()

export_list <- list(
  "MyePro_SPvsBM_regulons_Bin" = markers_MyePro_SPvsBM_bin,
  "MyePro_SPNACvsSP_regulons_Bin" = markers_MyePro_SPNACvsSP_bin,
  "MyePro_common_regulons_Bin" = common_MyePro_regulonBin_df
)

for (nm in names(export_list)) {
  wb <- wb_add_worksheet(wb, sheet = nm)
  wb <- wb_add_data(wb, sheet = nm, x = as.data.frame(export_list[[nm]]))
}

wb_save(wb, file = out_file)
message("Results saved to ", out_file)

# ---------------- End of Script ---------------------------------------------
