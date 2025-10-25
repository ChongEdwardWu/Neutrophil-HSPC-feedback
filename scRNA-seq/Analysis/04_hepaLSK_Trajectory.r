### ---------------- Section 0: Preparation -----------------------------------
rm(list = ls()); gc()

# Load personal R configuration (if applicable)
suppressMessages(source("path_to_.radian_profile"))
# .libPaths()

# Set working directory
workdir <- "path_to_data"
setwd(workdir)

# Create folders if not exist
if (!file.exists(file.path(workdir, "figures"))) dir.create(file.path(workdir, "figures"))
if (!file.exists(file.path(workdir, "results"))) dir.create(file.path(workdir, "results"))

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratWrappers)
  library(SeuratObject)
  library(patchwork)
  library(ggplot2)
  library(future)
  library(SingleR)
  library(BiocParallel)
  library(openxlsx2)
  library(monocle3)
  library(plotly)
})

set.seed(123)

# Optional: for interactive visualization
# httpgd::hgd()

# Load annotated Seurat object
seu <- readRDS("results/02_Annotation_LSK.rds")


### ---------------- Section 1: Trajectory Inference --------------------------
# Inspect current cell type annotations
table(seu$CellType_l2)

# Visualize group differences on UMAP
p01 <- DimPlot(
  seu,
  reduction = "umap",
  split.by = "group",
  group.by = "CellType_l2",
  label = TRUE
) + coord_fixed(ratio = 1)
p01


### ---------------- Section 2: Convert and run Slingshot ---------------------
# Convert to SingleCellExperiment for trajectory inference
sce <- as.SingleCellExperiment(seu, assay = "RNA")

# Load slingshot
library(slingshot)
library(ggplot2)

# Define clusters for slingshot
colData(sce)$slingshot_clusters <- droplevels(colData(sce)$CellType_l2)
table(colData(sce)$slingshot_clusters)

# Run trajectory inference
sce <- slingshot(
  sce,
  reducedDim = "UMAP",
  start.clus = c("HSC"),  # specify the root cluster
  # end.clus = c("EryPro","SP_MyePro","LymPro"),  # optional end points
  clusterLabels = colData(sce)$slingshot_clusters,
  omega = TRUE,
  approx_points = 1000
)


### ---------------- Section 3: Save trajectory results -----------------------
# Save the trajectory-inferred object for downstream analyses
saveRDS(sce, file = "results/Step4_LSKsce_for_slingshot.rds")

# ---------------- End of Script ----------------------------------------------
