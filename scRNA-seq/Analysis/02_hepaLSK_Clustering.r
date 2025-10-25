### ---------------- Section 0: Preparation -----------------------------------
rm(list = ls()); gc()

# Personal profile
suppressMessages(source("path_to_.radian_profile"))
# .libPaths()

# Working directory
workdir <- "path_to_data"
setwd(workdir)
if (!file.exists(file.path(workdir, "figures"))) dir.create(file.path(workdir, "figures"))
if (!file.exists(file.path(workdir, "results"))) dir.create(file.path(workdir, "results"))

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
  # Optional blocks:
  # library(GSEABase); library(AUCell)
  # library(clusterProfiler); library(org.Mm.eg.db)
})

# httpgd::hgd()

set.seed(123)
nworkers <- 20
plan("sequential")
plan("multicore", workers = nworkers)

species <- "mm"   # "mm" or "hs"
groups_order <- c("BM", "SP", "SP_NAC")  # <- adjust if your groups differ

# Input (integration output; update path if needed)
integrated_rds <- "results/01.5_pySCENIC2seurat.rds"
stopifnot(file.exists(integrated_rds))

seu <- readRDS(integrated_rds)
message("Loaded integrated object: cells=", ncol(seu), " genes=", nrow(seu))

### ---------------- Section 1: Cluster handle (@res=0.2) ---------------------
# Use the clustering you chose during integration

# ! resolution choice
res <- 0.2
# Plot 1: group difference in clusters (UMAP split)
DimPlot(
  seu,
  reduction = "umap",  split.by = "group",
  group.by = paste0("integrated_snn_res.", res), label = TRUE
) + coord_fixed(ratio = 1)

# ! set resolution
Idents(seu) <- seu$CellType_l1 <- factor(seu@meta.data[[paste0("integrated_snn_res.", res)]])
# seu@meta.data[,grep("integrated_snn_res.",colnames(seu@meta.data),value = TRUE)] <- NULL

# Quick visualization
# DimPlot(seu, reduction="umap", group.by="CellType_l1", label=TRUE) + coord_fixed(1)

# Ensure 'group' is well-ordered factor for downstream conserved DE
if ("group" %in% colnames(seu@meta.data)) {
  seu$group <- factor(seu$group, levels = intersect(groups_order, unique(seu$group)))
}

# Fallback CC phase: if CCphase missing but Phase present (from Seurat CC scoring)
if (!"CCphase" %in% colnames(seu@meta.data) && "Phase" %in% colnames(seu@meta.data)) {
  seu$CCphase <- seu$Phase
}

### ---------------- Section 2: Automated annotations -------------------------
DefaultAssay(seu) <- "RNA"
# 3.1 SingleR + ImmGen (Bioconductor; robust for mouse hematopoiesis)
#    Uses the RNA assay expression. This is independent from integrated embeddings.
#    If you already added ImmGen in QC step2, this simply re-computes (optional but consistent).
# immgen_rdata <- "path_to_ImmGen_reference.RData"
# if (file.exists(immgen_rdata)) {
#   load(immgen_rdata)  # loads 'immgen'
#   message("ImmGen reference loaded.")
#   DefaultAssay(seu) <- "RNA"
#   seu <- NormalizeData(seu, verbose = FALSE)
#   sce <- as.SingleCellExperiment(seu)
#   pred <- SingleR(
#     test            = sce,
#     ref             = immgen,
#     labels          = immgen$label.fine,
#     assay.type.test = 1,
#     BPPARAM         = MulticoreParam(nworkers)
#   )
#   seu$CellType_immgen <- pred$pruned.labels
#   # Optional diagnostic:
#   # plotScoreHeatmap(pred)
# } else {
#   warning("ImmGen RData not found: ", immgen_rdata, " — skipping SingleR ImmGen.")
# }

# 3.2 AUCell (optional): score cluster signatures (e.g., Izzo et al. Nat Genet, 2020)
# make gene sets cell type markers
# ClMk<- read.csv(file = "HSPC cluster marker genes.csv", )
# MkClusters <- unique(ClMk$Cluster.annotation)
# ClMarkers <- list()
# for (i in 1:length(MkClusters)) {
#  ClMarkers[[MkClusters[i]]]<- as.character(ClMk$Gene.symbol[ClMk$Cluster.annotation %in% MkClusters[i]])
# }
# lengths(ClMarkers)
# save(ClMarkers, file="HSPC cluster marker genes.RData")

hspc_gs_file <- "HSPC_cluster_marker_genes.RData"  # make sure file exists in 'workdir'
if (file.exists(hspc_gs_file)) {
  suppressPackageStartupMessages({ library(GSEABase); library(AUCell) })
  load(hspc_gs_file)  # expects a list 'ClMarkers': names are set names; values are SYMBOL vectors
  all.sets <- GeneSetCollection(lapply(names(ClMarkers), function(x) GeneSet(ClMarkers[[x]], setName=x)))
  rankings <- AUCell_buildRankings(GetAssayData(seu, assay="RNA", slot="counts"),
                                   plotStats=FALSE, verbose=FALSE, splitByBlocks=TRUE,
                                   BPPARAM=MulticoreParam(nworkers))
  cell.aucs <- AUCell_calcAUC(all.sets, rankings)
  auc_mat <- t(assay(cell.aucs))
  # Best-scoring set per cell:
  seu$annot_HSPC_AUCell <- colnames(auc_mat)[max.col(auc_mat, ties.method="first")]
} else {
  message("AUCell marker file not found: ", hspc_gs_file, " — AUCell block skipped.")
}

# Combine available auto-annotations for quick inspection (kept semantics)
auto_cols <- intersect(c("CellType_immgen","annot_HSPC_AUCell","CellType_mca","CellType_pansci"), colnames(seu@meta.data))
if (length(auto_cols) > 0) {
  seu$annot_auto <- factor(apply(seu@meta.data[, auto_cols, drop=FALSE], 1, function(v) paste(na.omit(as.character(v)), collapse="; ")))
}

# View cluster annotation
cluster_annot <- tibble(
  seu$seurat_clusters,
  seu$group,
  seu$annot_auto
) %>%
  set_names("Cluster", "Source", "CellType") %>%
  group_by(Cluster, CellType) %>%
  summarise(no.cell = n()) %>%
  group_by(Cluster) %>%
  mutate(
    total.no = sum(no.cell),
    perc = 100 * no.cell / total.no
  ) %>%
  arrange(Cluster, dplyr::desc(perc)) %>%
  top_n(n = 5, wt = perc)

# View and plot cluster distribution
source_cluster <- tibble(
  seu$seurat_clusters,
  seu$group
) %>%
  set_names("Cluster", "Group") %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n()) %>%
  group_by(Group) %>%
  mutate(
    total.no = sum(no.cell),
    perc = 100 * no.cell / total.no
  ) %>%
  dplyr::select(Cluster, Group, perc)

library(ggplot2)
ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") +
  # scale_fill_manual(values=cl) +
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("group") +
  ylab("%")

# View and plot CCphase distribution
source_CCphase <- tibble(
  seu$seurat_clusters,
  seu$CCphase
) %>%
  set_names("Cluster", "CCphase") %>%
  group_by(Cluster, CCphase) %>%
  summarise(no.cell = n()) %>%
  group_by(CCphase) %>%
  mutate(
    total.no = sum(no.cell),
    perc = 100 * no.cell / total.no
  ) %>%
  dplyr::select(Cluster, CCphase, perc)

# Export final results
save_file_name <- "results/Step2_cluster_celltype.xlsx"
wb <- wb_workbook()

wb <- wb_add_worksheet(wb, sheet = "cluster_annot")
wb <- wb_add_data(wb, sheet = "cluster_annot",
                  x = as.data.frame(cluster_annot))

source_cluster_wide <- source_cluster %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = perc)
wb <- wb_add_worksheet(wb, sheet = "source_cluster")
wb <- wb_add_data(wb, sheet = "source_cluster",
                  x = as.data.frame(source_cluster_wide))

source_CCphase_wide <- source_CCphase %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = perc)
wb <- wb_add_worksheet(wb, sheet = "source_CCphase")
wb <- wb_add_data(wb, sheet = "source_CCphase",
                  x = as.data.frame(source_CCphase_wide))

wb_save(wb, save_file_name, overwrite = TRUE)

### ---------------- Section 3: Feature genes dot plot -------------------------
Idents(seu) <- "seurat_clusters"
# Compact mouse hematopoiesis/immune lineage panel (from Kucinski, Iwo et al., Cell Stem Cell, Volume 31, Issue 2, 244 - 259.e)
lineage_panel_mm <- list(
  HSC                         = c("Procr", "Mecom"),
  `Megakaryocyte (progenitor)` = c("Pf4", "Vwf", "Itga2b", "Gp9", "Itgb3"),
  `Erythroid (progenitor)`     = c("Klf1", "Gata1", "Hba-a1"),
  `Neutrophil (progenitor)`    = c("Elane", "Prtn3", "Ms4a3"),
  `Basophil (progenitor)`      = c("Prss34", "Mcpt8"),
  `Eosinophil (progenitor)`    = c("Prg2", "Prg3", "Epx"),
  `Mast cell (progenitor)`     = c("Kit", "Cma1", "Gzmb"),
  `Monocyte/DC progenitor`     = c("Csf1r", "Mpo"),
  `Lymphoid progenitor`        = c("Il7r", "Dntt"),
  `T cell (progenitor)`        = c("Cd3e", "Cd8a", "Bcl11b"),
  `B cell (progenitor)`        = c("Ly6d", "Vpreb3", "Cd79a"),
  ILC                          = c("Gata3", "Id2", "Il7r"),
  Monocyte                     = c("Cd14", "Ctsg", "Ly6c1", "Ly6c2"),
  pDC                          = c("Irf8", "Siglech", "Ly6d"),
  `Ifn-activated cells`        = c("Ifitm3", "Irf7", "Isg15"),
  `Complement expressing`      = c("C1qa", "C1qb"),
  DCs                          = c("H2-Aa", "Xcr1", "Itgax")
)

# keep only genes present in object
gene_panel <- unique(unlist(lineage_panel_mm))
gene_panel <- intersect(gene_panel, rownames(seu))

DotPlot(
  seu,
  group.by = "seurat_clusters",
  assay   = "SCT",
  features = gene_panel
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title   = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        legend.direction = "vertical",
        legend.position  = "bottom") +
  scale_color_gradientn(colours = c('white',"#fde3d8","#ed684e","#d21e20","#981917"))

### ---------------- Section 4: Cluster merging / manual naming ---------------
seu$CellType <- NULL
# Annotate clusters
seu$CellType_l1 <- factor(dplyr::recode(
  seu$seurat_clusters,
  "1"  = "HSC",
  "2"  = "HSC",
  "3"  = "MyePro",
  "4"  = "LymPro",
  "5"  = "MyePro",
  "6"  = "LMPP",
  "7"  = "EryPro"
))
seu$CellType_l1 <- factor(seu$CellType_l1,
                       levels = c(
                         "HSC","MyePro","EryPro","LMPP","LymPro"
                       ))

table(seu$seurat_clusters, seu$CellType_l1)
table(seu$group, seu$CellType_l1)

# Quick UMAP by CellType
DimPlot(seu, reduction="umap", group.by="CellType_l1", label=TRUE) + coord_fixed(1)
DimPlot(seu, reduction="umap", group.by="CellType_l1", split.by = "group",label=TRUE) + coord_fixed(1)

# Annotate clusters
seu$CellType_l2 <- factor(dplyr::recode(
  seu$seurat_clusters,
  "1"  = "HSC",
  "2"  = "HSC_cycling",
  "3"  = "SP_MyePro",
  "4"  = "LymPro",
  "5"  = "BM_MyePro",
  "6"  = "LMPP",
  "7"  = "EryPro"
))

seu$CellType_l2 <- factor(seu$CellType_l2,
                       levels = c(
                         "HSC","HSC_cycling","BM_MyePro","SP_MyePro","EryPro","LMPP","LymPro"
                       ))

table(seu$seurat_clusters, seu$CellType_l2)
table(seu$group, seu$CellType_l2)

# Quick UMAP by CellType
DimPlot(seu, reduction="umap", group.by="CellType_l2", label=TRUE) + coord_fixed(1)

### ---------------- Section 5: Marker genes per cluster ----------------------
# Use RNA assay; exclude "bad" features to avoid artifacts (TCR/BCR/MT/HB/etc.)
Idents(seu) <- "CellType_l2"
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)

rna_genes <- rownames(seu[["RNA"]])

# Build mouse blacklist (programmatic). Extend as needed.

# Histone
hist_genes   <- grep("^Hist", rna_genes, ignore.case = TRUE, value = TRUE)
# Hemoglobin
hb_genes     <- grep("^Hb[ab]-|^HB(?!P)", rna_genes, perl = TRUE, value = TRUE)
# Mitochondrial
mt_genes     <- grep("^mt-|^MT-", rna_genes, ignore.case = TRUE, value = TRUE)
# Ribosomal
rps_genes    <- grep("^Rp[sl]|^RP[SL]", rna_genes, ignore.case = TRUE, value = TRUE)
# Riken (mouse)
rik_genes     <- grep("Rik$", rna_genes, value = TRUE)
# Alu elements
if (species == "hs") {
  alu_genes <- grep("^AL[0-9]+", rna_genes, value = TRUE)
} else {
  alu_genes <- character(0)
}
# Pseudogenes
pseudo_genes  <- grep("-(ps|ps[0-9]+)$", rna_genes, ignore.case = TRUE, value = TRUE)
# microRNAs
mir_genes    <- grep("^Mir", rna_genes, ignore.case = TRUE, value = TRUE)
# Gencode-style “Gm” (mouse models)
gencode_genes <- grep("^Gm[0-9]+", rna_genes, value = TRUE)

# Combine
bad_features <- unique(c(
  # intersect(cc_genes, rna_genes),
  # hist_genes, hb_genes, rps_genes, mt_genes,
  rik_genes, gencode_genes, pseudo_genes, alu_genes
))

features <- setdiff(rna_genes, bad_features)

# Markers by cluster
Idents(seu) <- seu$CellType_l2
seu <- PrepSCTFindMarkers(seu)
seu_markers <- FindAllMarkers(
  seu,
  only.pos       = FALSE,
  min.pct        = 0.1,
  logfc.threshold= 0,
  features       = features,
  assay = "RNA"
)

seu_markers <- seu_markers %>% group_by(cluster) %>% arrange(p_val_adj, .by_group = TRUE)
# top10_by_cluster <- seu_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)

# View and plot cluster distribution
source_cluster <- tibble(
  seu$CellType_l2,
  seu$group
) %>%
  set_names("Cluster", "Group") %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n()) %>%
  group_by(Group) %>%
  mutate(
    total.no = sum(no.cell),
    perc = 100 * no.cell / total.no
  ) %>%
  dplyr::select(Cluster, Group, perc)

# for exporting results
library(openxlsx2)

# Export final results
save_file_name <- "results/Step2_cluster_celltype.xlsx"
wb <- wb_workbook()

wb <- wb_add_worksheet(wb, sheet = "source_cluster")
wb <- wb_add_data(wb, sheet = "source_cluster",
                  x = as.data.frame(source_cluster))

# Export marker tables
wb <- wb_add_worksheet(wb, "Clus_markers")
wb <- wb_add_data(wb, "Clus_markers", as.data.frame(seu_markers))
wb_save(wb, save_file_name, overwrite = TRUE)

# Regulon-AUC Markers by cluster
seu_AUC <- FindAllMarkers(
  seu,
  assay = "AUC",
  only.pos       = FALSE,
  min.pct        = 0.1,
  logfc.threshold= 0
)
seu_AUC <- seu_AUC %>% group_by(cluster) %>% arrange(p_val_adj, .by_group = TRUE)
# Export marker tables
wb <- wb_add_worksheet(wb, "Clus_Regulon_AUC")
wb <- wb_add_data(wb, "Clus_Regulon_AUC", as.data.frame(seu_AUC))
wb_save(wb, save_file_name, overwrite = TRUE)

# Regulon-Bin Markers by cluster (positive markers only)
seu_Bin <- FindAllMarkers(
  seu,
  assay = "Bin",
  only.pos       = FALSE,
  min.pct        = 0.1,
  logfc.threshold= 0
)
seu_Bin <- seu_Bin %>% group_by(cluster) %>% arrange(p_val_adj, .by_group = TRUE)
# Export marker tables
wb <- wb_add_worksheet(wb, "Clus_Regulon_Bin")
wb <- wb_add_data(wb, "Clus_Regulon_Bin", as.data.frame(seu_Bin))
wb_save(wb, save_file_name, overwrite = TRUE)

### ---------------- Section 6: Cluster feature genes selection and plot ------
DefaultAssay(seu) <- "SCT"
feature_genes <- list(
    `commonLSK` = c("Kit", "Ly6a","Flt3"),
    `LTHSC` = c("Hlf", "Mecom", "Mpl", "Procr", "Fgd5"),
    `HSC_cycling` = c("Mpo", "Ctsg", "Mki67", "Top2a", "Cdk1"),
    `BM_MyePro` = c("Lyz2", "Prtn3","Apoe", "Selenop","Ctsd", "Ctsb"),
    `SP_MyePro` = c( "Clec12a","Ly86","Irf8","Batf3"),
    `EryPro` = c("Klf1", "Gata1", "Car1"),
    `LMPP` = c("Dntt", "Cd33","Blnk"),
    `LymPro` = c( "Il7r", "Tcf7", "Il2ra","Bcl11b","Themis")
)

gene_panel <- unique(unlist(feature_genes))
gene_panel <- intersect(gene_panel, rownames(seu))

DotPlot(
  seu,
  group.by = "CellType_l2",
  assay   = "SCT",
  features = gene_panel
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title   = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        legend.direction = "vertical",
        legend.position  = "bottom") +
  scale_color_gradientn(colours = c('white',"#fde3d8","#ed684e","#d21e20","#981917"))

### ---------------- Section 9: Cluster feature regulons selection and plot ---
feature_TF <- list(
  LTHSC      = c("Meis1-pos", "Sox6-pos", "Etv6-pos", "Fli1-pos", "Runx1-pos"),
  HSC_cycling  = c("E2f1-pos", "E2f2-pos","Brca1-pos"),
  BM_MyePro  = c("Hoxa9-pos", "Myc-pos","Trp53-pos"),
  SP_MyePro  = c("Spi1-pos", "Irf8-pos", "Atf3-pos",  "Relb-pos", "Foxo1-pos"),
  EryPro     = c("Gata1-pos", "Klf1-pos", "Nfia-pos"),
  LMPP       = c("Ikzf1-pos", "Tcf4-pos", "Lef1-pos"),
  LymPro     = c("Nfil3-pos", "Gata3-pos", "Ets1-pos")
)

tf_panel <- unique(unlist(feature_TF))
tf_panel <- intersect(tf_panel, rownames(seu[["Bin"]]))

DotPlot(
  seu,
  group.by = "CellType_l2",
  assay   = "Bin",
  features = tf_panel
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title   = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        legend.direction = "vertical",
        legend.position  = "bottom") +
  scale_color_gradientn(colours = c('white',"#fde3d8","#ed684e","#d21e20","#981917"))

# Save annotated object
saveRDS(seu, file = "results/02_Annotation_LSK.rds")
# seu <- readRDS("results/02_Annotation_LSK.rds")

# Optional: save session info
writeLines(capture.output(sessionInfo()), "session_info.txt")

# End of script
