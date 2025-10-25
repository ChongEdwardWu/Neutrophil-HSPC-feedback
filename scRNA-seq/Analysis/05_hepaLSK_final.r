### ---------------- Section 0: Preparation -----------------------------------
rm(list = ls()); gc()

# Load personal R configuration (if applicable)
suppressMessages(source("path_to_.radian_profile"))
# .libPaths()

# Set working directory
workdir <- "path_to_data"
setwd(workdir)

# Create directories for figures and results
if (!file.exists(file.path(workdir, "figures"))) dir.create(file.path(workdir, "figures"))
if (!file.exists(file.path(workdir, "results"))) dir.create(file.path(workdir, "results"))

# Load required packages
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
  library(ggrepel)
  # Optional blocks:
  # library(GSEABase); library(AUCell)
  # library(clusterProfiler); library(org.Mm.eg.db)
})

set.seed(123)
nworkers <- 8

# Species configuration
species <- "mm"  # or "mm" for mouse
msig_species <- if (species == "hs") "Homo sapiens" else "Mus musculus"

# Create directory for final figures
if (!file.exists(file.path(workdir, "figures/final"))) {
  dir.create(file.path(workdir, "figures/final"), recursive = TRUE)
}

# Load annotated Seurat object
seu <- readRDS("results/02_Annotation_LSK.rds")

# Optional: for interactive plotting
# httpgd::hgd()


### ---------------- Section 1: UMAP and cluster visualization -----------------
table(seu$CellType_l2)
Idents(seu) <- seu$CellType_l2

# Define color palette
colorpanel <- c(
  "HSC" = "#575bb7",
  "HSC_cycling" = "#8c51c1",
  "BM_MyePro" = "#5cadff",
  "SP_MyePro" = "#e13d2d",
  "EryPro" = "#54a546",
  "LMPP" = "#66cccc",
  "LymPro" = "#fa9645"
)

# UMAP (save explicitly)
p_umap <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "CellType_l2",
  label = FALSE
) +
  coord_fixed(ratio = 1.5) +
  scale_color_manual(values = colorpanel)

ggsave("figures/final/fig1_hepaLSK_UMAP_L2.png", plot = p_umap, width = 150, height = 150, units = "mm", dpi = 300)
ggsave("figures/final/fig1_hepaLSK_UMAP_L2.pdf", plot = p_umap, width = 150, height = 150, units = "mm", bg = "transparent")


### ---------------- Section 2: Marker and regulon plots -----------------------
DefaultAssay(seu) <- "SCT"

# Cell-type markers
feature_genes <- list(
  `commonLSK` = c("Kit", "Ly6a", "Flt3"),
  `LTHSC` = c("Hlf", "Mecom", "Mpl", "Procr", "Fgd5"),
  `HSC_cycling` = c("Mpo", "Ctsg", "Mki67", "Top2a", "Cdk1"),
  `BM_MyePro` = c("Lyz2", "Prtn3", "Apoe", "Selenop", "Ctsd", "Ctsb"),
  `SP_MyePro` = c("Clec12a", "Ly86", "Irf8", "Batf3"),
  `EryPro` = c("Klf1", "Gata1", "Car1"),
  `LMPP` = c("Dntt", "Cd33", "Blnk"),
  `LymPro` = c("Il7r", "Tcf7", "Il2ra", "Bcl11b", "Themis")
)
gene_panel <- intersect(unique(unlist(feature_genes)), rownames(seu))

p_feat <- DotPlot(seu, group.by = "CellType_l2", assay = "SCT", features = rev(gene_panel)) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "bottom"
  ) +
  scale_color_gradientn(colours = c("white", "#fde3d8", "#ed684e", "#d21e20", "#981917"))

ggsave("figures/final/fig2_hepaLSK_CellTypeFeat_dotplot.png", plot = p_feat, width = 100, height = 250, units = "mm", dpi = 300)
ggsave("figures/final/fig2_hepaLSK_CellTypeFeat_dotplot.pdf", plot = p_feat, width = 100, height = 250, units = "mm", bg = "transparent")

# Regulon features (Bin assay)
feature_TF <- list(
  LTHSC = c("Meis1-pos", "Sox6-pos", "Etv6-pos", "Fli1-pos", "Runx1-pos"),
  HSC_cycling = c("E2f1-pos", "E2f2-pos", "Brca1-pos"),
  BM_MyePro = c("Hoxa9-pos", "Myc-pos", "Trp53-pos"),
  SP_MyePro = c("Spi1-pos", "Irf8-pos", "Atf3-pos", "Relb-pos", "Foxo1-pos"),
  EryPro = c("Gata1-pos", "Klf1-pos", "Nfia-pos"),
  LMPP = c("Ikzf1-pos", "Tcf4-pos", "Lef1-pos"),
  LymPro = c("Nfil3-pos", "Gata3-pos", "Ets1-pos")
)
tf_panel <- intersect(unique(unlist(feature_TF)), rownames(seu[["Bin"]]))

p_reg <- DotPlot(seu, group.by = "CellType_l2", assay = "Bin", features = rev(tf_panel)) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "bottom"
  ) +
  scale_color_gradientn(colours = c("white", "#fde3d8", "#ed684e", "#d21e20", "#981917"))

ggsave("figures/final/fig3_hepaLSK_CellTypeRegulon_dotplot.png", plot = p_reg, width = 100, height = 250, units = "mm", dpi = 300)
ggsave("figures/final/fig3_hepaLSK_CellTypeRegulon_dotplot.pdf", plot = p_reg, width = 100, height = 250, units = "mm", bg = "transparent")


### ---------------- Section 3: Cluster distribution --------------------------
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
  )

p_abund <- ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") +
  scale_fill_manual(values = colorpanel) +
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("Group") +
  ylab("%")

ggsave("figures/final/fig4_hepaLSK_cluster_abundance.png", plot = p_abund, width = 100, height = 100, units = "mm", dpi = 300)
ggsave("figures/final/fig4_hepaLSK_cluster_abundance.pdf", plot = p_abund, width = 100, height = 100, units = "mm", bg = "transparent")


### ---------------- Section 4: UMAP by group ---------------------------------
p_split <- DimPlot(
  seu,
  reduction = "umap",
  split.by = "group",
  group.by = "CellType_l2",
  label = FALSE
) +
  coord_fixed(ratio = 1.5) +
  scale_color_manual(values = colorpanel)

ggsave("figures/final/fig5_hepaLSK_UMAP_split.png", plot = p_split, width = 250, height = 150, units = "mm", dpi = 300)
ggsave("figures/final/fig5_hepaLSK_UMAP_split.pdf", plot = p_split, width = 250, height = 150, units = "mm", bg = "transparent")


### ---------------- Section 5: Density visualization -------------------------
library(viridis)

theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white", fill = "black"),
      panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black")
    )
}

coord <- Embeddings(seu, reduction = "umap")[, 1:2]
meta <- seu@meta.data %>% tibble::rownames_to_column("ID") %>%
  left_join(data.frame(ID = rownames(coord), coord), by = "ID")

density_plot <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  geom_point(color = "white", size = 0.01) +
  scale_fill_viridis(option = "magma") +
  coord_fixed(ratio = 1.5) +
  facet_wrap(~group, nrow = 1) +
  theme_black()

ggsave("figures/final/fig6_hepaLSK_density.png", plot = density_plot, width = 400, height = 200, units = "mm", dpi = 300)
ggsave("figures/final/fig6_hepaLSK_density.pdf", plot = density_plot, width = 400, height = 200, units = "mm", bg = "transparent")


### ---------------- Section 6: MyePro DE regulon plots (Figures 7–9) ---------
# Load cluster-level regulon markers from Step-2 export
clus_regulons <- openxlsx2::wb_read(
  file  = "results/Step2_cluster_celltype.xlsx",
  sheet = "Clus_Regulon_Bin",
  col_names = TRUE
)

SP_MyePro_regulons_df <- clus_regulons |>
  dplyr::filter(cluster == "SP_MyePro") |>
  dplyr::filter(p_val_adj < 0.05) |>
  dplyr::mutate(pct_diff = pct.1 - pct.2) |>
  dplyr::filter(pct_diff > 0.1)

SP_MyePro_regulons_gene <- SP_MyePro_regulons_df$gene

# Load DE regulons (Bin assay) from Step-3 exports
# SP vs BM
MyePro_SPvsBM_DEregulons <- openxlsx2::wb_read(
  file  = "results/Step3_Group_Comparisons.xlsx",
  sheet = "MyePro_SPvsBM_regulons_Bin",
  col_names = TRUE
)
MyePro_SPvsBM_DEregulons_df <- MyePro_SPvsBM_DEregulons |>
  dplyr::filter(p_val_adj < 0.05) |>
  dplyr::mutate(pct_diff = pct.1 - pct.2) |>
  dplyr::filter(pct_diff > 0.1)
MyePro_SPvsBM_DEregulons_gene <- MyePro_SPvsBM_DEregulons_df$Regulon

# SP_NAC vs SP
MyePro_SPNACvsSP_DEregulons <- openxlsx2::wb_read(
  file  = "results/Step3_Group_Comparisons.xlsx",
  sheet = "MyePro_SPNACvsSP_regulons_Bin",
  col_names = TRUE
)
MyePro_SPNACvsSP_DEregulons_df <- MyePro_SPNACvsSP_DEregulons |>
  dplyr::filter(p_val_adj < 0.05) |>
  dplyr::mutate(pct_diff = pct.1 - pct.2) |>
  dplyr::filter(pct_diff < -0.1)
MyePro_SPNACvsSP_DEregulons_gene <- MyePro_SPNACvsSP_DEregulons_df$Regulon

# Intersections & Venn: SP_MyePro sets
sets <- list(
  regulons_cluster_pos = unique(na.omit(as.character(SP_MyePro_regulons_gene))),
  SP_vs_BM_Up          = unique(na.omit(as.character(MyePro_SPvsBM_DEregulons_gene))),
  SPNAC_vs_SP_Down     = unique(na.omit(as.character(MyePro_SPNACvsSP_DEregulons_gene)))
)

cat("Set sizes:\n")
print(vapply(sets, length, integer(1)))

# Triple and pairwise intersections
intersect_all       <- Reduce(intersect, sets)  # triple
intersect_Mk_SpBm   <- intersect(sets$regulons_cluster_pos, sets$SP_vs_BM_Up)
intersect_Mk_SpNAC  <- intersect(sets$regulons_cluster_pos, sets$SPNAC_vs_SP_Down)
intersect_SpBm_SpN  <- intersect(sets$SP_vs_BM_Up, sets$SPNAC_vs_SP_Down)

# Export lists (no date in filename)
venn_lists_xlsx <- file.path("results", "Final_SP_MyePro_Regulon_Venn.xlsx")
openxlsx2::write_xlsx(
  x = list(
    "regulons_cluster_pos_all"      = data.frame(Gene = sets$regulons_cluster_pos),
    "SP_vs_BM_Up_all"               = data.frame(Gene = sets$SP_vs_BM_Up),
    "SPNAC_vs_SP_Down_all"          = data.frame(Gene = sets$SPNAC_vs_SP_Down),
    "Overlap_Triple"                = data.frame(Gene = intersect_all),
    "Overlap_cluster_x_SPvsBMUp"    = data.frame(Gene = intersect_Mk_SpBm),
    "Overlap_cluster_x_SPNACdown"   = data.frame(Gene = intersect_Mk_SpNAC),
    "Overlap_SPvsBMUp_x_SPNACdown"  = data.frame(Gene = intersect_SpBm_SpN)
  ),
  file = venn_lists_xlsx
)
message("Gene lists written to: ", venn_lists_xlsx)

# 3-set Venn diagram (save explicitly)
suppressPackageStartupMessages(library(ggVennDiagram))
p_venn <- ggVennDiagram(
  sets,
  label_alpha = 0,
  label = "count",
  edge_size = 0.7
) +
  scale_fill_gradient(low = "#F8F9FA", high = "#6C9BD2") +
  theme(legend.position = "none")

ggsave("figures/final/fig7_SP_MyePro_regulon_VennDiagram.png", plot = p_venn,
       width = 100, height = 100, units = "mm", dpi = 300)
ggsave("figures/final/fig7_SP_MyePro_regulon_VennDiagram.pdf", plot = p_venn,
       width = 100, height = 100, units = "mm", bg = "transparent")

# Class mapping for regulons (for colored scatter)
class_tbl <- tibble::tribble(
  ~Regulon,    ~Class,
  "Spi1-pos",  "Myeloid Differentiation",
  "Stat1-pos", "Myeloid Activation",
  "Irf5-pos",  "Myeloid Activation",
  "Irf8-pos",  "Myeloid Differentiation",
  "Atf3-pos",  "Oxid and integrated stress response",
  "Foxo1-pos", "Oxid and integrated stress response",
  "Bcl6-pos",  "Myeloid Differentiation",
  "Etv3-pos",  "Myeloid Differentiation",
  "Fos-pos",   "Oxid and integrated stress response",
  "Junb-pos",  "Oxid and integrated stress response",
  "Relb-pos",  "Myeloid Activation",
  "Elf4-pos",  "Myeloid Differentiation",
  "Nr2f6-pos", "Myeloid Differentiation",
  "Hltf-pos",  "Oxid and integrated stress response",
  "Zbtb7a-pos","Myeloid Differentiation",
  "Nfkb2-pos", "Myeloid Differentiation",
  "Atf6-pos",  "Oxid and integrated stress response",
  "Nfatc2-pos","Myeloid Activation",
  "Hes7-pos",  "Others",
  "Irf7-pos",  "Myeloid Activation",
  "Usf2-pos",  "Oxid and integrated stress response",
  "Rel-pos",   "Myeloid Activation",
  "Nfkb1-pos", "Myeloid Activation",
  "Zfp78-pos", "Others",
  "Irf1-pos",  "Myeloid Activation",
  "Irf2-pos",  "Myeloid Activation"
)
class_levels <- c(
  "Oxid and integrated stress response",
  "Myeloid Differentiation",
  "Myeloid Activation",
  "Others"
)

# Figure 8: SP vs BM dot plot (MyePro regulon activity)
MyePro_SPvsBM_DEregulons_df_in <- MyePro_SPvsBM_DEregulons_df |>
  dplyr::filter(Regulon %in% intersect_all)

df_plot <- MyePro_SPvsBM_DEregulons_df_in %>%
  dplyr::mutate(
    Regulon_clean = sub("[-_]?pos$", "", Regulon),
    neglog10FDR   = pmin(-log10(p_val_adj), 300)
  ) %>%
  dplyr::left_join(class_tbl, by = "Regulon") %>%
  dplyr::mutate(
    Class = dplyr::coalesce(Class, "Others"),
    Class = factor(Class, levels = class_levels)
  )

suppressPackageStartupMessages(library(ggrepel))
p_mye_spbm <- ggplot(df_plot, aes(x = pct.2, y = pct.1)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +
  geom_point(aes(size = neglog10FDR, color = Class), alpha = 0.9) +
  geom_text_repel(aes(label = Regulon_clean), size = 3, min.segment.length = 0,
                  box.padding = 0.25, point.padding = 0.25, max.overlaps = Inf) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0.02)) +
  coord_equal() +
  scale_size_continuous(range = c(2, 5), limits = c(0, 300), breaks = c(100, 200, 300),
                        name = expression(-log[10]("FDR"))) +
  scale_color_manual(
    values = c(
      "Oxid and integrated stress response" = "#e13d2d",
      "Myeloid Differentiation"             = "#5cadff",
      "Myeloid Activation"                  = "#54a546",
      "Others"                              = "black"
    ),
    name = "Class"
  ) +
  labs(
    x = "Active fraction in BM MyePro",
    y = "Active fraction in SP MyePro",
    title = "SP vs BM (MyePro) regulon activity"
  ) +
  theme_bw(base_size = 6) +
  theme(
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", hjust = 0)
  )

ggsave("figures/final/fig8_MyePro_regulon_SPvsBM.png", plot = p_mye_spbm, width = 120, height = 100, units = "mm", dpi = 300)
ggsave("figures/final/fig8_MyePro_regulon_SPvsBM.pdf", plot = p_mye_spbm, width = 120, height = 100, units = "mm", bg = "transparent")

# Figure 9: SP_NAC vs SP dot plot (MyePro regulon activity)
MyePro_SPNACvsSP_DEregulons_df_in <- MyePro_SPNACvsSP_DEregulons_df |>
  dplyr::filter(Regulon %in% intersect_all)

df_plot2 <- MyePro_SPNACvsSP_DEregulons_df_in %>%
  dplyr::mutate(
    Regulon_clean = sub("[-_]?pos$", "", Regulon),
    neglog10FDR   = pmin(-log10(p_val_adj), 300)
  ) %>%
  dplyr::left_join(class_tbl, by = "Regulon") %>%
  dplyr::mutate(
    Class = dplyr::coalesce(Class, "Others"),
    Class = factor(Class, levels = class_levels)
  )

p_mye_spnac <- ggplot(df_plot2, aes(x = pct.2, y = pct.1)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +
  geom_point(aes(size = neglog10FDR, color = Class), alpha = 0.9) +
  geom_text_repel(aes(label = Regulon_clean), size = 3, min.segment.length = 0,
                  box.padding = 0.25, point.padding = 0.25, max.overlaps = Inf) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0.02)) +
  coord_equal() +
  scale_size_continuous(range = c(2, 5), limits = c(0, 300), breaks = c(100, 200, 300),
                        name = expression(-log[10]("FDR"))) +
  scale_color_manual(
    values = c(
      "Oxid and integrated stress response" = "#e13d2d",
      "Myeloid Differentiation"             = "#5cadff",
      "Myeloid Activation"                  = "#54a546",
      "Others"                              = "black"
    ),
    name = "Class"
  ) +
  labs(
    x = "Active fraction in SP MyePro",
    y = "Active fraction in NAC-treated SP MyePro",
    title = "SP_NAC vs SP (MyePro) regulon activity"
  ) +
  theme_bw(base_size = 6) +
  theme(
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", hjust = 0)
  )

ggsave("figures/final/fig9_MyePro_regulon_SPNACvsSP.png", plot = p_mye_spnac, width = 120, height = 100, units = "mm", dpi = 300)
ggsave("figures/final/fig9_MyePro_regulon_SPNACvsSP.pdf", plot = p_mye_spnac, width = 120, height = 100, units = "mm", bg = "transparent")


### ---------------- Section 7: Trajectory plot (Figure 10) -------------------
suppressPackageStartupMessages({
  library(slingshot)
  library(scater)
})

# Load Slingshot-ready SCE saved in Step-4 (no date suffix)
sce <- readRDS("results/Step4_LSKsce_for_slingshot.rds")

# Branch assignment and pseudotime summary
curveAssignments <- slingBranchID(sce)
sce$curveAssignments <- curveAssignments

pseudo.paths <- slingPseudotime(sce)
shared.pseudo <- rowMeans(pseudo.paths, na.rm = TRUE)

# UMAP with shared pseudotime; overlay principal curves
p_traj <- plotUMAP(sce, colour_by = I(shared.pseudo))
embedded <- embedCurves(sce, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  dfc <- data.frame(path$s[path$ord, ])
  p_traj <- p_traj + geom_path(data = dfc, aes(x = umap_1, y = umap_2), size = 1.2) +
    coord_fixed(ratio = 1.5)
}

ggsave("figures/final/fig10_hepaLSK_trajectory.png", plot = p_traj, width = 150, height = 150, units = "mm", dpi = 300)
ggsave("figures/final/fig10_hepaLSK_trajectory.pdf", plot = p_traj, width = 150, height = 150, units = "mm", bg = "transparent")


### ---------------- Section 8: Lineage-1 signatures vs pseudotime (Fig. 11–12)
suppressPackageStartupMessages({
  library(UCell)
  library(msigdbr)
  library(pheatmap)
})

# 0) Collect MSigDB gene sets (auto species via msig_species)
target_sets <- c(
  "REACTOME_SELECTIVE_AUTOPHAGY",
  "REACTOME_AUTOPHAGY",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_PEROXISOME",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_HEME_METABOLISM",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",
  "GOBP_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_REGULATION_OF_AUTOPHAGY",
  "GOBP_REACTIVE_NITROGEN_SPECIES_METABOLIC_PROCESS",
  "GOBP_PYRUVATE_METABOLIC_PROCESS",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
  "GOBP_CELLULAR_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "GOBP_AUTOPHAGY_OF_MITOCHONDRION",
  "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT"
)

m_df_reactome <- msigdbr(species = msig_species, category = "C2", subcategory = "CP:REACTOME")
m_df_go       <- msigdbr(species = msig_species, category = "C5", subcategory = "GO:BP")
m_df_h        <- msigdbr(species = msig_species, category = "H")
msig_all <- dplyr::bind_rows(m_df_reactome, m_df_go, m_df_h) |>
  dplyr::distinct(gs_name, gene_symbol, .keep_all = TRUE)
gs_tbl <- msig_all |> dplyr::filter(gs_name %in% target_sets)

missing_sets <- setdiff(target_sets, unique(gs_tbl$gs_name))
if (length(missing_sets)) {
  warning("These target sets were not found in msigdbr for the selected species: ",
          paste(missing_sets, collapse = ", "))
}

gs_list_raw <- split(gs_tbl$gene_symbol, gs_tbl$gs_name)
gs_list_raw <- lapply(gs_list_raw, function(v) unique(na.omit(v)))

# Intersect with available genes and keep reasonable sizes
min_genes <- 50
max_genes <- 500
gs_list <- lapply(gs_list_raw, function(v) intersect(v, rownames(seu)))
gs_list <- gs_list[vapply(gs_list, length, integer(1)) >= min_genes]
gs_list <- gs_list[vapply(gs_list, length, integer(1)) <= max_genes]
message("Gene sets retained for UCell: ", length(gs_list))

# 1) UCell scores (SCT assay)
seu <- UCell::AddModuleScore_UCell(seu, features = gs_list, assay = "SCT", ncores = 20)

# 2) Restrict to lineage-1 cells (high weight on lineage-1)
pt_mat <- slingPseudotime(sce)       # cells x lineages
w_mat  <- slingCurveWeights(sce)
pt1 <- pt_mat[, 1]; w1 <- w_mat[, 1]
keep_cells <- names(pt1)[!is.na(pt1) & w1 > 0.7]
pt1_vec    <- pt1[keep_cells]
stopifnot(all(keep_cells %in% colnames(seu)))

# 3) Correlate UCell scores with pseudotime (by group)
cor_by_group <- function(score_mat, pt, groups, min_n = 30) {
  sigs <- colnames(score_mat)
  out  <- list()
  for (g in levels(groups)) {
    idx <- names(pt)[groups == g]
    idx <- idx[!is.na(pt[idx])]
    if (length(idx) < min_n) next
    sc_g <- score_mat[idx, , drop = FALSE]
    pt_g <- pt[idx]
    rho  <- apply(sc_g, 2, function(x) suppressWarnings(cor(x, pt_g, method = "spearman")))
    pvl  <- apply(sc_g, 2, function(x) suppressWarnings(cor.test(x, pt_g, method = "spearman")$p.value))
    out[[g]] <- tibble(Signature = sigs, rho = as.numeric(rho), pval = as.numeric(pvl),
                       Group = g, N = length(idx))
  }
  bind_rows(out) |>
    dplyr::mutate(padj = p.adjust(pval, "BH")) |>
    dplyr::arrange(Group, dplyr::desc(abs(rho)))
}

sig_cols <- grep("_UCell$", colnames(seu@meta.data), value = TRUE)
score_mat_l1 <- as.matrix(seu@meta.data[keep_cells, sig_cols, drop = FALSE])
colnames(score_mat_l1) <- gsub("_UCell$", "", colnames(score_mat_l1))

stopifnot("group" %in% colnames(seu@meta.data))
grp <- factor(seu@meta.data[keep_cells, "group"], levels = c("BM", "SP", "SP_NAC"))
names(grp) <- keep_cells

cor_grp <- cor_by_group(score_mat_l1, pt = pt1_vec, groups = grp, min_n = 30)
cor_all <- {
  rho <- apply(score_mat_l1, 2, function(x) cor(x, pt1_vec, method = "spearman"))
  pvl <- apply(score_mat_l1, 2, function(x) cor.test(x, pt1_vec, method = "spearman")$p.value)
  tibble(Signature = colnames(score_mat_l1),
         rho = as.numeric(rho), pval = as.numeric(pvl),
         padj = p.adjust(pvl, "BH"), Group = "All", N = length(pt1_vec)) |>
    dplyr::arrange(dplyr::desc(abs(rho)))
}

# 4) Export correlations (no date suffix)
openxlsx2::write_xlsx(
  list("Cor_All" = cor_all, "Cor_byGroup" = cor_grp),
  file = "results/final_Lineage1_UCell_vsPseudotime.xlsx"
)
message("UCell ~ pseudotime correlation written: results/final_Lineage1_UCell_vsPseudotime.xlsx")

# 5a) Heatmap of rho (signatures × groups)
rho_mat <- cor_grp |>
  dplyr::select(Signature, Group, rho) |>
  tidyr::pivot_wider(names_from = Group, values_from = rho) |>
  as.data.frame()
rownames(rho_mat) <- rho_mat$Signature; rho_mat$Signature <- NULL
rho_mat <- as.matrix(rho_mat)

pheatmap(
  rho_mat, cluster_rows = TRUE, cluster_cols = FALSE,
  color = colorRampPalette(c("#1b9e77","white","#d95f02"))(100),
  display_numbers = FALSE, border_color = NA,
  main = "Spearman rho (UCell vs pseudotime, lineage-1)",
  filename = "figures/final/fig11_Lineage1_UCell_rho_heatmap.png",
  width = 250/25.4, height = 200/25.4
)
pheatmap(
  rho_mat, cluster_rows = TRUE, cluster_cols = FALSE,
  color = colorRampPalette(c("#1b9e77","white","#d95f02"))(100),
  display_numbers = FALSE, border_color = NA,
  main = "Spearman rho (UCell vs pseudotime, lineage-1)",
  filename = "figures/final/fig11_Lineage1_UCell_rho_heatmap.pdf",
  width = 250/25.4, height = 200/25.4
)

# 5b) Per-signature trends across pseudotime (Signature × Group facets)
common_cells <- Reduce(intersect, list(rownames(score_mat_l1), names(pt1_vec), names(grp)))
if (length(common_cells) < nrow(score_mat_l1)) {
  message("Harmonizing cells: ", length(common_cells), " in common")
}
score_mat_l1 <- score_mat_l1[common_cells, , drop = FALSE]
pt1_vec      <- pt1_vec[common_cells]
grp          <- grp[common_cells]

df_long <- as.data.frame(score_mat_l1) |>
  tibble::rownames_to_column("Cell") |>
  tidyr::pivot_longer(cols = -Cell, names_to = "Signature", values_to = "score") |>
  dplyr::mutate(
    Cell      = as.character(Cell),
    Signature = as.character(Signature),
    pt        = unname(pt1_vec[Cell]),
    group     = unname(grp[Cell])
  ) |>
  dplyr::filter(!is.na(pt), !is.na(group))

sig_ord <- sort(unique(df_long$Signature))
df_long$Signature <- factor(df_long$Signature, levels = sig_ord)

p_all <- ggplot(df_long, aes(pt, score)) +
  geom_point(color = "black", alpha = .20, size = .15) +
  geom_smooth(color = "#e13d2d", method = "gam", formula = y ~ s(x, k = 5), se = FALSE, linewidth = 0.6) +
  facet_grid(Signature ~ group, scales = "free_y") +
  labs(x = "Pseudotime (lineage-1)", y = "UCell score") +
  theme_bw(base_size = 8) +
  theme(
    strip.text.y = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave("figures/final/fig12_Lineage1_UCell_trends_grid.png",
       p_all, width = 180, height = 8 + 15 * length(sig_ord), units = "mm", dpi = 300)
ggsave("figures/final/fig12_Lineage1_UCell_trends_grid.pdf",
       p_all, width = 180, height = 8 + 15 * length(sig_ord), units = "mm")


### ---------------- Section 9: Save final object (overwrite) ------------------
saveRDS(seu, file = "results/Final_LSK.rds")
# End of script
