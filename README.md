# Neutrophil–HSPC-feedback

**Manuscript title**  
**A neutrophil–progenitor positive feedback loop via ROS–FOXO1–DRP1 signaling sustains splenic myelopoiesis and promotes immune evasion in cancer**

**Abstract**  
Splenic myelopoiesis supplies immunosuppressive myeloid cells in cancer, but its local regulatory mechanism remains unclear. Here we identified a neutrophil–hematopoietic stem and progenitor cell (HSPC) positive feedback loop that sustained splenic myelopoiesis. Tumor educated neutrophils became the dominant source of NOX2 dependent reactive oxygen species (ROS) in the spleen. Neutrophil ROS drove FOXO1 activation and DRP1 mediated mitochondrial fission in neighboring HSPCs, reprogramming them toward myeloid biased expansion and suppressive output. These progeny neutrophils, in turn, reinforced ROS stress and perpetuated the loop. Disrupting this circuit by limiting neutrophil ROS or inhibiting HSPC mitochondrial fission dampened splenic myelopoiesis, curtailed suppressive myeloid cell generation, and enhanced cytotoxic T cell activity, thereby restoring antitumor immunity and restraining tumor progression. Analyses of human spleens and cord-blood HSPCs support conservation of this neutrophil–ROS–mitochondrial axis. These findings reveal a self-propagating neutrophil–HSPC feedback circuit that sustains tumor-promoting splenic myelopoiesis and represents a tractable target to normalize tumor immunity.

---

## Repository overview

This repository hosts **lightly annotated, de-personalized scripts** to reproduce the single-cell RNA-seq analysis of the neutrophil–HSPC feedback loop described in the study.  
All absolute paths were anonymized to placeholders (e.g., `path_to_data`, `path_to_.radian_profile`), and date suffixes in output filenames were removed for portability.

---

### Directory layout

```
scRNA-seq/
├── Analysis/
│   ├── 01_hepaLSK_Data_Integration.R      # Integration (Seurat SCT-RPCA) + UMAP
│   ├── 01.5_hepaLSK_pySCENIC.R            # Export counts to LOOM; import pySCENIC AUC/Bin back to Seurat
│   ├── 02_hepaLSK_Clustering.r            # Clustering + hierarchical annotation + marker tables
│   ├── 03_hepaLSK_Group_Comparisons.r     # Group-wise DE (RNA/AUC/Bin) and intersections
│   ├── 04_hepaLSK_Trajectory.r            # Slingshot trajectory inference
│   └── 05_hepaLSK_final.r                 # Final figure panels and UCell pseudotime analysis
└── QC/
    ├── run_cellbender_batch.sh            # Batch CellBender remove-background
    ├── scRNA_QC_loop_CB_step1.R           # Pre-QC: normalization, HVG, DR/cluster, doublet detection
    └── scRNA_QC_loop_CB_step2.R           # Final filtering and export per-sample Seurat objects
```

---

### Pipeline overview

```
CellRanger → CellBender → QC step1 → QC step2 → Integration → (optional) pySCENIC ↔ R
      │          │             │           │           │
      └─ run_cellbender_batch.sh
                 │
      └─ scRNA_QC_loop_CB_step1.R
                 │
      └─ scRNA_QC_loop_CB_step2.R        (per-sample Seurat RDS)
                 │
      └─ 01_hepaLSK_Data_Integration.R   (integrated Seurat)
                 │
      └─ 01.5_hepaLSK_pySCENIC.R         (export LOOM → run pySCENIC → import AUC/Bin)
                 │
      └─ 02_hepaLSK_Clustering.r         (annotations & marker tables)
                 │
      └─ 03_hepaLSK_Group_Comparisons.r  (RNA/AUC/Bin DE; intersections)
                 │
      └─ 04_hepaLSK_Trajectory.r         (Slingshot trajectory)
                 │
      └─ 05_hepaLSK_final.r              (final figures & UCell trends)
```

---

## Quick start

1. **Clone this repo**
   ```bash
   git clone https://github.com/ChongEdwardWu/Neutrophil-HSPC-feedback.git
   cd Neutrophil-HSPC-feedback/scRNA-seq
   ```

2. **Prepare inputs**
   - Each sample: 10x `raw_feature_bc_matrix/` output.
   - Run **CellBender** to produce `<sample>_filtered.h5` under `03_cellbender/<sample>/`.

3. **Edit placeholders**
   Replace all placeholder paths in scripts (`path_to_data`, `path_to_.radian_profile`, etc.) with your actual directories.

4. **Run in order**
   ```text
   QC/run_cellbender_batch.sh
   QC/scRNA_QC_loop_CB_step1.R
   QC/scRNA_QC_loop_CB_step2.R
   Analysis/01_hepaLSK_Data_Integration.R
   Analysis/01.5_hepaLSK_pySCENIC.R
   Analysis/02_hepaLSK_Clustering.r
   Analysis/03_hepaLSK_Group_Comparisons.r
   Analysis/04_hepaLSK_Trajectory.r
   Analysis/05_hepaLSK_final.r
   ```

---

## Software dependencies

### CellBender
- Run through `QC/run_cellbender_batch.sh`
- Requires a Conda environment with `cellbender` installed and accessible via `conda activate cellbender`.

### R packages
- **Core:** `Seurat`, `SeuratWrappers`, `Matrix`, `patchwork`, `ggplot2`, `tidyverse`, `openxlsx2`
- **Bioconductor:** `SingleCellExperiment`, `scater`, `scran`, `scuttle`, `DropletUtils`, `BiocParallel`, `BiocSingular`, `scDblFinder`, `bluster`
- **Optional:** `SingleR`, `scMCA`, `Azimuth`, `UCell`, `msigdbr`, `pheatmap`, `ggVennDiagram`, `ggrepel`, `viridis`, `clustree`, `monocle3`, `slingshot`

### Python (for pySCENIC)
- `pyscenic`, `loompy`, `pandas`
- Requires ranking DBs, motif annotation (`.tbl`), and TF list (`.txt`) under a species-specific folder, configured in `run_SCENIC_mouse.py`.


---

## Notes

- **All absolute paths** have been replaced with placeholders to protect privacy.  
- **No timestamps** are appended to output files to improve reproducibility.  
- Each step uses `set.seed(123)` and `--seed 123` for deterministic results.  
- Optional tools (`Azimuth`, `SingleR`, `AUCell`, etc.) are included as commented blocks for user flexibility.

---

## Outputs summary

| Script | Main outputs |
|--------|---------------|
| `scRNA_QC_loop_CB_step1.R` | `01_QC_step1_<sample>.rds` + QC figures |
| `scRNA_QC_loop_CB_step2.R` | `02_QC_step2_<sample>_seu.rds` |
| `01_hepaLSK_Data_Integration.R` | `01_Integration_LSK.rds` |
| `01.5_hepaLSK_pySCENIC.R` | `01.5_pySCENIC2seurat.rds` + pySCENIC result tables |
| `02_hepaLSK_Clustering.r` | `02_Annotation_LSK.rds`, `Step2_cluster_celltype.xlsx` |
| `03_hepaLSK_Group_Comparisons.r` | `Step3_Group_Comparisons.xlsx` |
| `04_hepaLSK_Trajectory.r` | `Step4_LSKsce_for_slingshot.rds` |
| `05_hepaLSK_final.r` | Final PNG/PDF figures + `Final_LSK.rds` |

---

## Citation

If you use these scripts, please cite the corresponding manuscript above.  
Please also cite relevant software tools (Seurat, scran, pySCENIC, Slingshot, UCell, etc.) as appropriate.

---

## License

This repository is provided for academic and non-commercial use.  
For redistribution or commercial use, please open an issue to discuss licensing.

---

## Contact

For any questions or reproducibility requests, please email to wuchong5@mail(dot)sysu(dot)edu(dot)cn.

---

**Acknowledgements**  
We thank developers of Seurat, Bioconductor’s single-cell packages, pySCENIC, Slingshot, and UCell for their foundational tools enabling this analysis.
