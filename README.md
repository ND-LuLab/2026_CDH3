# CDH3 analysis code

Analysis and visualization code for studying CDH3 expression, molecular subtypes,
and signaling in prostate and breast cancer. Scripts are organized by cohort or
analysis topic. Analyses plotted in GraphPad Prism end with table exports.

## Code and figure panels

| Script | Analysis | Figure panels | Output |
|---|---|---|---|
| [01_tang_crpc_rnaseq.R](analysis/01_tang_crpc_rnaseq.R) | Tang CRPC expression and pathway signatures | 1B–C, S1B | Heatmap and GraphPad tables |
| [02_tang_regulatory_tables.R](analysis/02_tang_regulatory_tables.R) | CDH3–TF correlations and YAP/TAZ knockdown | 3A, 3I | Correlation and normalized-expression tables |
| [03_prostate_pcs_pam50.R](analysis/03_prostate_pcs_pam50.R) | Prostate cohort PCS/PAM50 classification | 1D–E, S1C | Sankeys and GraphPad tables |
| [04_brca_subtype_sankey.R](analysis/04_brca_subtype_sankey.R) | TCGA BRCA biomarker/PAM50 classification | 1F, S1F | Sankey and GraphPad tables |
| [05_song_prostate_single_cell.R](analysis/05_song_prostate_single_cell.R) | Song prostate single-cell analysis | 2A–D, S2A | Annotation, feature plots, dot plots and tables |
| [06_wong_prostate_single_cell.R](analysis/06_wong_prostate_single_cell.R) | Wong prostate single-cell analysis | S2B–E | Annotation, feature plots, dot plots and tables |
| [07_jee_prostate_single_cell.R](analysis/07_jee_prostate_single_cell.R) | Jee prostate single-cell analysis | S2F–H | Annotation, feature plot, dot plot and tables |
| [08_jin_breast_single_cell.R](analysis/08_jin_breast_single_cell.R) | Jin TNBC single-cell analysis | 2E–F | Annotation, feature plots and CDH3 dot plot |
| [09_xu_breast_single_cell.R](analysis/09_xu_breast_single_cell.R) | Xu breast single-cell analysis | S2I–J | Annotation, feature plots and CDH3 dot plot |
| [10_su2c_arsi_survival.R](analysis/10_su2c_arsi_survival.R) | SU2C time on first-line ARSI | S1G | Kaplan–Meier plots, risk tables and statistics |
| [11_primary_metastatic_comparisons.R](analysis/11_primary_metastatic_comparisons.R) | FHCRC and AURORA primary/metastatic comparisons | S1H–I | Paired/unpaired comparisons and plots |
| [12_prostate_public_subtype_tables.R](analysis/12_prostate_public_subtype_tables.R) | TCGA PRAD and metastatic FHCRC PAM50 classification | S1D–E | GraphPad tables |


Figure 3H was visualized in IGV using control, YAP and FOSL1 bigWig tracks for
DU145 and MSKPCa3 (hg38, chr16:68577374–68801708; track range 0–1).
No R visualization script is required for that panel.

## Usage

Edit the **Local settings** block at the top of the desired script to set input
paths and the output directory. Each script includes its supporting functions
and can be run independently. Relative paths resolve from the working directory.

Run from the repository root in a fresh R session:

```r
source("analysis/05_song_prostate_single_cell.R")
```

Or use a terminal:

```sh
Rscript analysis/03_prostate_pcs_pam50.R
Rscript analysis/10_su2c_arsi_survival.R
```

The regulatory analysis has two stages:

```sh
Rscript analysis/02_tang_regulatory_tables.R correlation
Rscript analysis/02_tang_regulatory_tables.R knockdown
```

Outputs are written to topic-specific folders under `results/`. Each analysis
records its R environment in `sessionInfo.txt`. Datasets and packages are not
downloaded automatically.

## Input data

External data and processed R objects are not included. Required input paths
are listed at the top of each script.

| Analysis | Inputs |
|---|---|
| Tang CRPC | Transformed gene-by-model matrix; model metadata with `group`; gene map with `ensembl_gene_id_version`, `external_gene_name`, `chromosome_name` |
| Regulatory correlations | Transformed expression, gene map and transcription-factor list with a `tf` column |
| Knockdown | Gene-by-sample counts; sample sheet with `sample` and `group`; gene map. Groups: scr, TAZ_KD, YAP_KD, YAP_TAZ_KD |
| Prostate cohort subtyping | Normalized expression with a `gene` column; PAM50 scores with sample IDs and Basal/LumA/LumB columns |
| TCGA BRCA | Biomarker/PAM50 clinical calls and expression z-scores |
| Single-cell | Processed Seurat objects with RNA expression, stored UMAP coordinates, `seurat_clusters` and annotation metadata |
| SU2C survival | Expression matrix, sample-to-patient clinical table and sheet 2 of `Abida2019_TableS1.xlsx` |
| Primary/metastatic comparisons | AURORA normalized expression; FHCRC standardized expression and clinical sample data |
| Public prostate subtypes | TCGA PRAD/FHCRC expression z-scores; FHCRC clinical sample annotations |

The single-cell cohorts are Song (GSE176031), Wong (GSE185344),
Jee (GSE221603), Jin (GSE263995) and Xu (GSE180286). Input filenames are
`song.rds`, `wong.rds`, `lee.rds` (Jee), `Jin.rds` and `Xu.rds`.

Single-cell scripts begin with processed objects and apply dataset-specific
cell annotations and epithelial-cell selection. Cluster IDs and UMAP coordinate
gates depend on those objects and must be reassessed if clustering is repeated.
Raw-read preprocessing and integration are not performed by these scripts.

The Tang workflow starts from transformed expression by default. Optional
DESeq2 processing requires verified raw counts; transformed expression must not
be passed to DESeq2 as counts. The prostate cohort script includes an optional
function for CEL RMA normalization and probe-to-gene aggregation. Sample
exclusions are specified explicitly.

## Dependencies

Install the packages needed for the relevant topic in a compatible R and
Bioconductor environment.

| Analysis | R packages |
|---|---|
| Tang CRPC | pheatmap, GSVA; optional DESeq2, SummarizedExperiment |
| Knockdown | edgeR, limma |
| Prostate cohort subtyping | singscore, ggplot2, ggalluvial, dplyr, svglite; optional oligo, Biobase and the platform design package |
| TCGA BRCA | ggplot2, ggalluvial, dplyr, svglite |
| Single-cell | Seurat, ggplot2, RColorBrewer, patchwork; prostate WNT plots also dplyr, ggh4x |
| SU2C survival | data.table, dplyr, ggplot2, survival, patchwork, readxl |
| Primary/metastatic comparisons | data.table, dplyr, ggplot2, ggsignif, tidyplots |
| Public prostate subtypes | genefu |

The survival plotting code uses ggplot2's `legend.position.inside` and
patchwork's `free()`. Survival and primary/metastatic PDF exports require Cairo
support. Seurat objects must be compatible with the installed Seurat version.

## GraphPad Prism

Import `*_GraphPad.csv` files as Column tables. Each column represents a group;
blank padding cells are missing values, not zeroes. Companion long-format tables
retain sample identifiers. Use the corresponding figure legend for statistical
comparisons, multiple-comparison adjustment, whiskers and error bars.

The regulatory knockdown workflow exports TMM-normalized CPM and logCPM.
Relative-to-control normalization and final Prism formatting are separate steps.
The correlation workflow requires an explicitly supplied transcription-factor
list. These outputs are analysis tables, not saved Prism projects.

## Analysis notes

- The Tang heatmap displays 40 CRPC models in AR/WNT/NE/SCL order and uses
  gene-wise standardized expression. GSVA uses a Gaussian kernel,
  `maxDiff = FALSE` and gene-set sizes of 5–500.
- Prostate PCS-associated scores use 37-gene singscore sets, cohort-wise score
  standardization and maximum-score assignments. PAM50 assignments in prostate
  cohorts are restricted to Basal, LumA and LumB. CDH3 expression is standardized
  across each cohort, not separately within each subtype.
- The survival endpoint is time on first-line ARSI. The analysis retains the
  first merged sample per patient, excludes median ties, and includes boundary
  ties in upper/lower tertile groups.
- Paired primary/metastatic tests use patient-level means; plotted points remain
  sample-level observations. Check plotting-library whisker defaults when
  matching a specified box-plot convention.
- Scripts accept external inputs and do not include a complete execution
  environment. Syntax has been checked with R 4.5.3; end-to-end execution and
  exact figure equivalence have not been validated with the external datasets.
