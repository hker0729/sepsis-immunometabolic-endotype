# Mitoxyperilysis 17-gene MPS program in sepsis

Analysis scripts and processed scores for the manuscript:

**"A 17-Gene Mitoxyperilysis Program Score Identifies an Immunometabolic
Endotype with Independent Prognostic Value in Sepsis"**
Kun He et al. (under review, Theranostics)

The MPS score integrates immune-sensing (TLR2/TLR4/TLR7/TLR8/MYD88/NLRP3/
TNFRSF1A/TNFRSF1B/NFE2L2), mTORC2-regulation (RICTOR/MTOR/RPTOR), and
death-execution (BAX/BAK1/BID/NINJ1) genes into a single ssGSEA score
quantifying the mitoxyperilysis program (mitochondrial dysfunction +
regulated cell death) in sepsis.

## Data sources (all public)

| Cohort | Accession | Use |
|---|---|---|
| GSE65682 | GEO, n=479 | Discovery: MPS score vs 28-day mortality (HR 0.772/SD, p=0.003) |
| GSE95233 | GEO, n=102 analyzed | External validation (OR 0.449, p=0.001, AUC 0.709) |
| E-MTAB-4421 (Davenport) | ArrayExpress, n=265 | SRS endotype x MPS interaction (OR 1.380, p=0.068) |
| GSE13904 | GEO, n=227 | Sepsis vs control discrimination (p=9.5e-7) |
| GSE167363 | GEO | Single-cell validation (21k PBMC, monocyte NINJ1 enrichment) |
| eQTLGen / IEU OpenGWAS (ieu-b-4980) | public | MR instruments and sepsis GWAS (NINJ1/BAK1/BAX) |

Expression matrices are available from NCBI GEO / ArrayExpress under the
accessions above. Processed score files are provided in `data/` (raw data
are public and reproducible from the scripts).

## Scripts (in suggested run order)

| Script | Purpose |
|---|---|
| `mps_pre_main.R` | Core 17-gene ssGSEA scoring in GSE65682 (n=479) |
| `mps_pre_modules.R` | Module decomposition: immune-sensing vs death-execution |
| `mps_pre_davenport.R` | Davenport SRS scoring in E-MTAB-4421 (z-score mean, NA-safe) |
| `mps_b_verify.R` | SRS x MPS interaction robustness (GSE65682) |
| `mps_c_95233.R` | GSE95233 external survival validation |
| `mps_d_13904.R` | GSE13904 sepsis vs control MPS comparison |
| `mps_mr_main.R` | MR: NINJ1/BAK1/BAX top cis-eQTL x sepsis GWAS (single-SNP Wald ratios; requires IEU OpenGWAS JWT) |
| `mps_mr_run.R` | MR run wrapper (JWT read from local file; token not redistributed) |
| `mps_sc_analysis.R` | Single-cell validation in GSE167363 (monocyte NINJ1, module scores) |
| `mps_strengthen.R` | Immune-axis adjustment (HR 0.816), SRS interaction, third cohort |
| `mps_supp_tables.R` | Supplementary Tables S1-S5 |
| `mps_verify_numbers.R` | Anti-fabrication verification of all manuscript numbers |
| `mps_verify2.R` | Remaining manuscript numbers (builds `data/mps_GSE65682_full.rds`) |
| `mps_cindex.R` | C-index baseline (0.731) and increment (+0.007) |
| `mps_fig1_fix.R` | Figure 1 (verified E-panel numbers: 0.772/0.449/1.380) |
| `mps_fig4_v2.R` | Figure 4: NINJ1 monocyte panels (group / T0-T6 / septic shock) |
| `mps_figfix.R` | Figures 1E/2/3 with verified numbers (C-index 0.731 -> 0.738) |

## Data files

| File | Description |
|---|---|
| `mps_GSE65682_full.rds` | GSE65682 phenotype + scores (main analysis frame) |
| `mps_GSE65682_scores.rds` | ssGSEA scores (17-gene + modules) |
| `mps_GSE95233_scores.rds` | GSE95233 scores (n=102 analyzed) |
| `mps_GSE13904_scores.rds` | GSE13904 scores |
| `mps_Davenport_scores.rds` | Davenport E-MTAB-4421 scores |
| `mps_Davenport_srs.rds` | Davenport SRS endotype labels |
| `mps_mr_wald.csv` | MR Wald-ratio results (Table S3 source) |
| `mps17_probe_map_GPL6947.csv` | Probe mapping for GPL6947 (GSE65682) |
| `mps_fig2_data.rds` | Figure 2 source data |
| `mps_eqtlgen_subset.csv` | eQTLGen cis-eQTL subset used as MR instruments |

## Key results

- GSE65682: MPS HR 0.772/SD (95% CI 0.651-0.915), p=0.003; adjusted for
  age/sex/immune score HR 0.816, p=0.026
- GSE95233: OR 0.449 (0.268-0.708) per SD, p=0.001, AUC 0.709
- Davenport: SRS2 stratum OR 1.570 (0.964-2.750), interaction p=0.308
- MR: NINJ1 OR 1.019 (F=219.9), BAK1 OR 1.020 (F=791.2), BAX OR 1.086
  (F=18.3); none significant for sepsis susceptibility

## Requirements

- R >= 4.1 (survival, pROC, GSVA/ssGSEA, Seurat for single-cell)
- IEU OpenGWAS JWT for `mps_mr_main.R` (free registration at
  https://api.opengwas.io)
- Local paths inside scripts point to the author's workspace; adjust the
  `OUT`/data paths to your own layout before running

## Contact

Kun He, Department of Critical Care Medicine, The First Hospital of
Putian, Putian, Fujian, China. hker0729@outlook.com
