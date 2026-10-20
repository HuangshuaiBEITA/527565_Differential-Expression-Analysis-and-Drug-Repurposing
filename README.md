# 527565_Differential-Expression-Analysis-and-Drug-Repurposing
#We established a data-driven framework bridging human transcriptomics with foundational organismal biology and discovered potential anti-aging drugs. The process is as follows：
Data-Driven Identification and Experimental Validation of Candesartan as a Longevity-Promoting Agent through AMPK–TOR–DAF-16 Signaling in Caenorhabditis elegans
==================================================================================================================

README: Computational Reproducibility and Analytical Workflow
================================================================================

--------------------------------------------------------------------------------
1. OVERVIEW AND GENERAL DESCRIPTION
--------------------------------------------------------------------------------
R scripts provides the complete computational 
code, and  intermediate datasets; Supplementary File S1 provides Detailed sample information;
Supplementary File S2 provides differential expression analyses ;
Supplementary File S3 provides analytical workflows required to reproduce the 
transcriptomic preprocessing, Aging-Promoting 
Immune Signature (APIS) construction, and Connectivity Map (CMap) query results 
presented in the manuscript.

--------------------------------------------------------------------------------
2. SYSTEM REQUIREMENTS AND COMPUTATIONAL ENVIRONMENT
--------------------------------------------------------------------------------
The computational pipeline was developed and executed under the following 
operating system and software specifications:

- Operating System: Windows 11 
- Base Language: R version version 4.5.1
- Core Bioconductor & CRAN Package Versions:
    * limma (version 3.64.3)
    * edgeR (version 4.6.3)
    * ggplot (version2_4.0.3)
    * data.table (version 1.18.2.1)
    * cmapR(version 1.20.0)
--------------------------------------------------------------------------------
3. DIRECTORY AND ARCHIVE STRUCTURE
--------------------------------------------------------------------------------
The unpacked archive contains the following directory tree:

Supplementary_Data_File_S1/
│
├── README.txt                     <- This file
│
├── 00_environment/
│   └── session_info.txt           <- Full R sessionInfo() details
│
├── 01_scripts/                    <- Custom R scripts (numbered by execution order)
│   ├── 01_raw_counts_limma_voom.R <- Pipeline for raw count GEO datasets
│   ├── 02_fpkm_limma_trend.R      <- Pipeline for FPKM GEO dataset (GSE135055)
│   ├── 03_APIS_construction.R     <- Integration and construction of the APIS signature
│   └── 04_cmap_query_prep.R       <- Preparation of UP/DOWN gene tags for CMap
│   └── 05_gsea_pathway_analysis. GSEA software (version 4.3.2, Broad Institute)  <- Functional enrichment  
│
├── 02_intermediate_data/          <- Processed inputs and intermediate expression matrices
│   ├── mapped_gene_matrices/      <- Probe-to-Symbol mapped expression matrices
│   └── APIS_score_matrix.rds      <- Final consolidated disease-specific APIS weights
│
└── 03_analytical_outputs/         <- Complete, unedited raw computational results
    ├── deg_results/               <- Per-dataset limma differential expression tables
    ├── cmap_results/              <- Unedited CMap ranking outputs and perturbagen scores
    └── gsea_results/              <- Full GSEA output tables and enrichment metrics

--------------------------------------------------------------------------------
4. STEP-BY-STEP EXECUTION WORKFLOW
--------------------------------------------------------------------------------
To reproduce the findings, users should set the working directory to the root 
of this folder in R (`setwd("path/to/Supplementary_Data_File_S1")`) and execute 
the scripts sequentially:

Step 1: Differential Expression Analysis 
- Script: `01_scripts/01_raw_counts_limma_voom.R`
- Description: Reads filtered raw counts for datasets GSE166925, GSE99248, 
  GSE104704, and GSE226151. Applies `edgeR::filterByExpr` and `edgeR::calcNormFactors` 
  (TMM). Employs `limma::voom` to model heteroscedasticity and fits linear models 
  incorporating covariates via design matrix `~ 0 + Disease + Sex + Age`.
- Output: Differential expression summary tables in `03_analytical_outputs/deg_results/`.

Step 2: Differential Expression Analysis for Normalized Data
- Script: `01_scripts/02_fpkm_limma_trend.R`
- Description: Applied exclusively to dataset GSE135055 (FPKM). Applies low-expression 
  filter (mean FPKM > 0.1), followed by log2-transformation. Linear modeling is 
  performed using `limma-trend` (`eBayes(..., trend=TRUE, robust=TRUE)`).
- Output: Differential expression summary tables in `03_analytical_outputs/deg_results/`.

Step 3: Signature Construction (APIS)
- Script: `01_scripts/03_APIS_construction.R`
- Description: Harmonizes directionality across datasets, integrates meta-effects, 
  and performs the Leave-One-Dataset-Out (LODO) sensitivity audit.
- Output: `02_intermediate_data/APIS_score_matrix.rds`.

Step 4: CMap Query Preparation and Target Identification
- Description: Extracts the top reversed gene sets (UP/DOWN regulated) based on the 
  APIS signature, format-compatible with the LINCS-L1000 CMap Touchstone database.

Step 5: Pathway Enrichment (GSEA)
- Description: Performs fast Gene Set Enrichment Analysis against Hallmark, KEGG, 
- Output: Enrichment dot plots and tabular summaries in `/gsea_results/`.

--------------------------------------------------------------------------------
5. THIRD-PARTY DATA AND CODE AVAILABILITY
--------------------------------------------------------------------------------
- Primary GEO Datasets: Publicly accessible from NCBI GEO under accessions:
  GSE166925, GSE99248, GSE104704, GSE226151, and GSE135055.
- CMap Touchstone Reference Database: Query scores and perturbagen mechanisms 
  depend on the CLUE platform (https://clue.io/).
- MSigDB Signatures: Downloaded from the Molecular Signatures Database 
  (http://www.gsea-msigdb.org/).

--------------------------------------------------------------------------------
6. LICENSE AND USAGE
--------------------------------------------------------------------------------
The custom analysis code included in this package is distributed under the 
MIT License and is provided "as is", without warranty of any kind. 
For questions or bug reports, please contact the corresponding author(s).
================================================================================
