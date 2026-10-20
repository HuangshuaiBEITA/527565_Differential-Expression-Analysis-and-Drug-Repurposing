# 527565_Differential-Expression-Analysis-and-Drug-Repurposing
#We established a data-driven framework bridging human transcriptomics with foundational organismal biology.
# ==============================================================================
# Supplementary Script for Differential Gene Expression (DGE) Analysis
# 
# Description: 
# This script performs covariate-adjusted DGE analysis using the 'limma' package.
# Part 1 covers the 'limma-voom' pipeline for raw count datasets.
# Part 2 covers the 'limma-trend' pipeline for FPKM/TPM normalized datasets.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Packages
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("limma", "edgeR", "ggplot2", "data.table")) # Uncomment to install

library(limma)
library(edgeR)
library(ggplot2)
library(data.table)

# ==============================================================================
# PART 1: 'limma-voom' Pipeline for Raw Count Data
# Applicable Datasets: GSE166925, GSE99248, GSE104704, GSE226151
# ==============================================================================

# 1.1 Define input and output file paths (Replace with actual file paths)
file_counts <- "path/to/raw_counts.tsv"      # Raw count matrix
file_meta   <- "path/to/metadata.csv"        # Clinical metadata (Batch, Sex, Age, Disease)
file_annot  <- "path/to/annotation.tsv"      # Gene annotation file (if applicable)
output_csv  <- "path/to/DEG_results_Covariate_Adjusted.csv"

# 1.2 Read and clean data
message("Reading raw counts and metadata...")
counts_raw <- fread(file_counts, data.table = FALSE)
rownames(counts_raw) <- counts_raw[, 1]
counts <- counts_raw[, -1]

meta <- read.csv(file_meta, row.names = 1)

# Align samples between count matrix and metadata
common_samples <- intersect(colnames(counts), rownames(meta))
if(length(common_samples) == 0) stop("Sample names do not match between counts and metadata.")
counts <- counts[, common_samples]
meta <- meta[common_samples, ]

# 1.3 Format metadata variables for linear modeling
meta$Disease <- factor(make.names(meta$Disease)) 
meta$Sex     <- factor(meta$Sex)                 
meta$Age     <- as.numeric(meta$Age)             

message("Disease group distribution:")
print(table(meta$Disease))

# 1.4 Construct Design Matrix (adjusting for Sex and Age)
design <- model.matrix(~ 0 + Disease + Sex + Age, data = meta)
colnames(design) <- gsub("^Disease", "", colnames(design))

# 1.5 Data filtering and voom transformation
dge <- DGEList(counts = counts, group = meta$Disease)
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = FALSE)

# 1.6 Differential Expression Analysis (limma)
fit <- lmFit(v, design)

# Define contrast (e.g., Patient vs Control)
dis_levels <- levels(meta$Disease)
contrast_str <- paste(dis_levels[1], "-", dis_levels[2]) 
message("Evaluating contrast: ", contrast_str)

contrast.matrix <- makeContrasts(contrasts = contrast_str, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DEG results
DEG_results <- topTable(fit2, sort.by = "P", n = Inf)
DEG_results$GeneID <- rownames(DEG_results)

# 1.7 Merge with gene annotations (Optional: Skip if annotation file is not needed)
if (file.exists(file_annot)) {
  message("Merging results with gene annotations...")
  annot <- fread(file_annot, header = TRUE, quote = "", stringsAsFactors = FALSE, data.table = FALSE)
  
  avail_cols <- intersect(c("GeneID", "Symbol", "Description"), colnames(annot))
  if(length(avail_cols) > 0) {
    DEG_results <- merge(DEG_results, annot[, avail_cols], by = "GeneID", all.x = TRUE)
  }
}

# 1.8 Sort and save results
DEG_results <- DEG_results[order(DEG_results$adj.P.Val), ]

# Reorder columns for readability
col_order <- c("GeneID", "Symbol", "Description", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
col_order <- intersect(col_order, colnames(DEG_results))
DEG_results <- DEG_results[, c(col_order, setdiff(colnames(DEG_results), col_order))]

write.csv(DEG_results, file = output_csv, row.names = FALSE, quote = FALSE)
message("Part 1 Analysis Complete. Results saved to: ", output_csv)


# ==============================================================================
# PART 2: 'limma-trend' Pipeline for FPKM/TPM Data
# Applicable Dataset: GSE135055
# ==============================================================================

# 2.1 Define input and output file paths
file_fpkm   <- "path/to/GSE135055_FPKM.txt"
file_meta   <- "path/to/GSE135055_metadata.csv"
file_annot  <- "path/to/annotation.tsv"
output_csv  <- "path/to/GSE135055_DEG_results_Trend.csv"

# 2.2 Read and clean FPKM data
message("Reading FPKM matrix and handling duplicated gene names...")
fpkm_raw <- fread(file_fpkm, data.table = FALSE)

# Detect and remove duplicate gene entries (retaining the first occurrence)
dup_idx <- duplicated(fpkm_raw[, 1])
if(any(dup_idx)) {
  message("Warning: Found ", sum(dup_idx), " duplicate gene names. Retaining the first occurrences.")
  fpkm_raw <- fpkm_raw[!dup_idx, ] 
}
rownames(fpkm_raw) <- fpkm_raw[, 1]
fpkm <- fpkm_raw[, -1] 

# Read and align metadata
meta_raw <- fread(file_meta, data.table = FALSE, encoding = "UTF-8") 
rownames(meta_raw) <- as.character(meta_raw[, 1])
meta <- meta_raw[, -1] 

common_samples <- intersect(colnames(fpkm), rownames(meta))
if(length(common_samples) == 0) stop("Sample names do not match between FPKM matrix and metadata.")

fpkm <- fpkm[, common_samples]
meta <- meta[common_samples, ]

# Format metadata variables
meta$Disease <- factor(make.names(meta$Disease)) 
meta$Sex     <- factor(meta$Sex)                 
meta$Age     <- as.numeric(meta$Age)             

# 2.3 Construct Design Matrix
design <- model.matrix(~ 0 + Disease + Sex + Age, data = meta)
colnames(design) <- gsub("^Disease", "", colnames(design))

# 2.4 Data Filtering and Log2 Transformation
message("Applying loose filtering criteria...")
# Strategy: Retain genes with FPKM > 0.1 in at least 2 samples
keep <- rowSums(fpkm > 0.1) >= 2 
fpkm_filtered <- fpkm[keep, ]
message("Genes before filtering: ", nrow(fpkm), " | Genes after filtering: ", nrow(fpkm_filtered))

# Log2 transformation: log2(FPKM + 1)
log_fpkm <- log2(fpkm_filtered + 1)

# 2.5 Differential Expression Analysis (limma-trend)
message("Running limma-trend analysis...")
fit <- lmFit(log_fpkm, design)

dis_levels <- levels(meta$Disease)
contrast_str <- paste(dis_levels[1], "-", dis_levels[2]) 
contrast.matrix <- makeContrasts(contrasts = contrast_str, levels = design) 

fit2 <- contrasts.fit(fit, contrast.matrix)
# Applying 'trend = TRUE' for log-transformed normalized data
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE) 

DEG_results <- topTable(fit2, sort.by = "P", n = Inf)
DEG_results$GeneID <- rownames(DEG_results)

# 2.6 Merge with gene annotations and save
if (file.exists(file_annot)) {
  annot <- fread(file_annot, header = TRUE, quote = "", stringsAsFactors = FALSE, data.table = FALSE)
  avail_cols <- intersect(c("GeneID", "Symbol", "Description"), colnames(annot))
  if(length(avail_cols) > 0) {
    DEG_results <- merge(DEG_results, annot[, avail_cols], by = "GeneID", all.x = TRUE)
  }
}

DEG_results <- DEG_results[order(DEG_results$P.Value), ]

col_order <- c("GeneID", "Symbol", "Description", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
col_order <- intersect(col_order, colnames(DEG_results)) 
DEG_results <- DEG_results[, c(col_order, setdiff(colnames(DEG_results), col_order))]

write.csv(DEG_results, file = output_csv, row.names = FALSE, quote = FALSE)
message("Part 2 Analysis Complete. Results saved to: ", output_csv)
