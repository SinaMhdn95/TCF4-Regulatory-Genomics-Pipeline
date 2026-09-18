####################################################################
# Doosparats Deconvolution of Transcriptional Networks Identifies
# TCF4 as a Master Regulator in Schizophrenia
# # Using KD_raw2 only
####################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Load libraries
library(GenomicRanges)  # handling genomic intervals and peak overlap analysis
library(rtracklayer)    # importing and exporting BED, GFF, GTF, BigWig, and other genomic files
library(ChIPseeker)     # peak annotation
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # Human hg38 transcript annotation database used by ChIPseeker
library(org.Hs.eg.db)  # Human gene annotation database, used to map Entrez IDs to gene symbols
library(clusterProfiler)  # For GO, KEGG, and pathway enrichment analysis
library(txdbmaker)     # chromosome naming style, such as converting 1 to chr1
library(GenomeInfoDb)  # enrichment plots, including dotplot and barplot
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationDbi)
library(GenomeInfoDb)
library(tibble)
library(pheatmap)
library(DESeq2)
library(limma)
library(ggplot2)


txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Import files
npc_peaks     <- import("NPC_hg38_summit_500bp.bed")
mcclay_peaks  <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
forrest_peaks <- import("Forrest_hg38.bed")

# Import RNA-seq TCF4 KD data
KD_raw2 <- read.table(
  "GSE128333_Data2.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

head(KD_raw2)


# Select samples
# Healthy Control Donor CD07
# Day 14 Glutamatergic Neurons
# 1. Select CD07 samples
cd07_counts <- KD_raw2[, c(
  "D14_7_C2",
  "D14_7_C3",
  "D14_7_T1",
  "D14_7_T2",
  "D14_7_T3"
)]

cd07_counts <- round(
  as.matrix(cd07_counts)
)

# 2. Create sample metadata
cd07_coldata <- data.frame(

  condition = c(
    "Ctrl",
    "Ctrl",
    "KD",
    "KD",
    "KD"
  ),

  batch = c(
    "Batch1",  # D14_7_C2
    "Batch2",  # D14_7_C3
    "Batch1",  # D14_7_T1
    "Batch2",  # D14_7_T2
    "Batch1"   # D14_7_T3
  )
)

rownames(cd07_coldata) <- colnames(cd07_counts)

cd07_coldata$condition <- factor(
  cd07_coldata$condition,
  levels = c("Ctrl", "KD")
)

cd07_coldata$batch <- factor(
  cd07_coldata$batch
)

all(rownames(cd07_coldata) == colnames(cd07_counts))

cd07_coldata

# 3. DESeq2 without batch correction
dds_cd07_no_batch <- DESeqDataSetFromMatrix(
  countData = cd07_counts,
  colData = cd07_coldata,
  design = ~ condition
)

dds_cd07_no_batch <- dds_cd07_no_batch[
  rowSums(counts(dds_cd07_no_batch)) >= 10,
]

dds_cd07_no_batch <- DESeq(
  dds_cd07_no_batch
)

cd07_res_no_batch <- results(
  dds_cd07_no_batch,
  contrast = c("condition", "KD", "Ctrl")
)

cd07_res_no_batch_df <- as.data.frame(
  cd07_res_no_batch
)

cd07_res_no_batch_df$ENSEMBL <- rownames(
  cd07_res_no_batch_df
)

cd07_res_no_batch_df$ENSEMBL_clean <- sub(
  "\\..*",
  "",
  cd07_res_no_batch_df$ENSEMBL
)

cd07_res_no_batch_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = cd07_res_no_batch_df$ENSEMBL_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

cd07_res_no_batch_df_clean <- cd07_res_no_batch_df[
  !is.na(cd07_res_no_batch_df$SYMBOL),
]

cd07_sig_no_batch <- subset(
  cd07_res_no_batch_df_clean,
  padj < 0.05
)


# 4. DESeq2 with batch correction
dds_cd07_batch <- DESeqDataSetFromMatrix(
  countData = cd07_counts,
  colData = cd07_coldata,
  design = ~ batch + condition
)

dds_cd07_batch <- dds_cd07_batch[
  rowSums(counts(dds_cd07_batch)) >= 10,
]

dds_cd07_batch <- DESeq(
  dds_cd07_batch
)

cd07_res_batch <- results(
  dds_cd07_batch,
  contrast = c("condition", "KD", "Ctrl")
)

cd07_res_batch_df <- as.data.frame(
  cd07_res_batch
)

cd07_res_batch_df$ENSEMBL <- rownames(
  cd07_res_batch_df
)

cd07_res_batch_df <- cd07_res_batch_df[
  order(cd07_res_batch_df$padj),
]

cd07_res_batch_df$ENSEMBL_clean <- sub(
  "\\..*",
  "",
  cd07_res_batch_df$ENSEMBL
)

cd07_res_batch_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = cd07_res_batch_df$ENSEMBL_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

cd07_res_batch_df_clean <- cd07_res_batch_df[
  !is.na(cd07_res_batch_df$SYMBOL),
]

View(cd07_res_batch_df_clean)

# 5. Create DEG lists after batch correction
cd07_sig <- subset(
  cd07_res_batch_df_clean,
  padj < 0.05
)

cd07_up <- subset(
  cd07_res_batch_df_clean,
  padj < 0.05 & log2FoldChange > 0
)

cd07_down <- subset(
  cd07_res_batch_df_clean,
  padj < 0.05 & log2FoldChange < 0
)

cd07_background_genes <- unique(
  cd07_res_batch_df_clean$SYMBOL
)

cd07_sig_genes <- unique(
  cd07_sig$SYMBOL
)

cd07_up_genes <- unique(
  cd07_up$SYMBOL
)

cd07_down_genes <- unique(
  cd07_down$SYMBOL
)


# 6. Compare DEG counts before and after batch correction
cd07_batch_summary <- data.frame(

  Analysis = c(
    "No_batch",
    "Batch_corrected"
  ),

  Significant_DEGs = c(
    nrow(cd07_sig_no_batch),
    length(cd07_sig_genes)
  ),

  Upregulated = c(
    sum(cd07_res_no_batch_df_clean$padj < 0.05 &
          cd07_res_no_batch_df_clean$log2FoldChange > 0,
        na.rm = TRUE),
    length(cd07_up_genes)
  ),

  Downregulated = c(
    sum(cd07_res_no_batch_df_clean$padj < 0.05 &
          cd07_res_no_batch_df_clean$log2FoldChange < 0,
        na.rm = TRUE),
    length(cd07_down_genes)
  )
)

cd07_batch_summary


# 7. PCA before batch correction
vsd_cd07 <- vst(
  dds_cd07_batch,
  blind = FALSE
)

cd07_vst_mat <- assay(
  vsd_cd07
)

cd07_pca_before <- prcomp(
  t(cd07_vst_mat),
  scale. = FALSE
)

percent_var_before <- round(
  100 * cd07_pca_before$sdev^2 /
    sum(cd07_pca_before$sdev^2),
  2
)

cd07_pca_before_df <- data.frame(
  Sample = rownames(cd07_pca_before$x),
  PC1 = cd07_pca_before$x[, 1],
  PC2 = cd07_pca_before$x[, 2],
  condition = cd07_coldata$condition,
  batch = cd07_coldata$batch
)

p_cd07_before <- ggplot(
  cd07_pca_before_df,
  aes(
    x = PC1,
    y = PC2,
    color = condition,
    shape = batch,
    label = Sample
  )
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.8,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    title = "CD07 PCA before batch correction",
    x = paste0("PC1: ", percent_var_before[1], "% variance"),
    y = paste0("PC2: ", percent_var_before[2], "% variance")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_cd07_before

# 8. Remove batch effect from VST matrix for PCA only
cd07_vst_batch_corrected <- removeBatchEffect(
  cd07_vst_mat,
  batch = cd07_coldata$batch,
  design = model.matrix(
    ~ condition,
    data = cd07_coldata
  )
)

# 9. PCA after batch correction
cd07_pca_after <- prcomp(
  t(cd07_vst_batch_corrected),
  scale. = FALSE
)

percent_var_after <- round(
  100 * cd07_pca_after$sdev^2 /
    sum(cd07_pca_after$sdev^2),
  2
)

cd07_pca_after_df <- data.frame(
  Sample = rownames(cd07_pca_after$x),
  PC1 = cd07_pca_after$x[, 1],
  PC2 = cd07_pca_after$x[, 2],
  condition = cd07_coldata$condition,
  batch = cd07_coldata$batch
)

p_cd07_after <- ggplot(
  cd07_pca_after_df,
  aes(
    x = PC1,
    y = PC2,
    color = condition,
    shape = batch,
    label = Sample
  )
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.8,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    title = "CD07 PCA after batch correction",
    x = paste0("PC1: ", percent_var_after[1], "% variance"),
    y = paste0("PC2: ", percent_var_after[2], "% variance")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_cd07_after


# 11. Final DEG summary after batch correction
data.frame(
  Category = c(
    "Background",
    "Significant",
    "Up",
    "Down"
  ),
  Gene_Count = c(
    length(cd07_background_genes),
    length(cd07_sig_genes),
    length(cd07_up_genes),
    length(cd07_down_genes)
  )
)

####################################################################
# Healthy Control Donor CD09
# Day 14 Glutamatergic Neurons
# TCF4 Knockdown vs Control
####################################################################
# Select CD09 samples
cd09_counts <- KD_raw2[, c(
  "D14_9_C1",
  "D14_9_C2",
  "D14_9_C3",
  "D14_9_T1",
  "D14_9_T2",
  "D14_9_T3"
)]

# Convert to integer matrix
cd09_counts <- round(
  as.matrix(cd09_counts)
)


# Create sample information
cd09_coldata <- data.frame(
  condition = c(
    "Ctrl",
    "Ctrl",
    "Ctrl",
    "KD",
    "KD",
    "KD"
  )
)

rownames(cd09_coldata) <- colnames(cd09_counts)

cd09_coldata$condition <- factor(
  cd09_coldata$condition,
  levels = c("Ctrl", "KD")
)

all(
  rownames(cd09_coldata) ==
    colnames(cd09_counts)
)


# DESeq2
dds_cd09 <- DESeqDataSetFromMatrix(
  countData = cd09_counts,
  colData = cd09_coldata,
  design = ~ condition
)

# Filter low-count genes
dds_cd09 <- dds_cd09[
  rowSums(counts(dds_cd09)) >= 10,
]

# Run DESeq2
dds_cd09 <- DESeq(
  dds_cd09
)


# Extract results
cd09_res <- results(
  dds_cd09,
  contrast = c(
    "condition",
    "KD",
    "Ctrl"
  )
)

cd09_res_df <- as.data.frame(
  cd09_res
)

cd09_res_df$ENSEMBL <- rownames(
  cd09_res_df
)

cd09_res_df <- cd09_res_df[
  order(cd09_res_df$padj),
]


# Annotate genes
cd09_res_df$ENSEMBL_clean <- sub(
  "\\..*",
  "",
  cd09_res_df$ENSEMBL
)

cd09_res_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = cd09_res_df$ENSEMBL_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

cd09_res_df_clean <- cd09_res_df[
  !is.na(cd09_res_df$SYMBOL),
]

View(cd09_res_df_clean)


# Create DEG lists
cd09_sig <- subset(
  cd09_res_df_clean,
  padj < 0.05
)

cd09_up <- subset(
  cd09_res_df_clean,
  padj < 0.05 &
    log2FoldChange > 0
)

cd09_down <- subset(
  cd09_res_df_clean,
  padj < 0.05 &
    log2FoldChange < 0
)


# Create gene lists
cd09_background_genes <- unique(
  cd09_res_df_clean$SYMBOL
)

cd09_sig_genes <- unique(
  cd09_sig$SYMBOL
)

cd09_up_genes <- unique(
  cd09_up$SYMBOL
)

cd09_down_genes <- unique(
  cd09_down$SYMBOL
)


# Summary
data.frame(
  Category = c(
    "Background",
    "Significant",
    "Up",
    "Down"
  ),

  Gene_Count = c(
    length(cd09_background_genes),
    length(cd09_sig_genes),
    length(cd09_up_genes),
    length(cd09_down_genes)
  )
)



# PCA
vsd_cd09 <- vst(
  dds_cd09,
  blind = FALSE
)

plotPCA(
  vsd_cd09,
  intgroup = "condition"
)


####################################################################
# Select CD09 samples
####################################################################
# CD09 is a healthy control hiPSC line.
# These samples are Day 14 glutamatergic neurons.
# C = control
# T = TCF4 knockdown

cd09_counts <- KD_raw2[, c(
  "D14_9_C1",
  "D14_9_C2",
  "D14_9_C3",
  "D14_9_T1",
  "D14_9_T2",
  "D14_9_T3"
)]

# edgeR requires raw integer counts.
# round() makes sure all values are integers.
cd09_counts <- round(
  as.matrix(cd09_counts)
)

## Check sample names
colnames(cd09_counts)

# Create CD09 sample metadata
# This tells edgeR which samples are Control and which are KD.
# The order must match the column order in cd09_counts.
cd09_sample_info <- data.frame(
  condition = c(
    "Ctrl",
    "Ctrl",
    "Ctrl",
    "KD",
    "KD",
    "KD"
  )
)

rownames(cd09_sample_info) <- colnames(cd09_counts)

# Set Ctrl as the reference group.
# This means logFC will be KD compared with Ctrl.
cd09_sample_info$condition <- factor(
  cd09_sample_info$condition,
  levels = c("Ctrl", "KD")
)

# Confirm metadata rows match count columns
all(rownames(cd09_sample_info) == colnames(cd09_counts))

cd09_sample_info


# Create edgeR DGEList object
# DGEList stores the count matrix and group information for edgeR.
cd09_dge <- DGEList(
  counts = cd09_counts,
  group = cd09_sample_info$condition
)



# Filter low-expression genes
# filterByExpr() removes genes with very low counts.
# This is the recommended edgeR filtering step.
# It keeps genes with enough expression for reliable testing.
cd09_keep_genes <- filterByExpr(
  cd09_dge,
  group = cd09_sample_info$condition
)

table(cd09_keep_genes)

cd09_dge <- cd09_dge[
  cd09_keep_genes,
  ,
  keep.lib.sizes = FALSE
]

# Normalize library sizes
# calcNormFactors() performs TMM normalization.
# This corrects for differences in sequencing depth and RNA composition.
cd09_dge <- calcNormFactors(
  cd09_dge
)

cd09_dge$samples


# Create design matrix
# The design matrix tells edgeR what comparison to test.
# Here we test the effect of condition:
# KD vs Ctrl
cd09_design <- model.matrix(
  ~ condition,
  data = cd09_sample_info
)

cd09_design

# The coefficient "conditionKD" represents:
# TCF4 knockdown compared with Control.


# Estimate dispersion
# Dispersion measures biological variability among replicates.
# edgeR estimates this before differential expression testing.
cd09_dge <- estimateDisp(
  cd09_dge,
  cd09_design
)



# Fit edgeR GLM model
# glmFit() fits a negative binomial model for each gene.
cd09_fit_lrt <- glmFit(
  cd09_dge,
  cd09_design
)


# Run edgeR likelihood ratio test
# glmLRT() tests whether KD is significantly different from Ctrl.
# coef = "conditionKD" tests the KD effect.
cd09_lrt <- glmLRT(
  cd09_fit_lrt,
  coef = "conditionKD"
)

# Extract edgeR LRT results
# topTags(..., n = Inf) returns all tested genes, not only the top genes.
cd09_edgeR_lrt_df <- topTags(
  cd09_lrt,
  n = Inf
)$table

# Add Ensembl IDs as a column.
cd09_edgeR_lrt_df$ENSEMBL <- rownames(
  cd09_edgeR_lrt_df
)

# Sort by FDR.
cd09_edgeR_lrt_df <- cd09_edgeR_lrt_df[
  order(cd09_edgeR_lrt_df$FDR),
]

head(cd09_edgeR_lrt_df)



# Annotate Ensembl IDs to gene symbols
# Remove Ensembl version numbers if present.
# Example: ENSG000001234.5 becomes ENSG000001234
cd09_edgeR_lrt_df$ENSEMBL_clean <- sub(
  "\\..*",
  "",
  cd09_edgeR_lrt_df$ENSEMBL
)

# Map Ensembl IDs to gene symbols.
cd09_edgeR_lrt_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = cd09_edgeR_lrt_df$ENSEMBL_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove genes that could not be mapped to symbols.
cd09_edgeR_lrt_df_clean <- cd09_edgeR_lrt_df[
  !is.na(cd09_edgeR_lrt_df$SYMBOL),
]

head(cd09_edgeR_lrt_df_clean)



# Create significant DEG tables
# Significant genes:
# FDR < 0.05
cd09_edgeR_lrt_sig <- subset(
  cd09_edgeR_lrt_df_clean,
  FDR < 0.05
)

# Upregulated after TCF4 knockdown:
# logFC > 0 means higher expression in KD than Ctrl.
cd09_edgeR_lrt_up <- subset(
  cd09_edgeR_lrt_df_clean,
  FDR < 0.05 & logFC > 0
)

# Downregulated after TCF4 knockdown:
# logFC < 0 means lower expression in KD than Ctrl.
cd09_edgeR_lrt_down <- subset(
  cd09_edgeR_lrt_df_clean,
  FDR < 0.05 & logFC < 0
)


# Create gene symbol lists
# Background genes = all genes tested by edgeR and successfully mapped.
cd09_edgeR_lrt_background_genes <- unique(
  cd09_edgeR_lrt_df_clean$SYMBOL
)

# Significant DEG symbols.
cd09_edgeR_lrt_sig_genes <- unique(
  cd09_edgeR_lrt_sig$SYMBOL
)

# Upregulated DEG symbols.
cd09_edgeR_lrt_up_genes <- unique(
  cd09_edgeR_lrt_up$SYMBOL
)

# Downregulated DEG symbols.
cd09_edgeR_lrt_down_genes <- unique(
  cd09_edgeR_lrt_down$SYMBOL
)


# Summary table
cd09_edgeR_lrt_summary <- data.frame(
  Category = c(
    "Background",
    "Significant",
    "Up",
    "Down"
  ),

  Gene_Count = c(
    length(cd09_edgeR_lrt_background_genes),
    length(cd09_edgeR_lrt_sig_genes),
    length(cd09_edgeR_lrt_up_genes),
    length(cd09_edgeR_lrt_down_genes)
  )
)

cd09_edgeR_lrt_summary


####################################################################
# CD09 edgeR LRT vs TCF4 Peak-Associated Gene Sets
# Goal:
# 1. Import TCF4 peak-associated gene lists
# 2. Compare each list with CD09 edgeR LRT DEGs
# 3. Run Fisher enrichment tests
####################################################################

# Import TCF4 peak-associated gene lists
# These CSV files contain genes associated with each TCF4 peak dataset.
# In your files, the gene symbols are stored in column "x".
npc_genes_hg38_df <- read.csv(
  "npc_genes_hg38.csv"
)

mcclay_genes_hg38_df <- read.csv(
  "mcclay_genes_hg38.csv"
)

forrest_genes_hg38_df <- read.csv(
  "forrest_genes_hg38.csv"
)

# Check column names
colnames(npc_genes_hg38_df)
colnames(mcclay_genes_hg38_df)
colnames(forrest_genes_hg38_df)

# Extract gene symbols from column x
peak_npc_genes <- unique(
  na.omit(npc_genes_hg38_df$x)
)

peak_mcclay_genes <- unique(
  na.omit(mcclay_genes_hg38_df$x)
)

peak_forrest_genes <- unique(
  na.omit(forrest_genes_hg38_df$x)
)

# Check gene list sizes
length(peak_npc_genes)
length(peak_mcclay_genes)
length(peak_forrest_genes)



# Prepare CD09 edgeR LRT DEG gene sets
# Background genes:
# All genes tested by edgeR LRT and successfully mapped to gene symbols.
cd09_background_genes <- unique(
  cd09_edgeR_lrt_df_clean$SYMBOL
)

# Significant genes:
# Genes with FDR < 0.05 after TCF4 knockdown.
cd09_sig_genes <- unique(
  cd09_edgeR_lrt_sig$SYMBOL
)

# Non-significant genes:
# Background genes that are not significant DEGs.
cd09_nonsig_genes <- setdiff(
  cd09_background_genes,
  cd09_sig_genes
)

# Upregulated genes:
# Significant genes with logFC > 0.
cd09_up_genes <- unique(
  cd09_edgeR_lrt_up$SYMBOL
)

# Downregulated genes:
# Significant genes with logFC < 0.
cd09_down_genes <- unique(
  cd09_edgeR_lrt_down$SYMBOL
)

# Check CD09 gene set sizes
data.frame(
  Category = c(
    "Background",
    "Significant",
    "Non_significant",
    "Upregulated",
    "Downregulated"
  ),
  Gene_Count = c(
    length(cd09_background_genes),
    length(cd09_sig_genes),
    length(cd09_nonsig_genes),
    length(cd09_up_genes),
    length(cd09_down_genes)
  )
)



# Restrict TCF4 gene sets to the CD09 background
# This is important because Fisher's exact test should use the same universe
# of genes that were actually tested in the CD09 edgeR analysis.
peak_npc_genes_cd09_bg <- intersect(
  peak_npc_genes,
  cd09_background_genes
)

peak_mcclay_genes_cd09_bg <- intersect(
  peak_mcclay_genes,
  cd09_background_genes
)

peak_forrest_genes_cd09_bg <- intersect(
  peak_forrest_genes,
  cd09_background_genes
)

# Check peak-associated genes present in CD09 background
data.frame(
  TCF4_Dataset = c(
    "NPC_peaks",
    "McClay_peaks",
    "Forrest_peaks"
  ),
  Peak_Genes_in_CD09_Background = c(
    length(peak_npc_genes_cd09_bg),
    length(peak_mcclay_genes_cd09_bg),
    length(peak_forrest_genes_cd09_bg)
  )
)


# Overall overlap: TCF4 peak genes vs CD09 significant DEGs
peakNPC_cd09_sig_overlap <- intersect(
  peak_npc_genes_cd09_bg,
  cd09_sig_genes
)

peakMcClay_cd09_sig_overlap <- intersect(
  peak_mcclay_genes_cd09_bg,
  cd09_sig_genes
)

peakForrest_cd09_sig_overlap <- intersect(
  peak_forrest_genes_cd09_bg,
  cd09_sig_genes
)

# Check raw overlap sizes
data.frame(
  TCF4_Dataset = c(
    "NPC_peaks",
    "McClay_peaks",
    "Forrest_peaks"
  ),
  Overlap_with_CD09_DEGs = c(
    length(peakNPC_cd09_sig_overlap),
    length(peakMcClay_cd09_sig_overlap),
    length(peakForrest_cd09_sig_overlap)
  )
)


# Fisher exact test: NPC peaks vs CD09 DEGs
# A = genes that are both CD09 DEGs and NPC TCF4 peak-associated genes
peakNPC_cd09_overlap <- intersect(
  peak_npc_genes_cd09_bg,
  cd09_sig_genes
)

# B = CD09 DEGs that are not NPC TCF4 peak genes
peakNPC_cd09_deg_only <- setdiff(
  cd09_sig_genes,
  peak_npc_genes_cd09_bg
)

# C = NPC TCF4 peak genes that are not CD09 DEGs
peakNPC_cd09_peak_only <- setdiff(
  peak_npc_genes_cd09_bg,
  cd09_sig_genes
)

# D = background genes that are neither
peakNPC_cd09_neither <- setdiff(
  cd09_background_genes,
  union(peak_npc_genes_cd09_bg, cd09_sig_genes)
)

peakNPC_cd09_fisher_table <- matrix(
  c(
    length(peakNPC_cd09_overlap),
    length(peakNPC_cd09_deg_only),
    length(peakNPC_cd09_peak_only),
    length(peakNPC_cd09_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(peakNPC_cd09_fisher_table) <- c(
  "CD09_DEG",
  "Not_CD09_DEG"
)

colnames(peakNPC_cd09_fisher_table) <- c(
  "NPC_TCF4_peak_gene",
  "No_NPC_TCF4_peak"
)

peakNPC_cd09_fisher_table

peakNPC_cd09_fisher <- fisher.test(
  peakNPC_cd09_fisher_table
)

peakNPC_cd09_fisher


# Fisher exact test: McClay peaks vs CD09 DEGs
peakMcClay_cd09_overlap <- intersect(
  peak_mcclay_genes_cd09_bg,
  cd09_sig_genes
)

peakMcClay_cd09_deg_only <- setdiff(
  cd09_sig_genes,
  peak_mcclay_genes_cd09_bg
)

peakMcClay_cd09_peak_only <- setdiff(
  peak_mcclay_genes_cd09_bg,
  cd09_sig_genes
)

peakMcClay_cd09_neither <- setdiff(
  cd09_background_genes,
  union(peak_mcclay_genes_cd09_bg, cd09_sig_genes)
)

peakMcClay_cd09_fisher_table <- matrix(
  c(
    length(peakMcClay_cd09_overlap),
    length(peakMcClay_cd09_deg_only),
    length(peakMcClay_cd09_peak_only),
    length(peakMcClay_cd09_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(peakMcClay_cd09_fisher_table) <- c(
  "CD09_DEG",
  "Not_CD09_DEG"
)

colnames(peakMcClay_cd09_fisher_table) <- c(
  "McClay_TCF4_peak_gene",
  "No_McClay_TCF4_peak"
)

peakMcClay_cd09_fisher_table

peakMcClay_cd09_fisher <- fisher.test(
  peakMcClay_cd09_fisher_table
)

peakMcClay_cd09_fisher


# Fisher exact test: Forrest peaks vs CD09 DEGs
peakForrest_cd09_overlap <- intersect(
  peak_forrest_genes_cd09_bg,
  cd09_sig_genes
)

peakForrest_cd09_deg_only <- setdiff(
  cd09_sig_genes,
  peak_forrest_genes_cd09_bg
)

peakForrest_cd09_peak_only <- setdiff(
  peak_forrest_genes_cd09_bg,
  cd09_sig_genes
)

peakForrest_cd09_neither <- setdiff(
  cd09_background_genes,
  union(peak_forrest_genes_cd09_bg, cd09_sig_genes)
)

peakForrest_cd09_fisher_table <- matrix(
  c(
    length(peakForrest_cd09_overlap),
    length(peakForrest_cd09_deg_only),
    length(peakForrest_cd09_peak_only),
    length(peakForrest_cd09_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(peakForrest_cd09_fisher_table) <- c(
  "CD09_DEG",
  "Not_CD09_DEG"
)

colnames(peakForrest_cd09_fisher_table) <- c(
  "Forrest_TCF4_peak_gene",
  "No_Forrest_TCF4_peak"
)

peakForrest_cd09_fisher_table

peakForrest_cd09_fisher <- fisher.test(
  peakForrest_cd09_fisher_table
)

peakForrest_cd09_fisher


# Summary table: CD09 DEGs vs TCF4 peak-associated genes
tcf4_cd09_fisher_summary <- data.frame(

  Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Background_Genes = c(
    length(cd09_background_genes),
    length(cd09_background_genes),
    length(cd09_background_genes)
  ),

  CD09_DEGs = c(
    length(cd09_sig_genes),
    length(cd09_sig_genes),
    length(cd09_sig_genes)
  ),

  TCF4_Peak_Genes_In_Background = c(
    length(peak_npc_genes_cd09_bg),
    length(peak_mcclay_genes_cd09_bg),
    length(peak_forrest_genes_cd09_bg)
  ),

  Overlap_Genes = c(
    length(peakNPC_cd09_overlap),
    length(peakMcClay_cd09_overlap),
    length(peakForrest_cd09_overlap)
  ),

  DEG_Only_Genes = c(
    length(peakNPC_cd09_deg_only),
    length(peakMcClay_cd09_deg_only),
    length(peakForrest_cd09_deg_only)
  ),

  Peak_Only_Genes = c(
    length(peakNPC_cd09_peak_only),
    length(peakMcClay_cd09_peak_only),
    length(peakForrest_cd09_peak_only)
  ),

  Neither_Genes = c(
    length(peakNPC_cd09_neither),
    length(peakMcClay_cd09_neither),
    length(peakForrest_cd09_neither)
  ),

  Odds_Ratio = c(
    as.numeric(peakNPC_cd09_fisher$estimate),
    as.numeric(peakMcClay_cd09_fisher$estimate),
    as.numeric(peakForrest_cd09_fisher$estimate)
  ),

  P_Value = c(
    peakNPC_cd09_fisher$p.value,
    peakMcClay_cd09_fisher$p.value,
    peakForrest_cd09_fisher$p.value
  )
)

tcf4_cd09_fisher_summary


# Upregulated and downregulated CD09 DEG overlaps
## NPC peaks
peakNPC_cd09_up_overlap <- intersect(
  peak_npc_genes_cd09_bg,
  cd09_up_genes
)

peakNPC_cd09_down_overlap <- intersect(
  peak_npc_genes_cd09_bg,
  cd09_down_genes
)

## McClay peaks
peakMcClay_cd09_up_overlap <- intersect(
  peak_mcclay_genes_cd09_bg,
  cd09_up_genes
)

peakMcClay_cd09_down_overlap <- intersect(
  peak_mcclay_genes_cd09_bg,
  cd09_down_genes
)

## Forrest peaks
peakForrest_cd09_up_overlap <- intersect(
  peak_forrest_genes_cd09_bg,
  cd09_up_genes
)

peakForrest_cd09_down_overlap <- intersect(
  peak_forrest_genes_cd09_bg,
  cd09_down_genes
)


# Fisher tests for UP and DOWN CD09 DEGs
## NPC peaks vs CD09 UP
peakNPC_cd09_up_fisher <- fisher.test(
  matrix(
    c(
      length(peakNPC_cd09_up_overlap),
      length(setdiff(cd09_up_genes, peak_npc_genes_cd09_bg)),
      length(setdiff(peak_npc_genes_cd09_bg, cd09_up_genes)),
      length(setdiff(
        cd09_background_genes,
        union(peak_npc_genes_cd09_bg, cd09_up_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

## NPC peaks vs CD09 DOWN
peakNPC_cd09_down_fisher <- fisher.test(
  matrix(
    c(
      length(peakNPC_cd09_down_overlap),
      length(setdiff(cd09_down_genes, peak_npc_genes_cd09_bg)),
      length(setdiff(peak_npc_genes_cd09_bg, cd09_down_genes)),
      length(setdiff(
        cd09_background_genes,
        union(peak_npc_genes_cd09_bg, cd09_down_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

## McClay peaks vs CD09 UP
peakMcClay_cd09_up_fisher <- fisher.test(
  matrix(
    c(
      length(peakMcClay_cd09_up_overlap),
      length(setdiff(cd09_up_genes, peak_mcclay_genes_cd09_bg)),
      length(setdiff(peak_mcclay_genes_cd09_bg, cd09_up_genes)),
      length(setdiff(
        cd09_background_genes,
        union(peak_mcclay_genes_cd09_bg, cd09_up_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

## McClay peaks vs CD09 DOWN
peakMcClay_cd09_down_fisher <- fisher.test(
  matrix(
    c(
      length(peakMcClay_cd09_down_overlap),
      length(setdiff(cd09_down_genes, peak_mcclay_genes_cd09_bg)),
      length(setdiff(peak_mcclay_genes_cd09_bg, cd09_down_genes)),
      length(setdiff(
        cd09_background_genes,
        union(peak_mcclay_genes_cd09_bg, cd09_down_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

## Forrest peaks vs CD09 UP
peakForrest_cd09_up_fisher <- fisher.test(
  matrix(
    c(
      length(peakForrest_cd09_up_overlap),
      length(setdiff(cd09_up_genes, peak_forrest_genes_cd09_bg)),
      length(setdiff(peak_forrest_genes_cd09_bg, cd09_up_genes)),
      length(setdiff(
        cd09_background_genes,
        union(peak_forrest_genes_cd09_bg, cd09_up_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

## Forrest peaks vs CD09 DOWN
peakForrest_cd09_down_fisher <- fisher.test(
  matrix(
    c(
      length(peakForrest_cd09_down_overlap),
      length(setdiff(cd09_down_genes, peak_forrest_genes_cd09_bg)),
      length(setdiff(peak_forrest_genes_cd09_bg, cd09_down_genes)),
      length(setdiff(
        cd09_background_genes,
        union(peak_forrest_genes_cd09_bg, cd09_down_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

# Summary table: UP and DOWN enrichment
tcf4_cd09_up_down_summary <- data.frame(

  Dataset = c(
    "NPC",
    "NPC",
    "McClay",
    "McClay",
    "Forrest",
    "Forrest"
  ),

  Direction = c(
    "Up",
    "Down",
    "Up",
    "Down",
    "Up",
    "Down"
  ),

  Overlap_Genes = c(
    length(peakNPC_cd09_up_overlap),
    length(peakNPC_cd09_down_overlap),
    length(peakMcClay_cd09_up_overlap),
    length(peakMcClay_cd09_down_overlap),
    length(peakForrest_cd09_up_overlap),
    length(peakForrest_cd09_down_overlap)
  ),

  Odds_Ratio = c(
    as.numeric(peakNPC_cd09_up_fisher$estimate),
    as.numeric(peakNPC_cd09_down_fisher$estimate),
    as.numeric(peakMcClay_cd09_up_fisher$estimate),
    as.numeric(peakMcClay_cd09_down_fisher$estimate),
    as.numeric(peakForrest_cd09_up_fisher$estimate),
    as.numeric(peakForrest_cd09_down_fisher$estimate)
  ),

  P_Value = c(
    peakNPC_cd09_up_fisher$p.value,
    peakNPC_cd09_down_fisher$p.value,
    peakMcClay_cd09_up_fisher$p.value,
    peakMcClay_cd09_down_fisher$p.value,
    peakForrest_cd09_up_fisher$p.value,
    peakForrest_cd09_down_fisher$p.value
  )
)

tcf4_cd09_up_down_summary



# Save CD09 overlap and Fisher results
dir.create(
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  tcf4_cd09_fisher_summary,
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_TCF4_overall_fisher_summary.csv",
  row.names = FALSE
)

write.csv(
  tcf4_cd09_up_down_summary,
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_TCF4_up_down_fisher_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakNPC_cd09_sig_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakNPC_all_overlap_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakMcClay_cd09_sig_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakMcClay_all_overlap_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakForrest_cd09_sig_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakForrest_all_overlap_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakNPC_cd09_up_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakNPC_up_overlap_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakNPC_cd09_down_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakNPC_down_overlap_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakMcClay_cd09_up_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakMcClay_up_overlap_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakMcClay_cd09_down_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakMcClay_down_overlap_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakForrest_cd09_up_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakForrest_up_overlap_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = peakForrest_cd09_down_overlap),
  "Decon_Doost/CD09_edgeR_LRT_TCF4_overlap/CD09_peakForrest_down_overlap_genes.csv",
  row.names = FALSE
)
