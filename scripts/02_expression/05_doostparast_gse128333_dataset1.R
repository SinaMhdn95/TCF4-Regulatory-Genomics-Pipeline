####################################################################
# Doosparats Deconvolution of Transcriptional Networks Identifies
# TCF4 as a Master Regulator in Schizophrenia
# Using KD_raw1 only
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

txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Import files
npc_peaks     <- import("NPC_hg38_summit_500bp.bed")
mcclay_peaks  <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
forrest_peaks <- import("Forrest_hg38.bed")

# Import RNA-seq TCF4 KD data
KD_raw1  <- read.table(
  "GSE128333_Data1.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

KD_raw2 <- read.table(
  "GSE128333_Data2.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

head(KD_raw1)
head(KD_raw2)

# STEP1: NPC Analysis (only using KD_raw1)
# Check the raw count matrix
dim(KD_raw1)
head(KD_raw1)
colnames(KD_raw1)

################################################################################
# Select only Day 3 samples
# Day 3 samples = hiPSC-derived neural progenitor cells (NPCs)
# Ctrl = control
# K = TCF4 knockdown
npc_counts_3d <- KD_raw1[, c(
  "D3_Ctrl1",
  "D3_Ctrl2",
  "D3_Ctrl3",
  "D3_K1",
  "D3_K2",
  "D3_K3"
)]

# Make sure counts are numeric integers
npc_counts_3d <- round(as.matrix(npc_counts_3d))

# Create sample information table
# This tells DESeq2 which samples are control and which are knockdown
npc_coldata <- data.frame(
  condition = c(
    "Ctrl",
    "Ctrl",
    "Ctrl",
    "KD",
    "KD",
    "KD"
  )
)

# Row names of coldata must exactly match column names of count matrix
rownames(npc_coldata) <- colnames(npc_counts_3d)

# Check that sample names match
all(rownames(npc_coldata) == colnames(npc_counts_3d))

# Set Ctrl as the reference group
# This means log2FoldChange will be KD compared with Ctrl
npc_coldata$condition <- factor(
  npc_coldata$condition,
  levels = c("Ctrl", "KD")
)

# Create DESeq2 object
dds_npc_3d <- DESeqDataSetFromMatrix(
  countData = npc_counts_3d,
  colData = npc_coldata,
  design = ~ condition
)

# Filter low-count genes
# Keeps genes with total count >= 10 across all NPC samples
dds_npc_3d <- dds_npc_3d[rowSums(counts(dds_npc_3d)) >= 10, ]

## Run differential expression analysis
dds_npc_3d <- DESeq(dds_npc_3d)

# Extract results
# This compares KD vs Ctrl
npc_res_3d <- results(
  dds_npc_3d,
  contrast = c("condition", "KD", "Ctrl")
)

## Convert results to data frame
npc_res_3d_df <- as.data.frame(npc_res_3d)

## Add Ensembl IDs as a column
npc_res_3d_df$ENSEMBL <- rownames(npc_res_3d_df)

## Sort by adjusted p-value
npc_res_3d_df <- npc_res_3d_df[order(npc_res_3d_df$padj), ]

## Check top genes
head(npc_res_3d_df)

## Save full results
write.csv(
  npc_res_3d_df,
  "Decon_Doost/GSE128333_NPC_TCF4_KD_vs_Ctrl_DESeq2_full_results.csv",
  row.names = FALSE
)

# Significant DEGs
npc_3d_sig <- subset(
  npc_res_3d_df,
  padj < 0.05
)

# Upregulated after TCF4 knockdown
npc_3d_up <- subset(
  npc_res_3d_df,
  padj < 0.05 & log2FoldChange > 0
)

## Downregulated after TCF4 knockdown
npc_3d_down <- subset(
  npc_res_3d_df,
  padj < 0.05 & log2FoldChange < 0
)

## Count significant genes
nrow(npc_3d_sig)
nrow(npc_3d_up)
nrow(npc_3d_down)

## Save significant gene lists
write.csv(npc_3d_sig, "GSE128333_NPC_TCF4_KD_significant_DEGs.csv", row.names = FALSE)
write.csv(npc_3d_up, "GSE128333_NPC_TCF4_KD_upregulated_DEGs.csv", row.names = FALSE)
write.csv(npc_3d_down, "GSE128333_NPC_TCF4_KD_downregulated_DEGs.csv", row.names = FALSE)


# Let's go back to our TCF4 datasets
# Annotate all three peak sets
# Annotate NPC peaks
npc_anno <- annotatePeak(
  npc_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

npc_genes <- unique(
  na.omit(as.data.frame(npc_anno)$SYMBOL)
)

length(npc_genes)


# Annotate McClay peaks
mcclay_anno <- annotatePeak(
  mcclay_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

mcclay_genes <- unique(
  na.omit(as.data.frame(mcclay_anno)$SYMBOL)
)

length(mcclay_genes)

# Annotate Forrest peaks
forrest_anno <- annotatePeak(
  forrest_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

forrest_genes <- unique(
  na.omit(as.data.frame(forrest_anno)$SYMBOL)
)

length(forrest_genes)


# Overlap and Fisher Test
# Background genes = all genes tested in DESeq2
## Remove Ensembl version numbers if present
## Your KD_raw1 IDs do not seem to have version numbers, but this is safe
npc_res_3d_df$ENSEMBL_clean <- sub(
  "\\..*",
  "",
  npc_res_3d_df$ENSEMBL
)

# Map Ensembl IDs to gene symbols
npc_res_3d_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = npc_res_3d_df$ENSEMBL_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Check result
head(npc_res_3d_df)

# How many genes mapped?
sum(!is.na(npc_res_3d_df$SYMBOL))

# Remove genes without symbol
npc_res_3d_df_clean <- npc_res_3d_df[
  !is.na(npc_res_3d_df$SYMBOL),
]

# Check clean table
head(npc_res_3d_df_clean)


# This includes significant and non-significant genes
kd_background_genes <- unique(npc_res_3d_df_clean$SYMBOL)

# Significant genes after TCF4 knockdown
kd_sig_genes <- unique(
  npc_res_3d_df_clean$SYMBOL[
    npc_res_3d_df_clean$padj < 0.05
  ]
)

# Restrict each TCF4 gene set to the DESeq2 background
npc_genes_bg <- intersect(npc_genes, kd_background_genes)
mcclay_genes_bg <- intersect(mcclay_genes, kd_background_genes)
forrest_genes_bg <- intersect(forrest_genes, kd_background_genes)


# Check gene set sizes
length(kd_background_genes)
length(kd_sig_genes)
length(npc_genes_bg)
length(mcclay_genes_bg)
length(forrest_genes_bg)


# Check Raw overlaps
npc_overlap_genes <- intersect(
  npc_genes_bg,
  kd_sig_genes
)

mcclay_overlap_genes <- intersect(
  mcclay_genes_bg,
  kd_sig_genes
)

forrest_overlap_genes <- intersect(
  forrest_genes_bg,
  kd_sig_genes
)

length(npc_overlap_genes)
length(mcclay_overlap_genes)
length(forrest_overlap_genes)

###### Fisher Test for NPC Peaks vs TCF4 Knockdown DEGs
# Find Genes that are both:
# 1) associated with NPC TCF4 peaks
# 2) significant after TCF4 knockdown
npc_overlap_genes <- intersect(
  npc_genes_bg,
  kd_sig_genes
)

# Find Genes that are significant after TCF4 knockdown
# but are NOT associated with NPC TCF4 peaks
npc_kd_only_genes <- setdiff(
  kd_sig_genes,
  npc_genes_bg
)

# Find Genes associated with NPC TCF4 peaks
# but are NOT significant after TCF4 knockdown
npc_peak_only_genes <- setdiff(
  npc_genes_bg,
  kd_sig_genes
)

# Background genes that are neither:
# 1) TCF4 peak-associated
# 2) significant after TCF4 knockdown
npc_neither_genes <- setdiff(
  kd_background_genes,
  union(npc_genes_bg, kd_sig_genes)
)

# Create Fisher table using counts
npc_fisher_table <- matrix(
  c(
    length(npc_overlap_genes),
    length(npc_kd_only_genes),
    length(npc_peak_only_genes),
    length(npc_neither_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(npc_fisher_table) <- c(
  "KD_DEG",
  "Not_KD_DEG"
)

colnames(npc_fisher_table) <- c(
  "NPC_TCF4_peak_gene",
  "No_NPC_TCF4_peak"
)

npc_fisher_table

# Run Fisher exact test
npc_fisher <- fisher.test(npc_fisher_table)

npc_fisher

##### Fisher Test for McClay Peaks vs TCF4 Knockdown DEGs
mcclay_overlap_genes <- intersect(
  mcclay_genes_bg,
  kd_sig_genes
)

mcclay_kd_only_genes <- setdiff(
  kd_sig_genes,
  mcclay_genes_bg
)

mcclay_peak_only_genes <- setdiff(
  mcclay_genes_bg,
  kd_sig_genes
)

mcclay_neither_genes <- setdiff(
  kd_background_genes,
  union(mcclay_genes_bg, kd_sig_genes)
)

mcclay_fisher_table <- matrix(
  c(
    length(mcclay_overlap_genes),
    length(mcclay_kd_only_genes),
    length(mcclay_peak_only_genes),
    length(mcclay_neither_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(mcclay_fisher_table) <- c(
  "KD_DEG",
  "Not_KD_DEG"
)

colnames(mcclay_fisher_table) <- c(
  "McClay_TCF4_peak_gene",
  "No_McClay_TCF4_peak"
)

mcclay_fisher_table

mcclay_fisher <- fisher.test(mcclay_fisher_table)

mcclay_fisher


##### Fisher Test for Forrest Peaks vs TCF4 Knockdown DEGs
forrest_overlap_genes <- intersect(
  forrest_genes_bg,
  kd_sig_genes
)

forrest_kd_only_genes <- setdiff(
  kd_sig_genes,
  forrest_genes_bg
)

forrest_peak_only_genes <- setdiff(
  forrest_genes_bg,
  kd_sig_genes
)

forrest_neither_genes <- setdiff(
  kd_background_genes,
  union(forrest_genes_bg, kd_sig_genes)
)

forrest_fisher_table <- matrix(
  c(
    length(forrest_overlap_genes),
    length(forrest_kd_only_genes),
    length(forrest_peak_only_genes),
    length(forrest_neither_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(forrest_fisher_table) <- c(
  "KD_DEG",
  "Not_KD_DEG"
)

colnames(forrest_fisher_table) <- c(
  "Forrest_TCF4_peak_gene",
  "No_Forrest_TCF4_peak"
)

forrest_fisher_table

forrest_fisher <- fisher.test(forrest_fisher_table)

forrest_fisher


# Summary Table for Fisher Enrichment Results
tcf4_kd_fisher_summary <- data.frame(
  Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Background_Genes = c(
    length(kd_background_genes),
    length(kd_background_genes),
    length(kd_background_genes)
  ),

  KD_DEGs = c(
    length(kd_sig_genes),
    length(kd_sig_genes),
    length(kd_sig_genes)
  ),

  TCF4_Peak_Genes_In_Background = c(
    length(npc_genes_bg),
    length(mcclay_genes_bg),
    length(forrest_genes_bg)
  ),

  Overlap_Genes = c(
    length(npc_overlap_genes),
    length(mcclay_overlap_genes),
    length(forrest_overlap_genes)
  ),

  KD_Only_Genes = c(
    length(npc_kd_only_genes),
    length(mcclay_kd_only_genes),
    length(forrest_kd_only_genes)
  ),

  Peak_Only_Genes = c(
    length(npc_peak_only_genes),
    length(mcclay_peak_only_genes),
    length(forrest_peak_only_genes)
  ),

  Neither_Genes = c(
    length(npc_neither_genes),
    length(mcclay_neither_genes),
    length(forrest_neither_genes)
  ),

  Odds_Ratio = c(
    as.numeric(npc_fisher$estimate),
    as.numeric(mcclay_fisher$estimate),
    as.numeric(forrest_fisher$estimate)
  ),

  P_Value = c(
    npc_fisher$p.value,
    mcclay_fisher$p.value,
    forrest_fisher$p.value
  )
)

tcf4_kd_fisher_summary



# Create Upregulated and Downregulated Gene Symbol Lists
npc_3d_sig <- subset(
  npc_res_3d_df_clean,
  padj < 0.05
)

npc_3d_up <- subset(
  npc_res_3d_df_clean,
  padj < 0.05 & log2FoldChange > 0
)

npc_3d_down <- subset(
  npc_res_3d_df_clean,
  padj < 0.05 & log2FoldChange < 0
)

# Extract gene symbols
kd_sig_genes <- unique(npc_3d_sig$SYMBOL)

kd_up_genes <- unique(npc_3d_up$SYMBOL)

kd_down_genes <- unique(npc_3d_down$SYMBOL)

# Check numbers
length(kd_sig_genes)
length(kd_up_genes)
length(kd_down_genes)


# NPC Peaks vs UPREGULATED genes
npc_up_overlap_genes <- intersect(
  npc_genes_bg,
  kd_up_genes
)

npc_up_only_genes <- setdiff(
  kd_up_genes,
  npc_genes_bg
)

npc_peak_not_up_genes <- setdiff(
  npc_genes_bg,
  kd_up_genes
)

npc_neither_up_genes <- setdiff(
  kd_background_genes,
  union(npc_genes_bg, kd_up_genes)
)

npc_up_fisher_table <- matrix(
  c(
    length(npc_up_overlap_genes),
    length(npc_up_only_genes),
    length(npc_peak_not_up_genes),
    length(npc_neither_up_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

npc_up_fisher <- fisher.test(npc_up_fisher_table)

npc_up_fisher


# NPC Peaks vs DOWNREGULATED genes
npc_down_overlap_genes <- intersect(
  npc_genes_bg,
  kd_down_genes
)

npc_down_only_genes <- setdiff(
  kd_down_genes,
  npc_genes_bg
)

npc_peak_not_down_genes <- setdiff(
  npc_genes_bg,
  kd_down_genes
)

npc_neither_down_genes <- setdiff(
  kd_background_genes,
  union(npc_genes_bg, kd_down_genes)
)

npc_down_fisher_table <- matrix(
  c(
    length(npc_down_overlap_genes),
    length(npc_down_only_genes),
    length(npc_peak_not_down_genes),
    length(npc_neither_down_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

npc_down_fisher <- fisher.test(npc_down_fisher_table)

npc_down_fisher


# McClay Peaks vs Upregulated Genes
mcclay_up_overlap_genes <- intersect(
  mcclay_genes_bg,
  kd_up_genes
)

mcclay_up_only_genes <- setdiff(
  kd_up_genes,
  mcclay_genes_bg
)

mcclay_peak_not_up_genes <- setdiff(
  mcclay_genes_bg,
  kd_up_genes
)

mcclay_neither_up_genes <- setdiff(
  kd_background_genes,
  union(mcclay_genes_bg, kd_up_genes)
)

mcclay_up_fisher_table <- matrix(
  c(
    length(mcclay_up_overlap_genes),
    length(mcclay_up_only_genes),
    length(mcclay_peak_not_up_genes),
    length(mcclay_neither_up_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

mcclay_up_fisher <- fisher.test(mcclay_up_fisher_table)


# McClay Peaks vs Downregulated Genes
mcclay_down_overlap_genes <- intersect(
  mcclay_genes_bg,
  kd_down_genes
)

mcclay_down_only_genes <- setdiff(
  kd_down_genes,
  mcclay_genes_bg
)

mcclay_peak_not_down_genes <- setdiff(
  mcclay_genes_bg,
  kd_down_genes
)

mcclay_neither_down_genes <- setdiff(
  kd_background_genes,
  union(mcclay_genes_bg, kd_down_genes)
)

mcclay_down_fisher_table <- matrix(
  c(
    length(mcclay_down_overlap_genes),
    length(mcclay_down_only_genes),
    length(mcclay_peak_not_down_genes),
    length(mcclay_neither_down_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

mcclay_down_fisher <- fisher.test(mcclay_down_fisher_table)


# Forrest Peaks vs Upregulated Genes
forrest_up_overlap_genes <- intersect(
  forrest_genes_bg,
  kd_up_genes
)

forrest_up_only_genes <- setdiff(
  kd_up_genes,
  forrest_genes_bg
)

forrest_peak_not_up_genes <- setdiff(
  forrest_genes_bg,
  kd_up_genes
)

forrest_neither_up_genes <- setdiff(
  kd_background_genes,
  union(forrest_genes_bg, kd_up_genes)
)

forrest_up_fisher_table <- matrix(
  c(
    length(forrest_up_overlap_genes),
    length(forrest_up_only_genes),
    length(forrest_peak_not_up_genes),
    length(forrest_neither_up_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

forrest_up_fisher <- fisher.test(forrest_up_fisher_table)


# Forrest Peaks vs Downregulated Genes
forrest_down_overlap_genes <- intersect(
  forrest_genes_bg,
  kd_down_genes
)

forrest_down_only_genes <- setdiff(
  kd_down_genes,
  forrest_genes_bg
)

forrest_peak_not_down_genes <- setdiff(
  forrest_genes_bg,
  kd_down_genes
)

forrest_neither_down_genes <- setdiff(
  kd_background_genes,
  union(forrest_genes_bg, kd_down_genes)
)

forrest_down_fisher_table <- matrix(
  c(
    length(forrest_down_overlap_genes),
    length(forrest_down_only_genes),
    length(forrest_peak_not_down_genes),
    length(forrest_neither_down_genes)
  ),
  nrow = 2,
  byrow = TRUE
)

forrest_down_fisher <- fisher.test(forrest_down_fisher_table)


# Compare Upregulated and Downregulated Enrichment
tcf4_up_down_comparison <- data.frame(

  Dataset = c(
    "NPC",
    "NPC",
    "McClay",
    "McClay",
    "Forrest",
    "Forrest"
  ),

  DEG_Direction = c(
    "Up",
    "Down",
    "Up",
    "Down",
    "Up",
    "Down"
  ),

  Overlap_Genes = c(
    length(npc_up_overlap_genes),
    length(npc_down_overlap_genes),

    length(mcclay_up_overlap_genes),
    length(mcclay_down_overlap_genes),

    length(forrest_up_overlap_genes),
    length(forrest_down_overlap_genes)
  ),

  Odds_Ratio = c(
    as.numeric(npc_up_fisher$estimate),
    as.numeric(npc_down_fisher$estimate),

    as.numeric(mcclay_up_fisher$estimate),
    as.numeric(mcclay_down_fisher$estimate),

    as.numeric(forrest_up_fisher$estimate),
    as.numeric(forrest_down_fisher$estimate)
  ),

  P_Value = c(
    npc_up_fisher$p.value,
    npc_down_fisher$p.value,

    mcclay_up_fisher$p.value,
    mcclay_down_fisher$p.value,

    forrest_up_fisher$p.value,
    forrest_down_fisher$p.value
  )
)

tcf4_up_down_comparison


# Genes shared by all 4 datasets
shared_4way_genes <- Reduce(
  intersect,
  list(
    npc_genes_bg,
    mcclay_genes_bg,
    forrest_genes_bg,
    kd_sig_genes
  )
)

length(shared_4way_genes)

head(shared_4way_genes)

#For the genes shared across all 4 datasets:
# Downregulated after KD
shared_4way_down <- intersect(
  shared_4way_genes,
  kd_down_genes
)

length(shared_4way_down)

# Upregulated after KD
shared_4way_up <- intersect(
  shared_4way_genes,
  kd_up_genes
)

length(shared_4way_up)


#### JUST FOR FUN #######
# Summary of Up vs Down Overlaps Across TCF4 Datasets

tcf4_overlap_direction_summary <- data.frame(

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
    length(npc_up_overlap_genes),
    length(npc_down_overlap_genes),

    length(mcclay_up_overlap_genes),
    length(mcclay_down_overlap_genes),

    length(forrest_up_overlap_genes),
    length(forrest_down_overlap_genes)
  ),

  Odds_Ratio = c(
    as.numeric(npc_up_fisher$estimate),
    as.numeric(npc_down_fisher$estimate),

    as.numeric(mcclay_up_fisher$estimate),
    as.numeric(mcclay_down_fisher$estimate),

    as.numeric(forrest_up_fisher$estimate),
    as.numeric(forrest_down_fisher$estimate)
  ),

  P_Value = c(
    npc_up_fisher$p.value,
    npc_down_fisher$p.value,

    mcclay_up_fisher$p.value,
    mcclay_down_fisher$p.value,

    forrest_up_fisher$p.value,
    forrest_down_fisher$p.value
  )
)

tcf4_overlap_direction_summary

#### Shared genes between datasets
# Shared DOWN targets

shared_down_all <- Reduce(
  intersect,
  list(
    npc_down_overlap_genes,
    mcclay_down_overlap_genes,
    forrest_down_overlap_genes
  )
)

length(shared_down_all)


# Shared UP targets
shared_up_all <- Reduce(
  intersect,
  list(
    npc_up_overlap_genes,
    mcclay_up_overlap_genes,
    forrest_up_overlap_genes
  )
)

length(shared_up_all)

# NPC vs McClay
npcUP_vs_mcclayUP <- intersect(
  npc_up_overlap_genes,
  mcclay_up_overlap_genes
)

npcUP_vs_mcclayDOWN <- intersect(
  npc_up_overlap_genes,
  mcclay_down_overlap_genes
)

npcDOWN_vs_mcclayUP <- intersect(
  npc_down_overlap_genes,
  mcclay_up_overlap_genes
)

npcDOWN_vs_mcclayDOWN <- intersect(
  npc_down_overlap_genes,
  mcclay_down_overlap_genes
)


# NPC vs Forrest
npcUP_vs_forrestUP <- intersect(
  npc_up_overlap_genes,
  forrest_up_overlap_genes
)

npcUP_vs_forrestDOWN <- intersect(
  npc_up_overlap_genes,
  forrest_down_overlap_genes
)

npcDOWN_vs_forrestUP <- intersect(
  npc_down_overlap_genes,
  forrest_up_overlap_genes
)

npcDOWN_vs_forrestDOWN <- intersect(
  npc_down_overlap_genes,
  forrest_down_overlap_genes
)


# McClay vs Forrest
mcclayUP_vs_forrestUP <- intersect(
  mcclay_up_overlap_genes,
  forrest_up_overlap_genes
)

mcclayUP_vs_forrestDOWN <- intersect(
  mcclay_up_overlap_genes,
  forrest_down_overlap_genes
)

mcclayDOWN_vs_forrestUP <- intersect(
  mcclay_down_overlap_genes,
  forrest_up_overlap_genes
)

mcclayDOWN_vs_forrestDOWN <- intersect(
  mcclay_down_overlap_genes,
  forrest_down_overlap_genes
)


# Pairwise Overlap Summary
pairwise_overlap_summary <- data.frame(

  Comparison = c(
    "NPC_UP_vs_McClay_UP",
    "NPC_UP_vs_McClay_DOWN",
    "NPC_DOWN_vs_McClay_UP",
    "NPC_DOWN_vs_McClay_DOWN",

    "NPC_UP_vs_Forrest_UP",
    "NPC_UP_vs_Forrest_DOWN",
    "NPC_DOWN_vs_Forrest_UP",
    "NPC_DOWN_vs_Forrest_DOWN",

    "McClay_UP_vs_Forrest_UP",
    "McClay_UP_vs_Forrest_DOWN",
    "McClay_DOWN_vs_Forrest_UP",
    "McClay_DOWN_vs_Forrest_DOWN"
  ),

  Shared_Genes = c(
    length(npcUP_vs_mcclayUP),
    length(npcUP_vs_mcclayDOWN),
    length(npcDOWN_vs_mcclayUP),
    length(npcDOWN_vs_mcclayDOWN),

    length(npcUP_vs_forrestUP),
    length(npcUP_vs_forrestDOWN),
    length(npcDOWN_vs_forrestUP),
    length(npcDOWN_vs_forrestDOWN),

    length(mcclayUP_vs_forrestUP),
    length(mcclayUP_vs_forrestDOWN),
    length(mcclayDOWN_vs_forrestUP),
    length(mcclayDOWN_vs_forrestDOWN)
  )
)

pairwise_overlap_summary[
  order(pairwise_overlap_summary$Shared_Genes,
        decreasing = TRUE),
]

####################################################################
# Day 14 Glutamatergic Neuron TCF4 Knockdown Analysis
# Dataset: GSE128333
# Goal: Compare glutamatergic DEGs with TCF4 peak-associated genes
####################################################################

# Select Day 14 glutamatergic samples
glut_counts_14d <- KD_raw1[, c(
  "D14_Ctrl1",
  "D14_Ctrl2",
  "D14_Ctrl3",
  "D14_K1",
  "D14_K2",
  "D14_K3"
)]

glut_counts_14d <- round(
  as.matrix(glut_counts_14d)
)


# Create sample metadata
glut_coldata <- data.frame(
  condition = c(
    "Ctrl",
    "Ctrl",
    "Ctrl",
    "KD",
    "KD",
    "KD"
  )
)

rownames(glut_coldata) <- colnames(glut_counts_14d)

glut_coldata$condition <- factor(
  glut_coldata$condition,
  levels = c("Ctrl", "KD")
)

all(rownames(glut_coldata) == colnames(glut_counts_14d))


# Run DESeq2
dds_glut_14d <- DESeqDataSetFromMatrix(
  countData = glut_counts_14d,
  colData = glut_coldata,
  design = ~ condition
)

dds_glut_14d <- dds_glut_14d[
  rowSums(counts(dds_glut_14d)) >= 10,
]

dds_glut_14d <- DESeq(
  dds_glut_14d
)

glut_res_14d <- results(
  dds_glut_14d,
  contrast = c("condition", "KD", "Ctrl")
)

glut_res_14d_df <- as.data.frame(glut_res_14d)

glut_res_14d_df$ENSEMBL <- rownames(glut_res_14d_df)

glut_res_14d_df <- glut_res_14d_df[
  order(glut_res_14d_df$padj),
]


# Annotate Ensembl IDs to gene symbols
glut_res_14d_df$ENSEMBL_clean <- sub(
  "\\..*",
  "",
  glut_res_14d_df$ENSEMBL
)

glut_res_14d_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = glut_res_14d_df$ENSEMBL_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

glut_res_14d_df_clean <- glut_res_14d_df[
  !is.na(glut_res_14d_df$SYMBOL),
]

head(glut_res_14d_df_clean)

# Create glut DEG gene lists
glut_sig <- subset(
  glut_res_14d_df_clean,
  padj < 0.05
)

glut_up <- subset(
  glut_res_14d_df_clean,
  padj < 0.05 & log2FoldChange > 0
)

glut_down <- subset(
  glut_res_14d_df_clean,
  padj < 0.05 & log2FoldChange < 0
)

glut_background_genes <- unique(
  glut_res_14d_df_clean$SYMBOL
)

glut_sig_genes <- unique(
  glut_sig$SYMBOL
)

glut_up_genes <- unique(
  glut_up$SYMBOL
)

glut_down_genes <- unique(
  glut_down$SYMBOL
)

data.frame(
  Category = c(
    "Background",
    "Significant_DEGs",
    "Upregulated",
    "Downregulated"
  ),
  Gene_Count = c(
    length(glut_background_genes),
    length(glut_sig_genes),
    length(glut_up_genes),
    length(glut_down_genes)
  )
)


# Restrict TCF4 peak-associated genes to glut background
peak_npc_genes_glut_bg <- intersect(
  npc_genes,
  glut_background_genes
)

peak_mcclay_genes_glut_bg <- intersect(
  mcclay_genes,
  glut_background_genes
)

peak_forrest_genes_glut_bg <- intersect(
  forrest_genes,
  glut_background_genes
)

data.frame(
  TCF4_Dataset = c(
    "NPC_peaks",
    "McClay_peaks",
    "Forrest_peaks"
  ),
  Peak_Genes_in_Glut_Background = c(
    length(peak_npc_genes_glut_bg),
    length(peak_mcclay_genes_glut_bg),
    length(peak_forrest_genes_glut_bg)
  )
)

# Overall overlaps: TCF4 peak genes vs glut DEGs
peakNPC_glut_sig_overlap <- intersect(
  peak_npc_genes_glut_bg,
  glut_sig_genes
)

peakMcClay_glut_sig_overlap <- intersect(
  peak_mcclay_genes_glut_bg,
  glut_sig_genes
)

peakForrest_glut_sig_overlap <- intersect(
  peak_forrest_genes_glut_bg,
  glut_sig_genes
)


# Fisher test: NPC peaks vs glut DEGs
peakNPC_glut_deg_only <- setdiff(
  glut_sig_genes,
  peak_npc_genes_glut_bg
)

peakNPC_glut_peak_only <- setdiff(
  peak_npc_genes_glut_bg,
  glut_sig_genes
)

peakNPC_glut_neither <- setdiff(
  glut_background_genes,
  union(peak_npc_genes_glut_bg, glut_sig_genes)
)

peakNPC_glut_fisher_table <- matrix(
  c(
    length(peakNPC_glut_sig_overlap),
    length(peakNPC_glut_deg_only),
    length(peakNPC_glut_peak_only),
    length(peakNPC_glut_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(peakNPC_glut_fisher_table) <- c(
  "Glut_DEG",
  "Not_Glut_DEG"
)

colnames(peakNPC_glut_fisher_table) <- c(
  "NPC_TCF4_peak_gene",
  "No_NPC_TCF4_peak"
)

peakNPC_glut_fisher <- fisher.test(
  peakNPC_glut_fisher_table
)

# Fisher test: McClay peaks vs glut DEGs
peakMcClay_glut_deg_only <- setdiff(
  glut_sig_genes,
  peak_mcclay_genes_glut_bg
)

peakMcClay_glut_peak_only <- setdiff(
  peak_mcclay_genes_glut_bg,
  glut_sig_genes
)

peakMcClay_glut_neither <- setdiff(
  glut_background_genes,
  union(peak_mcclay_genes_glut_bg, glut_sig_genes)
)

peakMcClay_glut_fisher_table <- matrix(
  c(
    length(peakMcClay_glut_sig_overlap),
    length(peakMcClay_glut_deg_only),
    length(peakMcClay_glut_peak_only),
    length(peakMcClay_glut_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(peakMcClay_glut_fisher_table) <- c(
  "Glut_DEG",
  "Not_Glut_DEG"
)

colnames(peakMcClay_glut_fisher_table) <- c(
  "McClay_TCF4_peak_gene",
  "No_McClay_TCF4_peak"
)

peakMcClay_glut_fisher <- fisher.test(
  peakMcClay_glut_fisher_table
)

# Fisher test: Forrest peaks vs glut DEGs
peakForrest_glut_deg_only <- setdiff(
  glut_sig_genes,
  peak_forrest_genes_glut_bg
)

peakForrest_glut_peak_only <- setdiff(
  peak_forrest_genes_glut_bg,
  glut_sig_genes
)

peakForrest_glut_neither <- setdiff(
  glut_background_genes,
  union(peak_forrest_genes_glut_bg, glut_sig_genes)
)

peakForrest_glut_fisher_table <- matrix(
  c(
    length(peakForrest_glut_sig_overlap),
    length(peakForrest_glut_deg_only),
    length(peakForrest_glut_peak_only),
    length(peakForrest_glut_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(peakForrest_glut_fisher_table) <- c(
  "Glut_DEG",
  "Not_Glut_DEG"
)

colnames(peakForrest_glut_fisher_table) <- c(
  "Forrest_TCF4_peak_gene",
  "No_Forrest_TCF4_peak"
)

peakForrest_glut_fisher <- fisher.test(
  peakForrest_glut_fisher_table
)

# Overall Fisher summary table
tcf4_glut_fisher_summary <- data.frame(
  Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Background_Genes = c(
    length(glut_background_genes),
    length(glut_background_genes),
    length(glut_background_genes)
  ),

  Glut_DEGs = c(
    length(glut_sig_genes),
    length(glut_sig_genes),
    length(glut_sig_genes)
  ),

  TCF4_Peak_Genes_In_Background = c(
    length(peak_npc_genes_glut_bg),
    length(peak_mcclay_genes_glut_bg),
    length(peak_forrest_genes_glut_bg)
  ),

  Overlap_Genes = c(
    length(peakNPC_glut_sig_overlap),
    length(peakMcClay_glut_sig_overlap),
    length(peakForrest_glut_sig_overlap)
  ),

  DEG_Only_Genes = c(
    length(peakNPC_glut_deg_only),
    length(peakMcClay_glut_deg_only),
    length(peakForrest_glut_deg_only)
  ),

  Peak_Only_Genes = c(
    length(peakNPC_glut_peak_only),
    length(peakMcClay_glut_peak_only),
    length(peakForrest_glut_peak_only)
  ),

  Neither_Genes = c(
    length(peakNPC_glut_neither),
    length(peakMcClay_glut_neither),
    length(peakForrest_glut_neither)
  ),

  Odds_Ratio = c(
    as.numeric(peakNPC_glut_fisher$estimate),
    as.numeric(peakMcClay_glut_fisher$estimate),
    as.numeric(peakForrest_glut_fisher$estimate)
  ),

  P_Value = c(
    peakNPC_glut_fisher$p.value,
    peakMcClay_glut_fisher$p.value,
    peakForrest_glut_fisher$p.value
  )
)

tcf4_glut_fisher_summary


# TCF4 vs GlutUP/GlutDOWN
# NPC Peaks vs Upregulated Glut DEGs
####################################################################
# TCF4 Peak Genes vs Glut UP / DOWN DEGs
####################################################################
## NPC peaks vs glut upregulated genes
peakNPC_glut_up_overlap <- intersect(
  peak_npc_genes_glut_bg,
  glut_up_genes
)

peakNPC_glut_down_overlap <- intersect(
  peak_npc_genes_glut_bg,
  glut_down_genes
)

peakNPC_glut_up_fisher <- fisher.test(
  matrix(
    c(
      length(peakNPC_glut_up_overlap),
      length(setdiff(glut_up_genes, peak_npc_genes_glut_bg)),
      length(setdiff(peak_npc_genes_glut_bg, glut_up_genes)),
      length(setdiff(
        glut_background_genes,
        union(peak_npc_genes_glut_bg, glut_up_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

peakNPC_glut_down_fisher <- fisher.test(
  matrix(
    c(
      length(peakNPC_glut_down_overlap),
      length(setdiff(glut_down_genes, peak_npc_genes_glut_bg)),
      length(setdiff(peak_npc_genes_glut_bg, glut_down_genes)),
      length(setdiff(
        glut_background_genes,
        union(peak_npc_genes_glut_bg, glut_down_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)


## McClay peaks vs glut up/down genes
peakMcClay_glut_up_overlap <- intersect(
  peak_mcclay_genes_glut_bg,
  glut_up_genes
)

peakMcClay_glut_down_overlap <- intersect(
  peak_mcclay_genes_glut_bg,
  glut_down_genes
)

peakMcClay_glut_up_fisher <- fisher.test(
  matrix(
    c(
      length(peakMcClay_glut_up_overlap),
      length(setdiff(glut_up_genes, peak_mcclay_genes_glut_bg)),
      length(setdiff(peak_mcclay_genes_glut_bg, glut_up_genes)),
      length(setdiff(
        glut_background_genes,
        union(peak_mcclay_genes_glut_bg, glut_up_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

peakMcClay_glut_down_fisher <- fisher.test(
  matrix(
    c(
      length(peakMcClay_glut_down_overlap),
      length(setdiff(glut_down_genes, peak_mcclay_genes_glut_bg)),
      length(setdiff(peak_mcclay_genes_glut_bg, glut_down_genes)),
      length(setdiff(
        glut_background_genes,
        union(peak_mcclay_genes_glut_bg, glut_down_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)


## Forrest peaks vs glut up/down genes
peakForrest_glut_up_overlap <- intersect(
  peak_forrest_genes_glut_bg,
  glut_up_genes
)

peakForrest_glut_down_overlap <- intersect(
  peak_forrest_genes_glut_bg,
  glut_down_genes
)

peakForrest_glut_up_fisher <- fisher.test(
  matrix(
    c(
      length(peakForrest_glut_up_overlap),
      length(setdiff(glut_up_genes, peak_forrest_genes_glut_bg)),
      length(setdiff(peak_forrest_genes_glut_bg, glut_up_genes)),
      length(setdiff(
        glut_background_genes,
        union(peak_forrest_genes_glut_bg, glut_up_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

peakForrest_glut_down_fisher <- fisher.test(
  matrix(
    c(
      length(peakForrest_glut_down_overlap),
      length(setdiff(glut_down_genes, peak_forrest_genes_glut_bg)),
      length(setdiff(peak_forrest_genes_glut_bg, glut_down_genes)),
      length(setdiff(
        glut_background_genes,
        union(peak_forrest_genes_glut_bg, glut_down_genes)
      ))
    ),
    nrow = 2,
    byrow = TRUE
  )
)

####################################################################
# Summary of Glut UP / DOWN Enrichment
####################################################################

tcf4_glut_up_down_summary <- data.frame(
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
    length(peakNPC_glut_up_overlap),
    length(peakNPC_glut_down_overlap),
    length(peakMcClay_glut_up_overlap),
    length(peakMcClay_glut_down_overlap),
    length(peakForrest_glut_up_overlap),
    length(peakForrest_glut_down_overlap)
  ),

  Odds_Ratio = c(
    as.numeric(peakNPC_glut_up_fisher$estimate),
    as.numeric(peakNPC_glut_down_fisher$estimate),
    as.numeric(peakMcClay_glut_up_fisher$estimate),
    as.numeric(peakMcClay_glut_down_fisher$estimate),
    as.numeric(peakForrest_glut_up_fisher$estimate),
    as.numeric(peakForrest_glut_down_fisher$estimate)
  ),

  P_Value = c(
    peakNPC_glut_up_fisher$p.value,
    peakNPC_glut_down_fisher$p.value,
    peakMcClay_glut_up_fisher$p.value,
    peakMcClay_glut_down_fisher$p.value,
    peakForrest_glut_up_fisher$p.value,
    peakForrest_glut_down_fisher$p.value
  )
)

tcf4_glut_up_down_summary


####################################################################
# Extract Helpful Overlap Gene Lists for Pathway Analysis
## Day 3 NPC: all significant overlap genes
pathway_d3_npc_peakNPC_all <- npc_overlap_genes
pathway_d3_npc_peakMcClay_all <- mcclay_overlap_genes
pathway_d3_npc_peakForrest_all <- forrest_overlap_genes

## Day 3 NPC: upregulated overlap genes
pathway_d3_npc_peakNPC_up <- npc_up_overlap_genes
pathway_d3_npc_peakMcClay_up <- mcclay_up_overlap_genes
pathway_d3_npc_peakForrest_up <- forrest_up_overlap_genes

## Day 3 NPC: downregulated overlap genes
pathway_d3_npc_peakNPC_down <- npc_down_overlap_genes
pathway_d3_npc_peakMcClay_down <- mcclay_down_overlap_genes
pathway_d3_npc_peakForrest_down <- forrest_down_overlap_genes


## Day 14 Glut: all significant overlap genes
pathway_d14_glut_peakNPC_all <- peakNPC_glut_sig_overlap
pathway_d14_glut_peakMcClay_all <- peakMcClay_glut_sig_overlap
pathway_d14_glut_peakForrest_all <- peakForrest_glut_sig_overlap

## Day 14 Glut: upregulated overlap genes
pathway_d14_glut_peakNPC_up <- peakNPC_glut_up_overlap
pathway_d14_glut_peakMcClay_up <- peakMcClay_glut_up_overlap
pathway_d14_glut_peakForrest_up <- peakForrest_glut_up_overlap

## Day 14 Glut: downregulated overlap genes
pathway_d14_glut_peakNPC_down <- peakNPC_glut_down_overlap
pathway_d14_glut_peakMcClay_down <- peakMcClay_glut_down_overlap
pathway_d14_glut_peakForrest_down <- peakForrest_glut_down_overlap


# Summary of Gene Lists for Pathway Analysis
pathway_gene_list_summary <- data.frame(
  Gene_List = c(
    "D3_NPC_peakNPC_all",
    "D3_NPC_peakMcClay_all",
    "D3_NPC_peakForrest_all",
    "D3_NPC_peakNPC_up",
    "D3_NPC_peakMcClay_up",
    "D3_NPC_peakForrest_up",
    "D3_NPC_peakNPC_down",
    "D3_NPC_peakMcClay_down",
    "D3_NPC_peakForrest_down",

    "D14_GLUT_peakNPC_all",
    "D14_GLUT_peakMcClay_all",
    "D14_GLUT_peakForrest_all",
    "D14_GLUT_peakNPC_up",
    "D14_GLUT_peakMcClay_up",
    "D14_GLUT_peakForrest_up",
    "D14_GLUT_peakNPC_down",
    "D14_GLUT_peakMcClay_down",
    "D14_GLUT_peakForrest_down"
  ),

  Gene_Count = c(
    length(pathway_d3_npc_peakNPC_all),
    length(pathway_d3_npc_peakMcClay_all),
    length(pathway_d3_npc_peakForrest_all),
    length(pathway_d3_npc_peakNPC_up),
    length(pathway_d3_npc_peakMcClay_up),
    length(pathway_d3_npc_peakForrest_up),
    length(pathway_d3_npc_peakNPC_down),
    length(pathway_d3_npc_peakMcClay_down),
    length(pathway_d3_npc_peakForrest_down),

    length(pathway_d14_glut_peakNPC_all),
    length(pathway_d14_glut_peakMcClay_all),
    length(pathway_d14_glut_peakForrest_all),
    length(pathway_d14_glut_peakNPC_up),
    length(pathway_d14_glut_peakMcClay_up),
    length(pathway_d14_glut_peakForrest_up),
    length(pathway_d14_glut_peakNPC_down),
    length(pathway_d14_glut_peakMcClay_down),
    length(pathway_d14_glut_peakForrest_down)
  )
)

pathway_gene_list_summary

####################################################################
# Save Gene Lists for Pathway Analysis
####################################################################
write.csv(pathway_gene_list_summary, "Decon_Doost/Pathway_gene_list_summary.csv", row.names = FALSE)

write.csv(data.frame(Gene = pathway_d3_npc_peakNPC_all), "Decon_Doost/D3_NPC_peakNPC_all_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d3_npc_peakMcClay_all), "Decon_Doost/D3_NPC_peakMcClay_all_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d3_npc_peakForrest_all), "Decon_Doost/D3_NPC_peakForrest_all_overlap_genes.csv", row.names = FALSE)

write.csv(data.frame(Gene = pathway_d3_npc_peakNPC_up), "Decon_Doost/D3_NPC_peakNPC_up_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d3_npc_peakMcClay_up), "Decon_Doost/D3_NPC_peakMcClay_up_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d3_npc_peakForrest_up), "Decon_Doost/D3_NPC_peakForrest_up_overlap_genes.csv", row.names = FALSE)

write.csv(data.frame(Gene = pathway_d3_npc_peakNPC_down), "Decon_Doost/D3_NPC_peakNPC_down_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d3_npc_peakMcClay_down), "Decon_Doost/D3_NPC_peakMcClay_down_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d3_npc_peakForrest_down), "Decon_Doost/D3_NPC_peakForrest_down_overlap_genes.csv", row.names = FALSE)

write.csv(data.frame(Gene = pathway_d14_glut_peakNPC_all), "Decon_Doost/D14_GLUT_peakNPC_all_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d14_glut_peakMcClay_all), "Decon_Doost/D14_GLUT_peakMcClay_all_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d14_glut_peakForrest_all), "Decon_Doost/D14_GLUT_peakForrest_all_overlap_genes.csv", row.names = FALSE)

write.csv(data.frame(Gene = pathway_d14_glut_peakNPC_up), "Decon_Doost/D14_GLUT_peakNPC_up_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d14_glut_peakMcClay_up), "Decon_Doost/D14_GLUT_peakMcClay_up_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d14_glut_peakForrest_up), "Decon_Doost/D14_GLUT_peakForrest_up_overlap_genes.csv", row.names = FALSE)

write.csv(data.frame(Gene = pathway_d14_glut_peakNPC_down), "Decon_Doost/D14_GLUT_peakNPC_down_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d14_glut_peakMcClay_down), "Decon_Doost/D14_GLUT_peakMcClay_down_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pathway_d14_glut_peakForrest_down), "Decon_Doost/D14_GLUT_peakForrest_down_overlap_genes.csv", row.names = FALSE)
