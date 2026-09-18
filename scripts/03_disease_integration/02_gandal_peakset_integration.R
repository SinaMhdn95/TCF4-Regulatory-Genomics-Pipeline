########################################################################
# Transcriptome-wide isoform-level dysregulation in ASD, SCZ, and BD
# https://pubmed.ncbi.nlm.nih.gov/30545856/
# This table comes from:
# NIHMS1007359-supplement-Table_S1.xlsx
########################################################################
# Load libraries
library(readxl)
library(dplyr)
library(org.Hs.eg.db)   # Human gene annotation and ID conversion
library(AnnotationDbi)  # Human gene annotation and ID conversion
library(clusterProfiler)  #Pathway enrichment analysis
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene


# Import dataset
npc_genes_hg38_df <- read.csv(
  "npc_genes_hg38.csv")

mcclay_genes_hg38_df <- read.csv(
  "mcclay_genes_hg38.csv")

forrest_genes_hg38_df <- read.csv(
  "forrest_genes_hg38.csv")

gandal_S1_DGE <-  read_xlsx("NIHMS1007359-supplement-Table_S1.xlsx",
                             sheet = "DGE")
gandal_S1_DTE <- read_xlsx("NIHMS1007359-supplement-Table_S1.xlsx",
                           sheet = "DTE")
gandal_S1_DTU <- read_xlsx("NIHMS1007359-supplement-Table_S1.xlsx",
                           sheet = "DTU")
gandal_S1_isoOnlyDE <- read_xlsx("NIHMS1007359-supplement-Table_S1.xlsx",
                           sheet = "isoformOnlyDE", skip = 1)


# Import hg19 peaks
npc_hg19_peaks <- import("NPC_ab21_idr_500bp_summit.bed")
mcclay_hg19_peaks <- import("McClay_TCF4_consensus_hg19.bed")
forrest_hg19_peaks <- import("Forrest_hg19_peaks.bed")

# dataset exploration
dim(gandal_S1_DGE)
dim(gandal_S1_DTE)
dim(gandal_S1_DTU)
dim(gandal_S1_isoOnlyDE)

colnames(gandal_S1_DGE)
colnames(gandal_S1_DTE)
colnames(gandal_S1_DTU)
colnames(gandal_S1_isoOnlyDE)

length(npc_hg19_peaks)
length(mcclay_hg19_peaks)
length(forrest_hg19_peaks)

# Create hg19 TCF4 peak-associated gene lists
# Annotate hg19 peaks
# tssRegion = c(-3000, 3000) means promoter is defined as
# 3 kb upstream and 3 kb downstream of the transcription start site.
npc_hg19_peak_anno <- annotatePeak(
  npc_hg19_peaks,
  TxDb = txdb_hg19,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

mcclay_hg19_peak_anno <- annotatePeak(
  mcclay_hg19_peaks,
  TxDb = txdb_hg19,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

forrest_hg19_peak_anno <- annotatePeak(
  forrest_hg19_peaks,
  TxDb = txdb_hg19,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

# Convert peak annotations to data frames
npc_hg19_peak_anno_df <- as.data.frame(
  npc_hg19_peak_anno
)

mcclay_hg19_peak_anno_df <- as.data.frame(
  mcclay_hg19_peak_anno
)

forrest_hg19_peak_anno_df <- as.data.frame(
  forrest_hg19_peak_anno
)

# Extract clean hg19 gene lists
# SYMBOL is the gene symbol column added by annoDb = "org.Hs.eg.db".
# We remove missing values and duplicated symbols
npc_genes_hg19 <- unique(
  na.omit(npc_hg19_peak_anno_df$SYMBOL)
)

mcclay_genes_hg19 <- unique(
  na.omit(mcclay_hg19_peak_anno_df$SYMBOL)
)

forrest_genes_hg19 <- unique(
  na.omit(forrest_hg19_peak_anno_df$SYMBOL)
)

# Check hg19 gene list sizes
tcf4_gene_list_hg19_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC_hg19",
    "McClay_hg19",
    "Forrest_hg19"
  ),

  Peak_Count = c(
    length(npc_hg19_peaks),
    length(mcclay_hg19_peaks),
    length(forrest_hg19_peaks)
  ),

  Gene_Count = c(
    length(npc_genes_hg19),
    length(mcclay_genes_hg19),
    length(forrest_genes_hg19)
  )
)

tcf4_gene_list_hg19_summary


# Save full hg19 peak annotation tables
write.csv(
  npc_hg19_peak_anno_df,
  "npc_hg19_peak_annotation.csv",
  row.names = FALSE
)

write.csv(
  mcclay_hg19_peak_anno_df,
  "mcclay_hg19_peak_annotation.csv",
  row.names = FALSE
)

write.csv(
  forrest_hg19_peak_anno_df,
  "forrest_hg19_peak_annotation.csv",
  row.names = FALSE
)

####################################################################
# Prepare Gandal Table S1 DGE Data
# Goal:
# 1. Clean the DGE table
# 2. Extract SCZ, ASD, and BD significant genes
# 3. Create SCZ, ASD, and BD upregulated and downregulated gene lists
# 4. Create the background gene set for Fisher tests
####################################################################
# Clean DGE table
gandal_S1_DGE_clean <- gandal_S1_DGE %>%
  filter(
    !is.na(gene_name)
  )

# Create background, SCZ, ASD, and BD datasets
gandal_SCZ_background_genes <- unique(
  gandal_S1_DGE_clean$gene_name
)


# SCZ #
# Significant SCZ DGE genes
# The paper uses FDR < 0.05 for DGE
gandal_S1_DGE_SCZ_sig_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$SCZ.fdr < 0.05
  ]
)

# SCZ upregulated genes
# log2FC > 0 means higher expression in SCZ compared with controls
gandal_S1_DGE_SCZ_up_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$SCZ.fdr < 0.05 &
      gandal_S1_DGE_clean$SCZ.log2FC > 0
  ]
)

# SCZ downregulated genes
# log2FC < 0 means lower expression in SCZ compared with controls
gandal_S1_DGE_SCZ_down_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$SCZ.fdr < 0.05 &
      gandal_S1_DGE_clean$SCZ.log2FC < 0
  ]
)

# Summarize SCZ DGE gene counts
gandal_SCZ_DGE_summary <- data.frame(
  Category = c(
    "Background",
    "SCZ_significant_DGE",
    "SCZ_upregulated",
    "SCZ_downregulated"
  ),
  Gene_Count = c(
    length(gandal_SCZ_background_genes),
    length(gandal_S1_DGE_SCZ_sig_genes),
    length(gandal_S1_DGE_SCZ_up_genes),
    length(gandal_S1_DGE_SCZ_down_genes)
  )
)

gandal_SCZ_DGE_summary


# ASD #
# Significant ASD DGE genes
# FDR < 0.05 means significantly differentially expressed in ASD.
gandal_S1_DGE_ASD_sig_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$ASD.fdr < 0.05
  ]
)

# ASD upregulated genes
# log2FC > 0 means higher expression in ASD compared with controls.
gandal_S1_DGE_ASD_up_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$ASD.fdr < 0.05 &
      gandal_S1_DGE_clean$ASD.log2FC > 0
  ]
)

# ASD downregulated genes
# log2FC < 0 means lower expression in ASD compared with controls.
gandal_S1_DGE_ASD_down_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$ASD.fdr < 0.05 &
      gandal_S1_DGE_clean$ASD.log2FC < 0
  ]
)


# Summarize ASD DGE gene counts
gandal_ASD_DGE_summary <- data.frame(

  Category = c(
    "Background",
    "ASD_significant_DGE",
    "ASD_upregulated",
    "ASD_downregulated"
  ),

  Gene_Count = c(
    length(gandal_SCZ_background_genes),
    length(gandal_S1_DGE_ASD_sig_genes),
    length(gandal_S1_DGE_ASD_up_genes),
    length(gandal_S1_DGE_ASD_down_genes)
  )
)

gandal_ASD_DGE_summary


# BD #
# Significant BD DGE genes
# FDR < 0.05 means significantly differentially expressed in BD.
gandal_S1_DGE_BD_sig_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$BD.fdr < 0.05
  ]
)

# BD upregulated genes
# log2FC > 0 means higher expression in BD compared with controls.
gandal_S1_DGE_BD_up_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$BD.fdr < 0.05 &
      gandal_S1_DGE_clean$BD.log2FC > 0
  ]
)

# BD downregulated genes
# log2FC < 0 means lower expression in BD compared with controls.
gandal_S1_DGE_BD_down_genes <- unique(
  gandal_S1_DGE_clean$gene_name[
    gandal_S1_DGE_clean$BD.fdr < 0.05 &
      gandal_S1_DGE_clean$BD.log2FC < 0
  ]
)


# Summarize BD DGE gene counts
gandal_BD_DGE_summary <- data.frame(

  Category = c(
    "Background",
    "BD_significant_DGE",
    "BD_upregulated",
    "BD_downregulated"
  ),

  Gene_Count = c(
    length(gandal_SCZ_background_genes),
    length(gandal_S1_DGE_BD_sig_genes),
    length(gandal_S1_DGE_BD_up_genes),
    length(gandal_S1_DGE_BD_down_genes)
  )
)

gandal_BD_DGE_summary


# Combined DGE Summary for ASD, SCZ, and BD
gandal_all_DGE_summary <- data.frame(

  Disorder = c(
    "ASD",
    "SCZ",
    "BD"
  ),

  Background_Genes = c(
    length(gandal_SCZ_background_genes),
    length(gandal_SCZ_background_genes),
    length(gandal_SCZ_background_genes)
  ),

  Significant_DGE = c(
    length(gandal_S1_DGE_ASD_sig_genes),
    length(gandal_S1_DGE_SCZ_sig_genes),
    length(gandal_S1_DGE_BD_sig_genes)
  ),

  Upregulated = c(
    length(gandal_S1_DGE_ASD_up_genes),
    length(gandal_S1_DGE_SCZ_up_genes),
    length(gandal_S1_DGE_BD_up_genes)
  ),

  Downregulated = c(
    length(gandal_S1_DGE_ASD_down_genes),
    length(gandal_S1_DGE_SCZ_down_genes),
    length(gandal_S1_DGE_BD_down_genes)
  )
)

gandal_all_DGE_summary

####################################################################
# Step 2: Import TCF4 Peak-Associated Gene Lists
# Goal:
# 1. Extract gene symbols
# 2. Remove missing values and duplicates
# 3. Prepare gene sets for overlap with Gandal DGE results
####################################################################

# Extract hg19 gene symbols
# These gene lists were created from hg19 peak annotation.
# unique() removes duplicated genes.
# na.omit() removes missing values.

tcf4_npc_genes <- unique(
  na.omit(npc_genes_hg19)
)

tcf4_mcclay_genes <- unique(
  na.omit(mcclay_genes_hg19)
)

tcf4_forrest_genes <- unique(
  na.omit(forrest_genes_hg19)
)


# Check TCF4 gene list sizes
tcf4_gene_list_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC_peaks_hg19",
    "McClay_peaks_hg19",
    "Forrest_peaks_hg19"
  ),

  Gene_Count = c(
    length(tcf4_npc_genes),
    length(tcf4_mcclay_genes),
    length(tcf4_forrest_genes)
  )
)

tcf4_gene_list_summary

####################################################################
# Step 3: Restrict TCF4 Gene Sets to Gandal DGE Background
# Goal:
# Use only genes that were actually tested in the Gandal DGE analysis
####################################################################
# The background is all genes in the Gandal DGE table with valid gene symbols
gandal_DGE_background_genes <- gandal_SCZ_background_genes


# Restrict each TCF4 gene set to the Gandal DGE background
tcf4_npc_genes_gandal_bg <- intersect(
  tcf4_npc_genes,
  gandal_DGE_background_genes
)

tcf4_mcclay_genes_gandal_bg <- intersect(
  tcf4_mcclay_genes,
  gandal_DGE_background_genes
)

tcf4_forrest_genes_gandal_bg <- intersect(
  tcf4_forrest_genes,
  gandal_DGE_background_genes
)


# Check how many TCF4 genes are present in Gandal background
tcf4_gandal_background_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC_peaks",
    "McClay_peaks",
    "Forrest_peaks"
  ),

  TCF4_Genes_in_Gandal_Background = c(
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg)
  )
)

tcf4_gandal_background_summary

####################################################################
# Step 4: SCZ DGE vs TCF4 Peak-Associated Gene Sets
# Goal:
# 1. Find overlap between SCZ DGE genes and each TCF4 gene set
# 2. Run Fisher exact test for enrichment
# 3. Create a clean summary table
####################################################################
# Raw overlaps: SCZ significant DGE genes
# NPC peaks vs SCZ DGE genes
gandal_SCZ_npc_overlap <- intersect(
  tcf4_npc_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

# McClay peaks vs SCZ DGE genes
gandal_SCZ_mcclay_overlap <- intersect(
  tcf4_mcclay_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

# Forrest peaks vs SCZ DGE genes
gandal_SCZ_forrest_overlap <- intersect(
  tcf4_forrest_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

# Check raw overlap sizes
gandal_SCZ_raw_overlap_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC_peaks",
    "McClay_peaks",
    "Forrest_peaks"
  ),

  SCZ_DGE_Overlap_Genes = c(
    length(gandal_SCZ_npc_overlap),
    length(gandal_SCZ_mcclay_overlap),
    length(gandal_SCZ_forrest_overlap)
  )
)

gandal_SCZ_raw_overlap_summary


# Fisher exact test: NPC peaks vs SCZ DGE genes
# A: genes that are both SCZ DGE and NPC TCF4 peak-associated
gandal_SCZ_npc_deg_and_peak <- intersect(
  tcf4_npc_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

# B: SCZ DGE genes that are not NPC TCF4 peak-associated
gandal_SCZ_npc_deg_only <- setdiff(
  gandal_S1_DGE_SCZ_sig_genes,
  tcf4_npc_genes_gandal_bg
)

# C: NPC TCF4 peak genes that are not SCZ DGE genes
gandal_SCZ_npc_peak_only <- setdiff(
  tcf4_npc_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

# D: background genes that are neither SCZ DGE nor NPC TCF4 peak-associated
gandal_SCZ_npc_neither <- setdiff(
  gandal_DGE_background_genes,
  union(
    tcf4_npc_genes_gandal_bg,
    gandal_S1_DGE_SCZ_sig_genes
  )
)

# Create 2 x 2 Fisher table
gandal_SCZ_npc_fisher_table <- matrix(
  c(
    length(gandal_SCZ_npc_deg_and_peak),
    length(gandal_SCZ_npc_deg_only),
    length(gandal_SCZ_npc_peak_only),
    length(gandal_SCZ_npc_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(gandal_SCZ_npc_fisher_table) <- c(
  "SCZ_DGE",
  "Not_SCZ_DGE"
)

colnames(gandal_SCZ_npc_fisher_table) <- c(
  "NPC_TCF4_peak_gene",
  "No_NPC_TCF4_peak"
)

gandal_SCZ_npc_fisher_table

gandal_SCZ_npc_fisher <- fisher.test(
  gandal_SCZ_npc_fisher_table
)

gandal_SCZ_npc_fisher


# Fisher exact test: McClay peaks vs SCZ DGE genes
gandal_SCZ_mcclay_deg_and_peak <- intersect(
  tcf4_mcclay_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

gandal_SCZ_mcclay_deg_only <- setdiff(
  gandal_S1_DGE_SCZ_sig_genes,
  tcf4_mcclay_genes_gandal_bg
)

gandal_SCZ_mcclay_peak_only <- setdiff(
  tcf4_mcclay_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

gandal_SCZ_mcclay_neither <- setdiff(
  gandal_DGE_background_genes,
  union(
    tcf4_mcclay_genes_gandal_bg,
    gandal_S1_DGE_SCZ_sig_genes
  )
)

gandal_SCZ_mcclay_fisher_table <- matrix(
  c(
    length(gandal_SCZ_mcclay_deg_and_peak),
    length(gandal_SCZ_mcclay_deg_only),
    length(gandal_SCZ_mcclay_peak_only),
    length(gandal_SCZ_mcclay_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(gandal_SCZ_mcclay_fisher_table) <- c(
  "SCZ_DGE",
  "Not_SCZ_DGE"
)

colnames(gandal_SCZ_mcclay_fisher_table) <- c(
  "McClay_TCF4_peak_gene",
  "No_McClay_TCF4_peak"
)

gandal_SCZ_mcclay_fisher_table

gandal_SCZ_mcclay_fisher <- fisher.test(
  gandal_SCZ_mcclay_fisher_table
)

gandal_SCZ_mcclay_fisher


# Fisher exact test: Forrest peaks vs SCZ DGE genes
gandal_SCZ_forrest_deg_and_peak <- intersect(
  tcf4_forrest_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

gandal_SCZ_forrest_deg_only <- setdiff(
  gandal_S1_DGE_SCZ_sig_genes,
  tcf4_forrest_genes_gandal_bg
)

gandal_SCZ_forrest_peak_only <- setdiff(
  tcf4_forrest_genes_gandal_bg,
  gandal_S1_DGE_SCZ_sig_genes
)

gandal_SCZ_forrest_neither <- setdiff(
  gandal_DGE_background_genes,
  union(
    tcf4_forrest_genes_gandal_bg,
    gandal_S1_DGE_SCZ_sig_genes
  )
)

gandal_SCZ_forrest_fisher_table <- matrix(
  c(
    length(gandal_SCZ_forrest_deg_and_peak),
    length(gandal_SCZ_forrest_deg_only),
    length(gandal_SCZ_forrest_peak_only),
    length(gandal_SCZ_forrest_neither)
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(gandal_SCZ_forrest_fisher_table) <- c(
  "SCZ_DGE",
  "Not_SCZ_DGE"
)

colnames(gandal_SCZ_forrest_fisher_table) <- c(
  "Forrest_TCF4_peak_gene",
  "No_Forrest_TCF4_peak"
)

gandal_SCZ_forrest_fisher_table

gandal_SCZ_forrest_fisher <- fisher.test(
  gandal_SCZ_forrest_fisher_table
)

gandal_SCZ_forrest_fisher


# Summary table: SCZ DGE enrichment in TCF4 peak-associated genes
gandal_SCZ_TCF4_fisher_summary <- data.frame(

  Disorder = c(
    "SCZ",
    "SCZ",
    "SCZ"
  ),

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Background_Genes = c(
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes)
  ),

  SCZ_DGE_Genes = c(
    length(gandal_S1_DGE_SCZ_sig_genes),
    length(gandal_S1_DGE_SCZ_sig_genes),
    length(gandal_S1_DGE_SCZ_sig_genes)
  ),

  TCF4_Genes_In_Background = c(
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg)
  ),

  Overlap_Genes = c(
    length(gandal_SCZ_npc_deg_and_peak),
    length(gandal_SCZ_mcclay_deg_and_peak),
    length(gandal_SCZ_forrest_deg_and_peak)
  ),

  DGE_Only_Genes = c(
    length(gandal_SCZ_npc_deg_only),
    length(gandal_SCZ_mcclay_deg_only),
    length(gandal_SCZ_forrest_deg_only)
  ),

  Peak_Only_Genes = c(
    length(gandal_SCZ_npc_peak_only),
    length(gandal_SCZ_mcclay_peak_only),
    length(gandal_SCZ_forrest_peak_only)
  ),

  Neither_Genes = c(
    length(gandal_SCZ_npc_neither),
    length(gandal_SCZ_mcclay_neither),
    length(gandal_SCZ_forrest_neither)
  ),

  Odds_Ratio = c(
    as.numeric(gandal_SCZ_npc_fisher$estimate),
    as.numeric(gandal_SCZ_mcclay_fisher$estimate),
    as.numeric(gandal_SCZ_forrest_fisher$estimate)
  ),

  P_Value = c(
    gandal_SCZ_npc_fisher$p.value,
    gandal_SCZ_mcclay_fisher$p.value,
    gandal_SCZ_forrest_fisher$p.value
  )
)

gandal_SCZ_TCF4_fisher_summary


####################################################################
# Directional SCZ DGE Enrichment in TCF4 Peak-Associated Genes
# Goal:
# 1. Test whether SCZ upregulated genes are enriched in TCF4 gene sets
# 2. Test whether SCZ downregulated genes are enriched in TCF4 gene sets
# 3. Compare NPC, Dr. McClay, and Forrest TCF4 datasets
####################################################################
# Function for Fisher exact test
# This function takes:
# deg_genes        = disease gene list, such as SCZ up or SCZ down genes
# peak_genes       = TCF4 peak-associated genes in Gandal background
# background_genes = all genes tested in the Gandal DGE analysis
#
# It returns:
# overlap counts, Fisher odds ratio, and p-value

run_fisher_overlap <- function(
    deg_genes,
    peak_genes,
    background_genes
) {

  # Keep only genes inside the background
  deg_genes <- intersect(
    deg_genes,
    background_genes
  )

  peak_genes <- intersect(
    peak_genes,
    background_genes
  )

  # A: genes that are both DEG and TCF4 peak-associated
  overlap_genes <- intersect(
    deg_genes,
    peak_genes
  )

  # B: DEG genes that are not TCF4 peak-associated
  deg_only_genes <- setdiff(
    deg_genes,
    peak_genes
  )

  # C: TCF4 peak-associated genes that are not DEG
  peak_only_genes <- setdiff(
    peak_genes,
    deg_genes
  )

  # D: background genes that are neither DEG nor TCF4 peak-associated
  neither_genes <- setdiff(
    background_genes,
    union(
      deg_genes,
      peak_genes
    )
  )

  # Make 2 x 2 Fisher table
  fisher_table <- matrix(
    c(
      length(overlap_genes),
      length(deg_only_genes),
      length(peak_only_genes),
      length(neither_genes)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(fisher_table) <- c(
    "DEG",
    "Not_DEG"
  )

  colnames(fisher_table) <- c(
    "TCF4_peak_gene",
    "No_TCF4_peak"
  )

  # Run Fisher exact test
  fisher_result <- fisher.test(
    fisher_table
  )

  # Return clean result
  return(
    list(
      overlap_genes = overlap_genes,
      deg_only_genes = deg_only_genes,
      peak_only_genes = peak_only_genes,
      neither_genes = neither_genes,
      fisher_table = fisher_table,
      fisher_result = fisher_result
    )
  )
}


# Run directional Fisher tests for SCZ upregulated genes
# NPC peaks vs SCZ upregulated genes
gandal_SCZ_up_npc_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_SCZ_up_genes,
  peak_genes = tcf4_npc_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

# McClay peaks vs SCZ upregulated genes
gandal_SCZ_up_mcclay_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_SCZ_up_genes,
  peak_genes = tcf4_mcclay_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

# Forrest peaks vs SCZ upregulated genes
gandal_SCZ_up_forrest_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_SCZ_up_genes,
  peak_genes = tcf4_forrest_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)



# Run directional Fisher tests for SCZ downregulated genes
# NPC peaks vs SCZ downregulated genes
gandal_SCZ_down_npc_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_SCZ_down_genes,
  peak_genes = tcf4_npc_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

# McClay peaks vs SCZ downregulated genes
gandal_SCZ_down_mcclay_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_SCZ_down_genes,
  peak_genes = tcf4_mcclay_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

# Forrest peaks vs SCZ downregulated genes
gandal_SCZ_down_forrest_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_SCZ_down_genes,
  peak_genes = tcf4_forrest_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)


# Create directional SCZ enrichment summary table
gandal_SCZ_directional_TCF4_summary <- data.frame(

  Disorder = c(
    "SCZ",
    "SCZ",
    "SCZ",
    "SCZ",
    "SCZ",
    "SCZ"
  ),

  Direction = c(
    "Up",
    "Up",
    "Up",
    "Down",
    "Down",
    "Down"
  ),

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest",
    "NPC",
    "McClay",
    "Forrest"
  ),

  Background_Genes = c(
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes)
  ),

  Directional_DGE_Genes = c(
    length(gandal_S1_DGE_SCZ_up_genes),
    length(gandal_S1_DGE_SCZ_up_genes),
    length(gandal_S1_DGE_SCZ_up_genes),
    length(gandal_S1_DGE_SCZ_down_genes),
    length(gandal_S1_DGE_SCZ_down_genes),
    length(gandal_S1_DGE_SCZ_down_genes)
  ),

  TCF4_Genes_In_Background = c(
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg),
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg)
  ),

  Overlap_Genes = c(
    length(gandal_SCZ_up_npc_fisher$overlap_genes),
    length(gandal_SCZ_up_mcclay_fisher$overlap_genes),
    length(gandal_SCZ_up_forrest_fisher$overlap_genes),
    length(gandal_SCZ_down_npc_fisher$overlap_genes),
    length(gandal_SCZ_down_mcclay_fisher$overlap_genes),
    length(gandal_SCZ_down_forrest_fisher$overlap_genes)
  ),

  DGE_Only_Genes = c(
    length(gandal_SCZ_up_npc_fisher$deg_only_genes),
    length(gandal_SCZ_up_mcclay_fisher$deg_only_genes),
    length(gandal_SCZ_up_forrest_fisher$deg_only_genes),
    length(gandal_SCZ_down_npc_fisher$deg_only_genes),
    length(gandal_SCZ_down_mcclay_fisher$deg_only_genes),
    length(gandal_SCZ_down_forrest_fisher$deg_only_genes)
  ),

  Peak_Only_Genes = c(
    length(gandal_SCZ_up_npc_fisher$peak_only_genes),
    length(gandal_SCZ_up_mcclay_fisher$peak_only_genes),
    length(gandal_SCZ_up_forrest_fisher$peak_only_genes),
    length(gandal_SCZ_down_npc_fisher$peak_only_genes),
    length(gandal_SCZ_down_mcclay_fisher$peak_only_genes),
    length(gandal_SCZ_down_forrest_fisher$peak_only_genes)
  ),

  Neither_Genes = c(
    length(gandal_SCZ_up_npc_fisher$neither_genes),
    length(gandal_SCZ_up_mcclay_fisher$neither_genes),
    length(gandal_SCZ_up_forrest_fisher$neither_genes),
    length(gandal_SCZ_down_npc_fisher$neither_genes),
    length(gandal_SCZ_down_mcclay_fisher$neither_genes),
    length(gandal_SCZ_down_forrest_fisher$neither_genes)
  ),

  Odds_Ratio = c(
    as.numeric(gandal_SCZ_up_npc_fisher$fisher_result$estimate),
    as.numeric(gandal_SCZ_up_mcclay_fisher$fisher_result$estimate),
    as.numeric(gandal_SCZ_up_forrest_fisher$fisher_result$estimate),
    as.numeric(gandal_SCZ_down_npc_fisher$fisher_result$estimate),
    as.numeric(gandal_SCZ_down_mcclay_fisher$fisher_result$estimate),
    as.numeric(gandal_SCZ_down_forrest_fisher$fisher_result$estimate)
  ),

  P_Value = c(
    gandal_SCZ_up_npc_fisher$fisher_result$p.value,
    gandal_SCZ_up_mcclay_fisher$fisher_result$p.value,
    gandal_SCZ_up_forrest_fisher$fisher_result$p.value,
    gandal_SCZ_down_npc_fisher$fisher_result$p.value,
    gandal_SCZ_down_mcclay_fisher$fisher_result$p.value,
    gandal_SCZ_down_forrest_fisher$fisher_result$p.value
  )
)

gandal_SCZ_directional_TCF4_summary

####################################################################
# ASD and BD DGE vs TCF4 Peak-Associated Gene Sets
# Goal:
# 1. Run overall enrichment for ASD and BD DGE genes
# 2. Run directional enrichment for ASD and BD up/down genes
# 3. Create one clean combined summary table
####################################################################

# Overall Fisher tests for ASD DGE genes
# NPC peaks vs ASD DGE genes
gandal_ASD_npc_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_sig_genes,
  peak_genes = tcf4_npc_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

# McClay peaks vs ASD DGE genes
gandal_ASD_mcclay_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_sig_genes,
  peak_genes = tcf4_mcclay_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

# Forrest peaks vs ASD DGE genes
gandal_ASD_forrest_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_sig_genes,
  peak_genes = tcf4_forrest_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)


# Overall Fisher tests for BD DGE genes
# NPC peaks vs BD DGE genes
gandal_BD_npc_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_sig_genes,
  peak_genes = tcf4_npc_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

# McClay peaks vs BD DGE genes
gandal_BD_mcclay_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_sig_genes,
  peak_genes = tcf4_mcclay_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

# Forrest peaks vs BD DGE genes
gandal_BD_forrest_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_sig_genes,
  peak_genes = tcf4_forrest_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)


# Overall ASD and BD enrichment summary table
gandal_ASD_BD_TCF4_fisher_summary <- data.frame(

  Disorder = c(
    "ASD",
    "ASD",
    "ASD",
    "BD",
    "BD",
    "BD"
  ),

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest",
    "NPC",
    "McClay",
    "Forrest"
  ),

  Background_Genes = c(
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes),
    length(gandal_DGE_background_genes)
  ),

  DGE_Genes = c(
    length(gandal_S1_DGE_ASD_sig_genes),
    length(gandal_S1_DGE_ASD_sig_genes),
    length(gandal_S1_DGE_ASD_sig_genes),
    length(gandal_S1_DGE_BD_sig_genes),
    length(gandal_S1_DGE_BD_sig_genes),
    length(gandal_S1_DGE_BD_sig_genes)
  ),

  TCF4_Genes_In_Background = c(
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg),
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg)
  ),

  Overlap_Genes = c(
    length(gandal_ASD_npc_fisher$overlap_genes),
    length(gandal_ASD_mcclay_fisher$overlap_genes),
    length(gandal_ASD_forrest_fisher$overlap_genes),
    length(gandal_BD_npc_fisher$overlap_genes),
    length(gandal_BD_mcclay_fisher$overlap_genes),
    length(gandal_BD_forrest_fisher$overlap_genes)
  ),

  DGE_Only_Genes = c(
    length(gandal_ASD_npc_fisher$deg_only_genes),
    length(gandal_ASD_mcclay_fisher$deg_only_genes),
    length(gandal_ASD_forrest_fisher$deg_only_genes),
    length(gandal_BD_npc_fisher$deg_only_genes),
    length(gandal_BD_mcclay_fisher$deg_only_genes),
    length(gandal_BD_forrest_fisher$deg_only_genes)
  ),

  Peak_Only_Genes = c(
    length(gandal_ASD_npc_fisher$peak_only_genes),
    length(gandal_ASD_mcclay_fisher$peak_only_genes),
    length(gandal_ASD_forrest_fisher$peak_only_genes),
    length(gandal_BD_npc_fisher$peak_only_genes),
    length(gandal_BD_mcclay_fisher$peak_only_genes),
    length(gandal_BD_forrest_fisher$peak_only_genes)
  ),

  Neither_Genes = c(
    length(gandal_ASD_npc_fisher$neither_genes),
    length(gandal_ASD_mcclay_fisher$neither_genes),
    length(gandal_ASD_forrest_fisher$neither_genes),
    length(gandal_BD_npc_fisher$neither_genes),
    length(gandal_BD_mcclay_fisher$neither_genes),
    length(gandal_BD_forrest_fisher$neither_genes)
  ),

  Odds_Ratio = c(
    as.numeric(gandal_ASD_npc_fisher$fisher_result$estimate),
    as.numeric(gandal_ASD_mcclay_fisher$fisher_result$estimate),
    as.numeric(gandal_ASD_forrest_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_npc_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_mcclay_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_forrest_fisher$fisher_result$estimate)
  ),

  P_Value = c(
    gandal_ASD_npc_fisher$fisher_result$p.value,
    gandal_ASD_mcclay_fisher$fisher_result$p.value,
    gandal_ASD_forrest_fisher$fisher_result$p.value,
    gandal_BD_npc_fisher$fisher_result$p.value,
    gandal_BD_mcclay_fisher$fisher_result$p.value,
    gandal_BD_forrest_fisher$fisher_result$p.value
  )
)

gandal_ASD_BD_TCF4_fisher_summary


# Directional Fisher tests for ASD upregulated genes
gandal_ASD_up_npc_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_up_genes,
  peak_genes = tcf4_npc_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

gandal_ASD_up_mcclay_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_up_genes,
  peak_genes = tcf4_mcclay_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

gandal_ASD_up_forrest_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_up_genes,
  peak_genes = tcf4_forrest_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)


# Directional Fisher tests for ASD downregulated genes
gandal_ASD_down_npc_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_down_genes,
  peak_genes = tcf4_npc_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

gandal_ASD_down_mcclay_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_down_genes,
  peak_genes = tcf4_mcclay_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

gandal_ASD_down_forrest_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_ASD_down_genes,
  peak_genes = tcf4_forrest_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)


# Directional Fisher tests for BD upregulated genes
gandal_BD_up_npc_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_up_genes,
  peak_genes = tcf4_npc_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

gandal_BD_up_mcclay_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_up_genes,
  peak_genes = tcf4_mcclay_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

gandal_BD_up_forrest_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_up_genes,
  peak_genes = tcf4_forrest_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)


# Directional Fisher tests for BD downregulated genes
gandal_BD_down_npc_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_down_genes,
  peak_genes = tcf4_npc_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

gandal_BD_down_mcclay_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_down_genes,
  peak_genes = tcf4_mcclay_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)

gandal_BD_down_forrest_fisher <- run_fisher_overlap(
  deg_genes = gandal_S1_DGE_BD_down_genes,
  peak_genes = tcf4_forrest_genes_gandal_bg,
  background_genes = gandal_DGE_background_genes
)


# Directional ASD and BD enrichment summary table
gandal_ASD_BD_directional_TCF4_summary <- data.frame(

  Disorder = c(
    "ASD", "ASD", "ASD",
    "ASD", "ASD", "ASD",
    "BD", "BD", "BD",
    "BD", "BD", "BD"
  ),

  Direction = c(
    "Up", "Up", "Up",
    "Down", "Down", "Down",
    "Up", "Up", "Up",
    "Down", "Down", "Down"
  ),

  TCF4_Dataset = c(
    "NPC", "McClay", "Forrest",
    "NPC", "McClay", "Forrest",
    "NPC", "McClay", "Forrest",
    "NPC", "McClay", "Forrest"
  ),

  Background_Genes = rep(
    length(gandal_DGE_background_genes),
    12
  ),

  Directional_DGE_Genes = c(
    length(gandal_S1_DGE_ASD_up_genes),
    length(gandal_S1_DGE_ASD_up_genes),
    length(gandal_S1_DGE_ASD_up_genes),
    length(gandal_S1_DGE_ASD_down_genes),
    length(gandal_S1_DGE_ASD_down_genes),
    length(gandal_S1_DGE_ASD_down_genes),
    length(gandal_S1_DGE_BD_up_genes),
    length(gandal_S1_DGE_BD_up_genes),
    length(gandal_S1_DGE_BD_up_genes),
    length(gandal_S1_DGE_BD_down_genes),
    length(gandal_S1_DGE_BD_down_genes),
    length(gandal_S1_DGE_BD_down_genes)
  ),

  TCF4_Genes_In_Background = c(
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg),
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg),
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg),
    length(tcf4_npc_genes_gandal_bg),
    length(tcf4_mcclay_genes_gandal_bg),
    length(tcf4_forrest_genes_gandal_bg)
  ),

  Overlap_Genes = c(
    length(gandal_ASD_up_npc_fisher$overlap_genes),
    length(gandal_ASD_up_mcclay_fisher$overlap_genes),
    length(gandal_ASD_up_forrest_fisher$overlap_genes),
    length(gandal_ASD_down_npc_fisher$overlap_genes),
    length(gandal_ASD_down_mcclay_fisher$overlap_genes),
    length(gandal_ASD_down_forrest_fisher$overlap_genes),
    length(gandal_BD_up_npc_fisher$overlap_genes),
    length(gandal_BD_up_mcclay_fisher$overlap_genes),
    length(gandal_BD_up_forrest_fisher$overlap_genes),
    length(gandal_BD_down_npc_fisher$overlap_genes),
    length(gandal_BD_down_mcclay_fisher$overlap_genes),
    length(gandal_BD_down_forrest_fisher$overlap_genes)
  ),

  Odds_Ratio = c(
    as.numeric(gandal_ASD_up_npc_fisher$fisher_result$estimate),
    as.numeric(gandal_ASD_up_mcclay_fisher$fisher_result$estimate),
    as.numeric(gandal_ASD_up_forrest_fisher$fisher_result$estimate),
    as.numeric(gandal_ASD_down_npc_fisher$fisher_result$estimate),
    as.numeric(gandal_ASD_down_mcclay_fisher$fisher_result$estimate),
    as.numeric(gandal_ASD_down_forrest_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_up_npc_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_up_mcclay_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_up_forrest_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_down_npc_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_down_mcclay_fisher$fisher_result$estimate),
    as.numeric(gandal_BD_down_forrest_fisher$fisher_result$estimate)
  ),

  P_Value = c(
    gandal_ASD_up_npc_fisher$fisher_result$p.value,
    gandal_ASD_up_mcclay_fisher$fisher_result$p.value,
    gandal_ASD_up_forrest_fisher$fisher_result$p.value,
    gandal_ASD_down_npc_fisher$fisher_result$p.value,
    gandal_ASD_down_mcclay_fisher$fisher_result$p.value,
    gandal_ASD_down_forrest_fisher$fisher_result$p.value,
    gandal_BD_up_npc_fisher$fisher_result$p.value,
    gandal_BD_up_mcclay_fisher$fisher_result$p.value,
    gandal_BD_up_forrest_fisher$fisher_result$p.value,
    gandal_BD_down_npc_fisher$fisher_result$p.value,
    gandal_BD_down_mcclay_fisher$fisher_result$p.value,
    gandal_BD_down_forrest_fisher$fisher_result$p.value
  )
)

gandal_ASD_BD_directional_TCF4_summary

# PERMUTATION for DEG
# Gene-set permutation analysis for Gandal DGE and TCF4-bound genes
# Prepare the common Gandal background
gandal_DGE_background_genes <- unique(
  na.omit(
    gandal_DGE_background_genes
  )
)

# Create a table describing all 27 comparisons
# We have:
# 3 disorders: ASD, SCZ and BD
# 3 directions: All, Up and Down
# 3 TCF4 datasets: NPC, McClay and Forrest
# Therefore:3 × 3 × 3 = 27 permutation tests
permutation_comparisons <- tibble(

  Disorder = rep(
    c(
      "ASD",
      "SCZ",
      "BD"
    ),
    each = 9
  ),

  Direction = rep(
    rep(
      c(
        "All",
        "Up",
        "Down"
      ),
      each = 3
    ),
    times = 3
  ),

  TCF4_Dataset = rep(
    c(
      "NPC",
      "McClay",
      "Forrest"
    ),
    times = 9
  ),

  # Each element of this column contains one Gandal DEG gene list
  DGE_Genes = c(

    # ASD gene lists
    rep(
      list(gandal_S1_DGE_ASD_sig_genes),
      times = 3
    ),

    rep(
      list(gandal_S1_DGE_ASD_up_genes),
      times = 3
    ),

    rep(
      list(gandal_S1_DGE_ASD_down_genes),
      times = 3
    ),

    # SCZ gene lists
    rep(
      list(gandal_S1_DGE_SCZ_sig_genes),
      times = 3
    ),

    rep(
      list(gandal_S1_DGE_SCZ_up_genes),
      times = 3
    ),

    rep(
      list(gandal_S1_DGE_SCZ_down_genes),
      times = 3
    ),

    # BD gene lists
    rep(
      list(gandal_S1_DGE_BD_sig_genes),
      times = 3
    ),

    rep(
      list(gandal_S1_DGE_BD_up_genes),
      times = 3
    ),

    rep(
      list(gandal_S1_DGE_BD_down_genes),
      times = 3
    )
  ),

  # The NPC, McClay and Forrest TCF4 gene lists are repeated
  # once for every disorder and direction.
  TCF4_Genes = rep(
    list(
      tcf4_npc_genes_gandal_bg,
      tcf4_mcclay_genes_gandal_bg,
      tcf4_forrest_genes_gandal_bg
    ),
    times = 9
  )
)


# Inspect the comparison table
permutation_comparisons %>%
  select(
    Disorder,
    Direction,
    TCF4_Dataset
  )


# Set the number of permutations: A value of 10,000 provides a minimum possible empirical P value of approximately 0.0001.
number_of_permutations <- 10000

# Set the random seed
# set.seed() makes the random results reproducible.
# Running the code again with the same seed should give the same results.
set.seed(12345)


# Create empty result columns
number_of_comparisons <- nrow(
  permutation_comparisons
)

background_gene_count <- integer(
  number_of_comparisons
)

DGE_gene_count <- integer(
  number_of_comparisons
)

TCF4_gene_count <- integer(
  number_of_comparisons
)

observed_overlap <- integer(
  number_of_comparisons
)

mean_random_overlap <- numeric(
  number_of_comparisons
)

median_random_overlap <- numeric(
  number_of_comparisons
)

random_overlap_SD <- numeric(
  number_of_comparisons
)

random_overlap_95_lower <- numeric(
  number_of_comparisons
)

random_overlap_95_upper <- numeric(
  number_of_comparisons
)

fold_enrichment <- numeric(
  number_of_comparisons
)

empirical_P_value <- numeric(
  number_of_comparisons
)

fisher_odds_ratio <- numeric(
  number_of_comparisons
)

fisher_P_value <- numeric(
  number_of_comparisons
)


# Run the 27 permutation analyses
# This is the only loop used in the analysis.
# Each repetition processes one disorder, direction and TCF4 dataset.
for (
  comparison_number in seq_len(number_of_comparisons)
) {

  # Extract the current Gandal DEG gene list

  # [[ ]] extracts the complete gene vector stored in one table cell.
  current_DGE_genes <- permutation_comparisons$DGE_Genes[[comparison_number]]


  # Extract the current TCF4-associated gene list
  current_TCF4_genes <- permutation_comparisons$TCF4_Genes[[comparison_number]]


  # Restrict the DEG genes to the common Gandal background
  # intersect() keeps genes found in both objects.
  current_DGE_genes <- intersect(
    unique(
      na.omit(
        current_DGE_genes
      )
    ),
    gandal_DGE_background_genes
  )


  # Restrict the TCF4 genes to the same background
  current_TCF4_genes <- intersect(
    unique(
      na.omit(
        current_TCF4_genes
      )
    ),
    gandal_DGE_background_genes
  )


  # Count the real overlap
  current_observed_overlap <- length(
    intersect(
      current_DGE_genes,
      current_TCF4_genes
    )
  )


  # Generate the random overlap distribution
  # rhyper() generates the overlap expected when a random TCF4 gene
  # set of the same size is selected without replacement.
  # nn = number of random permutations
  # m = number of DEG genes in the background
  # n = number of non-DEG genes in the background
  # k = number of randomly selected genes, equal to the size of the
  #     real TCF4 gene set
  # This is mathematically equivalent to repeatedly selecting a
  # random TCF4 gene set and counting its overlap with the DEG list.
  current_random_overlaps <- rhyper(
    nn = number_of_permutations,
    m = length(current_DGE_genes),
    n = length(gandal_DGE_background_genes) -
      length(current_DGE_genes),
    k = length(current_TCF4_genes)
  )


  # Calculate the empirical one-sided enrichment P value
  # The question is:
  # How often was the random overlap at least as large as the
  # observed overlap?
  # The added 1 prevents the empirical P value from being exactly zero.
  current_empirical_P_value <- (
    sum(
      current_random_overlaps >= current_observed_overlap
    ) + 1
  ) / (
    number_of_permutations + 1
  )


  # Construct the corresponding Fisher table
  #                         DEG         Not DEG
  # TCF4-associated          A              B
  # Not TCF4-associated      C              D

  current_fisher_table <- matrix(
    c(
      current_observed_overlap,

      length(current_TCF4_genes) -
        current_observed_overlap,

      length(current_DGE_genes) -
        current_observed_overlap,

      length(gandal_DGE_background_genes) -
        length(
          union(
            current_DGE_genes,
            current_TCF4_genes
          )
        )
    ),

    nrow = 2,
    byrow = TRUE
  )


  rownames(
    current_fisher_table
  ) <- c(
    "TCF4_Associated",
    "Not_TCF4_Associated"
  )


  colnames(
    current_fisher_table
  ) <- c(
    "DGE",
    "Not_DGE"
  )


  # Run a one-sided Fisher enrichment test
  # alternative = "greater" asks whether TCF4-associated genes
  # contain more DEGs than expected.
  current_fisher_result <- fisher.test(
    current_fisher_table,
    alternative = "greater"
  )


  # Store the basic gene counts
  background_gene_count[
    comparison_number
  ] <- length(
    gandal_DGE_background_genes
  )

  DGE_gene_count[
    comparison_number
  ] <- length(
    current_DGE_genes
  )

  TCF4_gene_count[
    comparison_number
  ] <- length(
    current_TCF4_genes
  )


  # Store the observed overlap
  observed_overlap[
    comparison_number
  ] <- current_observed_overlap


  # Store summaries of the random overlap distribution
  mean_random_overlap[
    comparison_number
  ] <- mean(
    current_random_overlaps
  )

  median_random_overlap[
    comparison_number
  ] <- median(
    current_random_overlaps
  )

  random_overlap_SD[
    comparison_number
  ] <- sd(
    current_random_overlaps
  )

  random_overlap_95_lower[
    comparison_number
  ] <- quantile(
    current_random_overlaps,
    probabilities = 0.025,
    names = FALSE
  )

  random_overlap_95_upper[
    comparison_number
  ] <- quantile(
    current_random_overlaps,
    probabilities = 0.975,
    names = FALSE
  )


  # Calculate fold enrichment
  # A value greater than 1 means the observed overlap was larger
  # than the mean random overlap.
  fold_enrichment[
    comparison_number
  ] <- current_observed_overlap /
    mean(current_random_overlaps)


  # Store the permutation P value
  empirical_P_value[
    comparison_number
  ] <- current_empirical_P_value


  # Store the corresponding Fisher results
  fisher_odds_ratio[
    comparison_number
  ] <- as.numeric(
    current_fisher_result$estimate
  )

  fisher_P_value[
    comparison_number
  ] <- current_fisher_result$p.value
}


# Create the complete permutation result table
gene_permutation_results <- data.frame(

  Disorder = permutation_comparisons$Disorder,

  Direction = permutation_comparisons$Direction,

  TCF4_Dataset = permutation_comparisons$TCF4_Dataset,

  Background_Genes = background_gene_count,

  DGE_Genes = DGE_gene_count,

  TCF4_Genes = TCF4_gene_count,

  Observed_Overlap = observed_overlap,

  Mean_Random_Overlap = mean_random_overlap,

  Median_Random_Overlap = median_random_overlap,

  Random_Overlap_SD = random_overlap_SD,

  Random_95_Percent_Lower = random_overlap_95_lower,

  Random_95_Percent_Upper = random_overlap_95_upper,

  Fold_Enrichment = fold_enrichment,

  Fisher_Odds_Ratio = fisher_odds_ratio,

  Fisher_P_Value = fisher_P_value,

  Empirical_P_Value = empirical_P_value
)


# Correct for the 27 comparisons
# p.adjust() applies the Benjamini-Hochberg multiple-testing correction.
gene_permutation_results <- gene_permutation_results %>%
  mutate(

    Fisher_FDR = p.adjust(
      Fisher_P_Value,
      method = "BH"
    ),

    Empirical_FDR = p.adjust(
      Empirical_P_Value,
      method = "BH"
    )
  )


# Arrange results from strongest to weakest empirical evidence
gene_permutation_results <- gene_permutation_results %>%
  arrange(
    Empirical_FDR,
    Empirical_P_Value
  )


# Display the complete result table
gene_permutation_results









########################################################################
# TCF4 TRANSCRIPT / ISOFORM ANALYSIS USING THE GANDAL DATASET
#
# steps:
#   1. Reads Gandal differential transcript expression (DTE) results.
#   2. Reads Gandal differential transcript usage (DTU) results.
#   3. Creates a promoter around each transcript start site (TSS).
#   4. Tests whether TCF4 peaks overlap transcript promoters.
#   5. Combines transcripts that share the same TSS.
#   6. Tests DTE and DTU enrichment.
#   7. Performs within-gene permutation analysis.
#   8. Creates candidate tables in the R environment.
#
# This code DOES NOT save, delete, move, or modify any files.
#
# Important:
#   DTE = change in the amount of an individual transcript.
#   DTU = change in the relative usage of a transcript within its gene.
#
# Genome build:
#   Gandal transcript coordinates = hg19 / GENCODE v19
#   TCF4 peak files used below     = hg19
########################################################################


# Import files. Set TCF4_DATA_DIR to the local authorized data directory.
data_root <- Sys.getenv("TCF4_DATA_DIR", unset = "data/raw")
#paste0() joins pieces of text together without adding spaces.
#also used to construct a long file path
gandal_file <- file.path(
  data_root,
  "expression",
  "gandal_table_s1.xlsx"
)

npc_peak_file <- file.path(
  data_root,
  "peaks",
  "NPC_ab21_idr_500bp_summit.bed"
)

mcclay_peak_file <- file.path(
  data_root,
  "peaks",
  "McClay_TCF4_consensus_hg19.bed"
)

forrest_peak_file <- file.path(
  data_root,
  "peaks",
  "Forrest_hg19_peaks.bed"
)


# Read Gandal's DTE and DTU results
gandal_DTE <- read_xlsx(
  gandal_file,
  sheet = "DTE",
  na = c("", "NA")
)

gandal_DTU <- read_xlsx(
  gandal_file,
  sheet = "DTU",
  na = c("", "NA")
)


# Check dataset
dim(gandal_DTE)
dim(gandal_DTU)

colnames(gandal_DTE)
colnames(gandal_DTU)

head(gandal_DTE)
head(gandal_DTU)



# 5. CHECK THE NUMBER OF SIGNIFICANT DTE TRANSCRIPTS
# FDR < 0.05 defines significant differential transcript expression.
DTE_count_summary <- data.frame(
  Disorder = c(
    "ASD",
    "SCZ",
    "BD"
  ),

  Significant_DTE_Transcripts = c(
    sum(
      !is.na(gandal_DTE$ASD.fdr) &
        gandal_DTE$ASD.fdr < 0.05
    ),

    sum(
      !is.na(gandal_DTE$SCZ.fdr) &
        gandal_DTE$SCZ.fdr < 0.05
    ),

    sum(
      !is.na(gandal_DTE$BD.fdr) &
        gandal_DTE$BD.fdr < 0.05
    )
  ),

  Significant_DTE_Genes = c(
    length(
      unique(
        gandal_DTE$ensembl_gene_id[
          !is.na(gandal_DTE$ASD.fdr) &
            gandal_DTE$ASD.fdr < 0.05
        ]
      )
    ),

    length(
      unique(
        gandal_DTE$ensembl_gene_id[
          !is.na(gandal_DTE$SCZ.fdr) &
            gandal_DTE$SCZ.fdr < 0.05
        ]
      )
    ),

    length(
      unique(
        gandal_DTE$ensembl_gene_id[
          !is.na(gandal_DTE$BD.fdr) &
            gandal_DTE$BD.fdr < 0.05
        ]
      )
    )
  )
)

DTE_count_summary


# CHECK THE NUMBER OF SIGNIFICANT DTU TRANSCRIPTS
# FDR < 0.05 defines significant differential transcript usage.
# DTU.<disorder>.Value is the DTU effect size.
# Positive value = increased relative transcript usage.
# Negative value = decreased relative transcript usage.
DTU_count_summary <- data.frame(
  Disorder = c(
    "ASD",
    "SCZ",
    "BD"
  ),

  Significant_DTU_Transcripts = c(
    sum(
      !is.na(gandal_DTU$DTU.ASD.FDR) &
        gandal_DTU$DTU.ASD.FDR < 0.05
    ),

    sum(
      !is.na(gandal_DTU$DTU.SCZ.FDR) &
        gandal_DTU$DTU.SCZ.FDR < 0.05
    ),

    sum(
      !is.na(gandal_DTU$DTU.BD.FDR) &
        gandal_DTU$DTU.BD.FDR < 0.05
    )
  ),

  Significant_DTU_Genes = c(
    length(
      unique(
        gandal_DTU$ensembl_gene_id[
          !is.na(gandal_DTU$DTU.ASD.FDR) &
            gandal_DTU$DTU.ASD.FDR < 0.05
        ]
      )
    ),

    length(
      unique(
        gandal_DTU$ensembl_gene_id[
          !is.na(gandal_DTU$DTU.SCZ.FDR) &
            gandal_DTU$DTU.SCZ.FDR < 0.05
        ]
      )
    ),

    length(
      unique(
        gandal_DTU$ensembl_gene_id[
          !is.na(gandal_DTU$DTU.BD.FDR) &
            gandal_DTU$DTU.BD.FDR < 0.05
        ]
      )
    )
  )
)

DTU_count_summary



# CLEAN THE DTE TABLE
# Use the Ensembl transcript ID as the unique transcript identifier.
# Required columns:
#   ensembl_transcript_id
#   ensembl_gene_id
#   chromosome_name
#   transcript_start
#   transcript_end
#   strand

gandal_DTE_clean <- gandal_DTE %>%
  filter(
    !is.na(ensembl_transcript_id),
    !is.na(ensembl_gene_id),
    !is.na(chromosome_name),
    !is.na(transcript_start),
    !is.na(transcript_end),
    !is.na(strand)
  ) %>%

  distinct(
    ensembl_transcript_id,
    .keep_all = TRUE
  ) %>%

  mutate(
    chromosome = as.character(chromosome_name),

    chromosome = if_else(
      grepl("^chr", chromosome),
      chromosome,
      paste0("chr", chromosome)
    ),

    strand_character = case_when(
      strand == 1  ~ "+",
      strand == -1 ~ "-",
      TRUE         ~ "*"
    ),

    transcript_TSS = case_when(
      strand == 1  ~ as.integer(transcript_start),
      strand == -1 ~ as.integer(transcript_end),
      TRUE         ~ as.integer(transcript_start)
    )
  )


# KEEP STANDARD HUMAN CHROMOSOMES
standard_human_chromosomes <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY"
)

gandal_DTE_clean <- gandal_DTE_clean %>%
  filter(
    chromosome %in% standard_human_chromosomes
  )


# CLEAN THE DTU TABLE
gandal_DTU_clean <- gandal_DTU %>%
  filter(
    !is.na(ensembl_transcript_id)
  ) %>%

  distinct(
    ensembl_transcript_id,
    .keep_all = TRUE
  )


# JOIN DTE AND DTU RESULTS
# The resulting table contains:
#   - transcript coordinates
#   - DTE effect sizes and FDR values
#   - DTU effect sizes and FDR values
transcript_data <- gandal_DTE_clean %>%
  select(
    ensembl_transcript_id,
    ensembl_gene_id,
    external_gene_id,
    external_transcript_id,
    chromosome,
    transcript_start,
    transcript_end,
    strand,
    strand_character,
    transcript_TSS,
    transcript_length,

    ASD.log2FC,
    ASD.fdr,

    SCZ.log2FC,
    SCZ.fdr,

    BD.log2FC,
    BD.fdr
  ) %>%

  left_join(
    gandal_DTU_clean %>%
      select(
        ensembl_transcript_id,

        DTU.ASD.Value,
        DTU.ASD.FDR,

        DTU.SCZ.Value,
        DTU.SCZ.FDR,

        DTU.BD.Value,
        DTU.BD.FDR
      ),

    by = "ensembl_transcript_id"
  )

View(transcript_data)


# CREATE CLEAR SIGNIFICANCE LABELS
# "All"  = all significant transcripts
# "Up"   = significant transcripts with positive effect
# "Down" = significant transcripts with negative effect
transcript_data <- transcript_data %>%
  mutate(

    # ASD DTE labels
    DTE_ASD_All = (
      !is.na(ASD.fdr) &
        ASD.fdr < 0.05
    ),

    DTE_ASD_Up = (
      !is.na(ASD.fdr) &
        ASD.fdr < 0.05 &
        ASD.log2FC > 0
    ),

    DTE_ASD_Down = (
      !is.na(ASD.fdr) &
        ASD.fdr < 0.05 &
        ASD.log2FC < 0
    ),


    # SCZ DTE labels
    DTE_SCZ_All = (
      !is.na(SCZ.fdr) &
        SCZ.fdr < 0.05
    ),

    DTE_SCZ_Up = (
      !is.na(SCZ.fdr) &
        SCZ.fdr < 0.05 &
        SCZ.log2FC > 0
    ),

    DTE_SCZ_Down = (
      !is.na(SCZ.fdr) &
        SCZ.fdr < 0.05 &
        SCZ.log2FC < 0
    ),

    # BD DTE labels
    DTE_BD_All = (
      !is.na(BD.fdr) &
        BD.fdr < 0.05
    ),

    DTE_BD_Up = (
      !is.na(BD.fdr) &
        BD.fdr < 0.05 &
        BD.log2FC > 0
    ),

    DTE_BD_Down = (
      !is.na(BD.fdr) &
        BD.fdr < 0.05 &
        BD.log2FC < 0
    ),

    # ASD DTU labels
    DTU_ASD_All = (
      !is.na(DTU.ASD.FDR) &
        DTU.ASD.FDR < 0.05
    ),

    DTU_ASD_Up = (
      !is.na(DTU.ASD.FDR) &
        DTU.ASD.FDR < 0.05 &
        DTU.ASD.Value > 0
    ),

    DTU_ASD_Down = (
      !is.na(DTU.ASD.FDR) &
        DTU.ASD.FDR < 0.05 &
        DTU.ASD.Value < 0
    ),


    # SCZ DTU labels
    DTU_SCZ_All = (
      !is.na(DTU.SCZ.FDR) &
        DTU.SCZ.FDR < 0.05
    ),

    DTU_SCZ_Up = (
      !is.na(DTU.SCZ.FDR) &
        DTU.SCZ.FDR < 0.05 &
        DTU.SCZ.Value > 0
    ),

    DTU_SCZ_Down = (
      !is.na(DTU.SCZ.FDR) &
        DTU.SCZ.FDR < 0.05 &
        DTU.SCZ.Value < 0
    ),


    # BD DTU labels
    DTU_BD_All = (
      !is.na(DTU.BD.FDR) &
        DTU.BD.FDR < 0.05
    ),

    DTU_BD_Up = (
      !is.na(DTU.BD.FDR) &
        DTU.BD.FDR < 0.05 &
        DTU.BD.Value > 0
    ),

    DTU_BD_Down = (
      !is.na(DTU.BD.FDR) &
        DTU.BD.FDR < 0.05 &
        DTU.BD.Value < 0
    )
  )



# CONVERT TRANSCRIPTS TO GENOMIC RANGES
transcript_ranges <- GRanges(
  seqnames = transcript_data$chromosome,

  ranges = IRanges(
    start = as.integer(transcript_data$transcript_start),
    end = as.integer(transcript_data$transcript_end)
  ),
  strand = transcript_data$strand_character   # The strand determines where transcription begins
                                              # For the positive strand: TSS = transcript_start
                                              # For the negative strand: TSS = transcript_end
)



# CREATE A PROMOTER AROUND EACH TRANSCRIPT TSS
# Here, the transcript promoter is defined as approximately:
# 1 kb ← TSS → 1 kb
transcript_promoters_1kb <- promoters(
  transcript_ranges,
  upstream = 1000,
  downstream = 1000
)

# IMPORT TCF4 PEAK FILES
npc_peaks_hg19 <- import(
  npc_peak_file
)

mcclay_peaks_hg19 <- import(
  mcclay_peak_file
)

forrest_peaks_hg19 <- import(
  forrest_peak_file
)


# CORRECT CHROMOSOME X IN THE MCCLAY PEAK FILE
# The McClay peak file uses chr23 for chromosome X
# Change chr23 to chrX in the R object.
# We need to do this because R considers chr23 and chrX different chromosomes. Therefore, it will not recognize the overlap.
if ("chr23" %in% seqlevels(mcclay_peaks_hg19)) {  # seqlevels() returns the chromosome names present in the peak object
                                                  # return "chr1" "chr2" "chr3" ... "chr22" "chr23"

  mcclay_peaks_hg19 <- renameSeqlevels(
    mcclay_peaks_hg19,
    c("chr23" = "chrX")
  )
}


# REMOVE RANDOM AND NONSTANDARD CHROMOSOMES
# keepStandardChromosomes(): this is a function from the GenomeInfoDb package.
# It examines the chromosome names in a genomic object and retains the standard chromosomes.
# It will remove for example this chromosome names: chr1_random chrUn_gl000220
npc_peaks_hg19 <- keepStandardChromosomes(
  npc_peaks_hg19,
  species = "Homo_sapiens",   # Use the definition of standard chromosomes for humans
  pruning.mode = "coarse"     # When removing a nonstandard chromosome, also remove all peaks located on that chromosome.
)

mcclay_peaks_hg19 <- keepStandardChromosomes(
  mcclay_peaks_hg19,
  species = "Homo_sapiens",
  pruning.mode = "coarse"
)

forrest_peaks_hg19 <- keepStandardChromosomes(
  forrest_peaks_hg19,
  species = "Homo_sapiens",
  pruning.mode = "coarse"
)



# CHECK THE NUMBER OF PEAKS
TCF4_peak_count_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Peak_Count = c(
    length(npc_peaks_hg19),
    length(mcclay_peaks_hg19),
    length(forrest_peaks_hg19)
  )
)

TCF4_peak_count_summary


# DETERMINE WHICH TRANSCRIPT PROMOTERS OVERLAP TCF4 PEAKS
# TRUE  = transcript promoter overlaps at least one TCF4 peak
# FALSE = no overlap
transcript_data$Promoter_Bound_NPC <- overlapsAny(       # overlapsAny(): compares two sets of genomic regions
  transcript_promoters_1kb,  # contains the ±1 kb promoter region for every transcript
  npc_peaks_hg19,            # contains the genomic positions of the NPC TCF4 peaks.
  ignore.strand = TRUE       # means that R ignores the + or − strand when determining overlap.
)

transcript_data$Promoter_Bound_McClay <- overlapsAny(
  transcript_promoters_1kb,
  mcclay_peaks_hg19,
  ignore.strand = TRUE
)

transcript_data$Promoter_Bound_Forrest <- overlapsAny(
  transcript_promoters_1kb,
  forrest_peaks_hg19,
  ignore.strand = TRUE
)



# SUMMARIZE TRANSCRIPT PROMOTER BINDING
transcript_promoter_binding_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_Tested_Transcripts = c(
    nrow(transcript_data),
    nrow(transcript_data),
    nrow(transcript_data)
  ),

  Promoter_Bound_Transcripts = c(
    sum(transcript_data$Promoter_Bound_NPC),
    sum(transcript_data$Promoter_Bound_McClay),
    sum(transcript_data$Promoter_Bound_Forrest)
  ),

  Promoter_Bound_Genes = c(
    length(
      unique(
        transcript_data$ensembl_gene_id[
          transcript_data$Promoter_Bound_NPC
        ]
      )
    ),

    length(
      unique(
        transcript_data$ensembl_gene_id[
          transcript_data$Promoter_Bound_McClay
        ]
      )
    ),

    length(
      unique(
        transcript_data$ensembl_gene_id[
          transcript_data$Promoter_Bound_Forrest
        ]
      )
    )
  )
)

transcript_promoter_binding_summary



# CREATE A TSS CLUSTER IDENTIFIER
# Some transcripts share the same TSS.
# If two transcripts share the same gene, chromosome, strand and TSS,
# one TCF4 peak cannot distinguish between those transcripts.
# Therefore, transcripts sharing an identical TSS are combined.
# this code Give transcripts the same label when they belong to the same gene and begin at exactly the same genomic position.
transcript_data <- transcript_data %>%
  mutate(
    TSS_Cluster_ID = paste(   # paste() combines the four pieces of information into one text label
                              # two transcripts receive the same cluster ID only when all four values match
      ensembl_gene_id,
      chromosome,
      strand,
      transcript_TSS,
      sep = ":"
    )
  )



# 21. COLLAPSE TRANSCRIPTS THAT SHARE THE SAME TSS
tss_cluster_data <- transcript_data %>%
  group_by(                              # group_by() places transcripts into the same group when all these values match.
    TSS_Cluster_ID,
    ensembl_gene_id,
    external_gene_id,
    chromosome,
    strand,
    transcript_TSS
  ) %>%

  summarise(                             # takes all the transcript rows in each group and produces one summary row.

    # Transcript information
    Number_of_Transcripts = n(),         # n() counts the number of transcript rows in each TSS cluster.

    Ensembl_Transcript_IDs = paste(      # paste(..., collapse = ";") This combines multiple transcript IDs into one text value,
                                         # separated by semicolons. like ENST00000354452 and ENST00000356073 convering to ENST00000354452;ENST00000356073
      sort(
        unique(
          ensembl_transcript_id
        )
      ),
      collapse = ";"
    ),

    Transcript_Names = paste(            # for example TCF4-002, TCF4-201, and NA are converted to TCF4-002;TCF4-201
      sort(
        unique(
          external_transcript_id[
            !is.na(external_transcript_id)
          ]
        )
      ),
      collapse = ";"
    ),

    # TCF4 promoter binding
    Promoter_Bound_NPC = any(            # any() asks: Is at least one transcript in this TSS cluster labeled TRUE?
      Promoter_Bound_NPC
    ),

    Promoter_Bound_McClay = any(
      Promoter_Bound_McClay
    ),

    Promoter_Bound_Forrest = any(
      Promoter_Bound_Forrest
    ),


    # DTE significance
    DTE_ASD_All = any(DTE_ASD_All),      # This asks: Does at least one transcript in this TSS cluster have significant DTE in schizophrenia?
    DTE_ASD_Up = any(DTE_ASD_Up),        # These ask: Does at least one transcript increase?
    DTE_ASD_Down = any(DTE_ASD_Down),

    DTE_SCZ_All = any(DTE_SCZ_All),
    DTE_SCZ_Up = any(DTE_SCZ_Up),
    DTE_SCZ_Down = any(DTE_SCZ_Down),

    DTE_BD_All = any(DTE_BD_All),
    DTE_BD_Up = any(DTE_BD_Up),
    DTE_BD_Down = any(DTE_BD_Down),


    # DTU significance
    DTU_ASD_All = any(DTU_ASD_All),
    DTU_ASD_Up = any(DTU_ASD_Up),
    DTU_ASD_Down = any(DTU_ASD_Down),

    DTU_SCZ_All = any(DTU_SCZ_All),
    DTU_SCZ_Up = any(DTU_SCZ_Up),
    DTU_SCZ_Down = any(DTU_SCZ_Down),

    DTU_BD_All = any(DTU_BD_All),
    DTU_BD_Up = any(DTU_BD_Up),
    DTU_BD_Down = any(DTU_BD_Down),


    # Smallest DTE FDR in each TSS cluster
    Minimum_DTE_ASD_FDR = if (
      all(is.na(ASD.fdr))
    ) {
      NA_real_
    } else {
      min(ASD.fdr, na.rm = TRUE)
    },

    Minimum_DTE_SCZ_FDR = if (
      all(is.na(SCZ.fdr))
    ) {
      NA_real_
    } else {
      min(SCZ.fdr, na.rm = TRUE)
    },

    Minimum_DTE_BD_FDR = if (
      all(is.na(BD.fdr))
    ) {
      NA_real_
    } else {
      min(BD.fdr, na.rm = TRUE)
    },


    # Smallest DTU FDR in each TSS cluster (find the strongest trascript in a cluster)
    Minimum_DTU_ASD_FDR = if (
      all(is.na(DTU.ASD.FDR))         # Check whether all FDR values are missing
    ) {
      NA_real_
    } else {
      min(DTU.ASD.FDR, na.rm = TRUE)
    },

    Minimum_DTU_SCZ_FDR = if (
      all(is.na(DTU.SCZ.FDR))
    ) {
      NA_real_
    } else {
      min(DTU.SCZ.FDR, na.rm = TRUE)
    },

    Minimum_DTU_BD_FDR = if (
      all(is.na(DTU.BD.FDR))
    ) {
      NA_real_
    } else {
      min(DTU.BD.FDR, na.rm = TRUE)
    },

    .groups = "drop"       # Removing the grouping
  )



# SUMMARIZE TSS CLUSTER BINDING
tss_cluster_binding_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_TSS_Clusters = c(
    nrow(tss_cluster_data),
    nrow(tss_cluster_data),
    nrow(tss_cluster_data)
  ),

  Promoter_Bound_TSS_Clusters = c(
    sum(tss_cluster_data$Promoter_Bound_NPC),
    sum(tss_cluster_data$Promoter_Bound_McClay),
    sum(tss_cluster_data$Promoter_Bound_Forrest)
  ),

  Promoter_Bound_Genes = c(
    length(
      unique(
        tss_cluster_data$ensembl_gene_id[
          tss_cluster_data$Promoter_Bound_NPC
        ]
      )
    ),

    length(
      unique(
        tss_cluster_data$ensembl_gene_id[
          tss_cluster_data$Promoter_Bound_McClay
        ]
      )
    ),

    length(
      unique(
        tss_cluster_data$ensembl_gene_id[
          tss_cluster_data$Promoter_Bound_Forrest
        ]
      )
    )
  )
)

tss_cluster_binding_summary



# PREPARE LABELS FOR ALL FISHER TESTS
TCF4_dataset_names <- c(
  "NPC",
  "McClay",
  "Forrest"
)

TCF4_binding_columns <- c(
  "Promoter_Bound_NPC",
  "Promoter_Bound_McClay",
  "Promoter_Bound_Forrest"
)

expression_features <- c(
  "DTE",
  "DTU"
)

disorder_names <- c(
  "ASD",
  "SCZ",
  "BD"
)

direction_names <- c(
  "All",
  "Up",
  "Down"
)


# RUN GLOBAL TSS-CLUSTER FISHER TESTS
# Question is: Are significant DTE or DTU TSS clusters more likely to have a
# TCF4 peak at their promoter than nonsignificant TSS clusters?
# alternative = "greater" tests enrichment rather than depletion.
global_enrichment_results <- data.frame()  # Create an empty result table: every Fisher result will later be added as a new row.

# This code automatically repeats the same Fisher enrichment test for every combination of:
# Because we have:  - 3 TCF4 datasets
#                   - 2 expression analyses
#                   - 3 disorders
#                   - 3 directions
# so we should do 54 fisher tests

for (                                         # First loop: select the TCF4 dataset
  TCF4_dataset_number in seq_along(TCF4_dataset_names)
) {

  current_TCF4_dataset <- TCF4_dataset_names[
    TCF4_dataset_number
  ]

  current_binding_column <- TCF4_binding_columns[
    TCF4_dataset_number
  ]

  current_bound_status <- tss_cluster_data[[current_binding_column]]


  for (                                        # Second loop: select DTE or DTU
    current_feature in expression_features
  ) {

    for (                                      # Third loop: select the disorder
      current_disorder in disorder_names
    ) {

      for (                                    # Fourth loop: select the direction
        current_direction in direction_names
      ) {

        current_significance_column <- paste(
          current_feature,
          current_disorder,
          current_direction,
          sep = "_"
        )

        current_significant_status <- tss_cluster_data[[current_significance_column]]



        # Construct the Fisher table
        #                         Significant    Not significant
        # TCF4 bound                   A               B
        # Not TCF4 bound               C               D
        # sum(): counts only TRUE values
        current_fisher_table <- matrix(         # A: TCF4-bound and significant
          c(
            sum(
              current_bound_status &
                current_significant_status
            ),

            sum(                                # B: TCF4-bound but not significant
              current_bound_status &
                !current_significant_status
            ),

            sum(                               # C: Not TCF4-bound but significant
              !current_bound_status &
                current_significant_status
            ),

            sum(                               # D: Neither TCF4-bound nor significant
              !current_bound_status &
                !current_significant_status
            )
          ),

          nrow = 2,
          byrow = TRUE
        )

        rownames(
          current_fisher_table
        ) <- c(
          "TCF4_Promoter_Bound",
          "Not_TCF4_Promoter_Bound"
        )

        colnames(
          current_fisher_table
        ) <- c(
          "Significant",
          "Not_Significant"
        )


        # Run Fisher's exact test
        current_fisher_result <- fisher.test(
          current_fisher_table,
          alternative = "greater"
        )


        # Add the result to the result table
        global_enrichment_results <- bind_rows(
          global_enrichment_results,

          data.frame(
            TCF4_Dataset = current_TCF4_dataset,

            Expression_Analysis = current_feature,

            Disorder = current_disorder,

            Direction = current_direction,

            Total_TSS_Clusters = length(
              current_bound_status
            ),

            TCF4_Bound_TSS_Clusters = sum(
              current_bound_status
            ),

            Significant_TSS_Clusters = sum(
              current_significant_status
            ),

            TCF4_Bound_and_Significant_TSS_Clusters = sum(
              current_bound_status &
                current_significant_status
            ),

            Odds_Ratio = as.numeric(
              current_fisher_result$estimate
            ),

            Confidence_Interval_Lower = as.numeric(
              current_fisher_result$conf.int[1]
            ),

            P_Value = current_fisher_result$p.value
          )
        )
      }
    }
  }
}



# CORRECT GLOBAL FISHER P VALUES FOR MULTIPLE TESTING
global_enrichment_results <- global_enrichment_results %>%
  mutate(
    FDR = p.adjust(       # p.adjust(): corrects P values for multiple comparisons
      P_Value,            # correct all values in the P_Value column
      method = "BH"
    )
  ) %>%

  arrange(
    FDR,
    P_Value
  )



# VIEW THE GLOBAL ENRICHMENT RESULTS
global_enrichment_results

global_enrichment_results %>%
  filter(
    Disorder == "SCZ"
  )

global_enrichment_results %>%
  filter(
    Disorder == "SCZ",
    Expression_Analysis == "DTE"
  )

global_enrichment_results %>%
  filter(
    Disorder == "SCZ",
    Expression_Analysis == "DTU"
  )


# WITHIN-GENE PERMUTATION ANALYSIS

# Why do this? Transcripts from the same gene are related and are not independent
# The permutation asks: Within genes that have multiple TSS clusters, are the significant
# clusters more often TCF4-bound than expected if significance were
# randomly distributed among that gene's TSS clusters?
# we want to perform 10,000 permutations
number_of_permutations <- 10000     #create 10,000 random results for each comparison

within_gene_permutation_results <- data.frame()  #Create an empty result table

permutation_test_number <- 0  #This counts how many comparisons have been performed.


for (                                                   # First loop:Select the TCF4 dataset
  TCF4_dataset_number in seq_along(TCF4_dataset_names)
) {

  current_TCF4_dataset <- TCF4_dataset_names[    # select TCF4 dataset
    TCF4_dataset_number
  ]

  current_binding_column <- TCF4_binding_columns[  # select the corresponding binding column
    TCF4_dataset_number
  ]

  current_bound_status <- tss_cluster_data[[ current_binding_column]]

  for (                                      # Second loop: select DTE or DTU
    current_feature in expression_features
  ) {

    for (                                    # Third loop: select the disorder
      current_disorder in disorder_names
    ) {

      for (                                  # Fourth loop: select the direction (All, UP, or DOWN)
        current_direction in direction_names
      ) {

        permutation_test_number <- permutation_test_number + 1

        current_significance_column <- paste(  #create a label like "DTE_SCZ_Down"
          current_feature,
          current_disorder,
          current_direction,
          sep = "_"
        )

        current_significant_status <- tss_cluster_data[[current_significance_column]]  # extract the labeled column (like True, False, True,etc.)


        # Create the within-gene permutation input table
        # This section starts with one row for every TSS cluster
        current_permutation_input <- data.frame(
          Gene_ID = tss_cluster_data$ensembl_gene_id,

          Bound = current_bound_status,

          Significant = current_significant_status
        ) %>%

          group_by(  # group the TSS clusters by gene now
            Gene_ID
          ) %>%

          summarise(
            Number_of_TSS_Clusters = n(),  # Count the TSS clusters in each gene

            Number_of_Bound_TSS_Clusters = sum( # Count TCF4-bound clusters
              Bound
            ),

            Number_of_Significant_TSS_Clusters = sum(   #Count significant clusters
              Significant
            ),

            Observed_Bound_and_Significant = sum(  #Count the observed overlap
              Bound &
                Significant
            ),

            .groups = "drop"
          ) %>%

          filter(       # define our criteria:Keep only informative genes

            # The gene must have at least two TSS clusters
            Number_of_TSS_Clusters >= 2,

            # The gene must have both bound and unbound TSS clusters
            Number_of_Bound_TSS_Clusters > 0,

            Number_of_Bound_TSS_Clusters <
              Number_of_TSS_Clusters,


            # The gene must have both significant and nonsignificant TSS clusters
            Number_of_Significant_TSS_Clusters > 0,

            Number_of_Significant_TSS_Clusters <
              Number_of_TSS_Clusters
          )


        # Run the permutations only if informative genes are available
        if (
          nrow(current_permutation_input) > 0
        ) {

          set.seed(20260723 + permutation_test_number) # makes the random results reproducible
                                                       # Adding permutation_test_number gives each of the 54 comparisons a different seed
          permuted_overlap_totals <- numeric(          # Create storage for the random overlaps: Because: number_of_permutations = 10000
                                                       # this creates 10,000 zeros: 0,0,0,0,0,0,0,0, etc. (Each position represents one complete permutation)
            number_of_permutations
          )



          # Process every informative gene: The inner gene loop examines each informative gene separately
          # For every gene, randomly distribute the significant label among its TSS clusters
          #
          # rhyper() performs the random sampling while preserving:
          #   - number of TSS clusters
          #   - number of bound clusters
          #   - number of significant clusters
          # this function asks If the one significant label were randomly placed among these total TSS clusters,
          # how often would it fall on one of the TCF4-bound clusters?
          # (This is mathematically equivalent to randomly moving the significance label while keeping TCF4 binding fixed.)

          for (
            current_gene_row in seq_len(
              nrow(
                current_permutation_input
              )
            )
          ) {

            total_clusters_in_gene <-
              current_permutation_input$Number_of_TSS_Clusters[
                current_gene_row
              ]

            significant_clusters_in_gene <-
              current_permutation_input$Number_of_Significant_TSS_Clusters[
                current_gene_row
              ]

            nonsignificant_clusters_in_gene <-
              total_clusters_in_gene -
              significant_clusters_in_gene

            bound_clusters_in_gene <-
              current_permutation_input$Number_of_Bound_TSS_Clusters[
                current_gene_row
              ]


            permuted_overlap_totals <-
              permuted_overlap_totals +
              rhyper(
                number_of_permutations,           #generate 10,000 random results

                m = significant_clusters_in_gene,

                n = nonsignificant_clusters_in_gene,

                k = bound_clusters_in_gene
              )
          }


          # Calculate the observed overlap
          observed_overlap_total <- sum(       #This adds the real bound-and-significant overlaps across all informative genes.
            current_permutation_input$
              Observed_Bound_and_Significant
          )


          # Calculate the empirical permutation P value
          # The +1 correction prevents a P value of exactly zero.
          empirical_P_value <- (
            sum(                               #This counts how many random results were equal to or greater than the real overlap.
              permuted_overlap_totals >=
                observed_overlap_total) +1) / (number_of_permutations +1)



          # Add the permutation result to the result table
          within_gene_permutation_results <- bind_rows(
            within_gene_permutation_results,

            data.frame(
              TCF4_Dataset = current_TCF4_dataset,

              Expression_Analysis = current_feature,

              Disorder = current_disorder,

              Direction = current_direction,

              Informative_Genes = nrow(
                current_permutation_input
              ),

              Informative_TSS_Clusters = sum(
                current_permutation_input$
                  Number_of_TSS_Clusters
              ),

              Observed_TCF4_Bound_and_Significant = observed_overlap_total,

              Mean_Random_TCF4_Bound_and_Significant = mean(
                permuted_overlap_totals
              ),

              Standard_Deviation_of_Random_Overlap = sd(
                permuted_overlap_totals
              ),

              Empirical_P_Value = empirical_P_value
            )
          )

        } else {

          # Add an NA result when no informative genes are available
          # It's like handling comparisons with no informative genes
          within_gene_permutation_results <- bind_rows(
            within_gene_permutation_results,

            data.frame(
              TCF4_Dataset = current_TCF4_dataset,

              Expression_Analysis = current_feature,

              Disorder = current_disorder,

              Direction = current_direction,

              Informative_Genes = 0,

              Informative_TSS_Clusters = 0,

              Observed_TCF4_Bound_and_Significant = NA_real_,

              Mean_Random_TCF4_Bound_and_Significant = NA_real_,

              Standard_Deviation_of_Random_Overlap = NA_real_,

              Empirical_P_Value = NA_real_
            )
          )
        }
      }
    }
  }
}


# CORRECT PERMUTATION P VALUES FOR MULTIPLE TESTING
within_gene_permutation_results <- within_gene_permutation_results %>%
  mutate(
    Empirical_FDR = p.adjust(
      Empirical_P_Value,
      method = "BH"
    )
  ) %>%

  arrange(
    Empirical_FDR,
    Empirical_P_Value
  )


# VIEW ALL PERMUTATION RESULTS
within_gene_permutation_results


# VIEW THE SCHIZOPHRENIA PERMUTATION RESULTS
within_gene_permutation_results %>%
  filter(
    Disorder == "SCZ"
  )



# VIEW SCHIZOPHRENIA DTE PERMUTATION RESULTS
within_gene_permutation_results %>%
  filter(
    Disorder == "SCZ",
    Expression_Analysis == "DTE"
  )



# VIEW SCHIZOPHRENIA DTU PERMUTATION RESULTS
# These results are the more direct test of isoform usage
within_gene_permutation_results %>%
  filter(
    Disorder == "SCZ",
    Expression_Analysis == "DTU"
  )



# FIND INDIVIDUAL TRANSCRIPTS THAT:
#   1. Have a promoter overlapping at least one TCF4 peak
#   2. Have significant DTE or DTU in ASD, SCZ or BD
candidate_TCF4_bound_transcripts <- transcript_data %>%
  filter(
    Promoter_Bound_NPC |
      Promoter_Bound_McClay |
      Promoter_Bound_Forrest,

    DTE_ASD_All |
      DTE_SCZ_All |
      DTE_BD_All |
      DTU_ASD_All |
      DTU_SCZ_All |
      DTU_BD_All
  ) %>%

  mutate(
    Best_FDR = pmin(
      ASD.fdr,
      SCZ.fdr,
      BD.fdr,
      DTU.ASD.FDR,
      DTU.SCZ.FDR,
      DTU.BD.FDR,
      na.rm = TRUE
    )
  ) %>%

  arrange(
    Best_FDR
  ) %>%

  select(
    ensembl_gene_id,
    external_gene_id,

    ensembl_transcript_id,
    external_transcript_id,

    chromosome,
    transcript_TSS,
    transcript_start,
    transcript_end,
    strand,

    Promoter_Bound_NPC,
    Promoter_Bound_McClay,
    Promoter_Bound_Forrest,

    ASD.log2FC,
    ASD.fdr,

    SCZ.log2FC,
    SCZ.fdr,

    BD.log2FC,
    BD.fdr,

    DTU.ASD.Value,
    DTU.ASD.FDR,

    DTU.SCZ.Value,
    DTU.SCZ.FDR,

    DTU.BD.Value,
    DTU.BD.FDR,

    Best_FDR
  )



# 34. VIEW THE CANDIDATE TRANSCRIPTS
candidate_TCF4_bound_transcripts

head(
  candidate_TCF4_bound_transcripts,
  50
)


# VIEW NPC-BOUND SCHIZOPHRENIA DTE TRANSCRIPTS
NPC_bound_SCZ_DTE_transcripts <- transcript_data %>%
  filter(
    Promoter_Bound_NPC,
    DTE_SCZ_All
  ) %>%

  arrange(
    SCZ.fdr
  ) %>%

  select(
    ensembl_gene_id,
    external_gene_id,

    ensembl_transcript_id,
    external_transcript_id,

    chromosome,
    transcript_TSS,

    SCZ.log2FC,
    SCZ.fdr,

    Promoter_Bound_NPC
  )

NPC_bound_SCZ_DTE_transcripts


# VIEW NPC-BOUND SCHIZOPHRENIA DTU TRANSCRIPTS
# This is the more direct list for transcript-usage changes.
NPC_bound_SCZ_DTU_transcripts <- transcript_data %>%
  filter(
    Promoter_Bound_NPC,
    DTU_SCZ_All
  ) %>%

  arrange(
    DTU.SCZ.FDR
  ) %>%

  select(
    ensembl_gene_id,
    external_gene_id,

    ensembl_transcript_id,
    external_transcript_id,

    chromosome,
    transcript_TSS,

    DTU.SCZ.Value,
    DTU.SCZ.FDR,

    Promoter_Bound_NPC
  )

NPC_bound_SCZ_DTU_transcripts
list(NPC_bound_SCZ_DTU_transcripts$external_gene_id)
View(NPC_bound_SCZ_DTU_transcripts)
write.csv(NPC_bound_SCZ_DTU_transcripts, "NPC_bound_SCZ_DTU_transcripts.csv")

# Identify genes contributing to NPC-bound SCZ DTU-Down enrichment
# For every gene, calculate:
# - total number of TSS clusters;
# - number of NPC TCF4-bound TSS clusters;
# - number of SCZ DTU-Down TSS clusters;
# - observed overlap between binding and DTU Down.
NPC_SCZ_DTU_Down_gene_contributions <- tss_cluster_data %>%
  group_by(
    ensembl_gene_id,
    external_gene_id
  ) %>%

  summarise(

    # Count all TSS clusters belonging to the gene.
    Total_TSS_Clusters = n(),

    # Count NPC TCF4-bound TSS clusters.
    NPC_Bound_TSS_Clusters = sum(
      Promoter_Bound_NPC,
      na.rm = TRUE
    ),

    # Count TSS clusters containing at least one significant
    # SCZ DTU-Down transcript.
    SCZ_DTU_Down_TSS_Clusters = sum(
      DTU_SCZ_Down,
      na.rm = TRUE
    ),

    # Count TSS clusters that are both NPC-bound and SCZ DTU Down.
    Observed_Overlap = sum(
      Promoter_Bound_NPC &
        DTU_SCZ_Down,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  # Retain exactly the genes eligible for the within-gene permutation.
  filter(

    # The gene must contain at least two TSS clusters.
    Total_TSS_Clusters >= 2,

    # It must contain bound and unbound TSS clusters.
    NPC_Bound_TSS_Clusters > 0,
    NPC_Bound_TSS_Clusters < Total_TSS_Clusters,

    # It must contain significant and nonsignificant DTU-Down clusters.
    SCZ_DTU_Down_TSS_Clusters > 0,
    SCZ_DTU_Down_TSS_Clusters < Total_TSS_Clusters
  ) %>%

  mutate(

    # Expected overlap within each gene under random alignment.
    #
    # For example, if a gene contains:
    # 10 total TSS clusters,
    # 4 bound clusters,
    # 2 DTU-Down clusters,
    #
    # expected overlap = 4 × 2 / 10 = 0.8
    Expected_Overlap = (
      NPC_Bound_TSS_Clusters *
        SCZ_DTU_Down_TSS_Clusters
    ) / Total_TSS_Clusters,

    # Positive values mean the gene contributed more overlap
    # than expected under random within-gene alignment.
    Excess_Overlap = Observed_Overlap -
      Expected_Overlap,

    # Values above 1 indicate more overlap than expected.
    Observed_Expected_Ratio = Observed_Overlap /
      Expected_Overlap
  ) %>%

  # Rank genes by how much observed overlap exceeds expectation.
  arrange(
    desc(Excess_Overlap),
    desc(Observed_Expected_Ratio)
  )


NPC_SCZ_DTU_Down_gene_contributions

# This should return 158 informative genes.
nrow(NPC_SCZ_DTU_Down_gene_contributions)

# This should return 1,576 TSS clusters.
sum(NPC_SCZ_DTU_Down_gene_contributions$Total_TSS_Clusters)

# This should return 100 observed overlap clusters.
sum(NPC_SCZ_DTU_Down_gene_contributions$Observed_Overlap)

# This should be close to the random mean of 79.25.
sum(NPC_SCZ_DTU_Down_gene_contributions$Expected_Overlap)

# First, retain all TSS clusters belonging to the 158 informative genes.
#
# semi_join() keeps rows from tss_cluster_data whose gene identifiers
# are present in the informative-gene table
NPC_SCZ_DTU_Down_informative_TSS_clusters <- tss_cluster_data %>%
  semi_join(
    NPC_SCZ_DTU_Down_gene_contributions %>%
      select(
        ensembl_gene_id
      ),

    by = "ensembl_gene_id"
  )

# Now retain the TSS clusters that created the observed overlap:
# NPC-bound and significant SCZ DTU Down.
NPC_SCZ_DTU_Down_overlap_TSS_clusters <-
  NPC_SCZ_DTU_Down_informative_TSS_clusters %>%

  filter(
    Promoter_Bound_NPC,
    DTU_SCZ_Down
  ) %>%

  left_join(
    NPC_SCZ_DTU_Down_gene_contributions %>%
      select(
        ensembl_gene_id,
        Expected_Overlap,
        Observed_Overlap,
        Excess_Overlap,
        Observed_Expected_Ratio
      ),

    by = "ensembl_gene_id"
  ) %>%

  arrange(
    desc(Excess_Overlap),
    Minimum_DTU_SCZ_FDR
  ) %>%

  select(
    ensembl_gene_id,
    external_gene_id,

    TSS_Cluster_ID,
    chromosome,
    strand,
    transcript_TSS,

    Number_of_Transcripts,
    Ensembl_Transcript_IDs,
    Transcript_Names,

    Promoter_Bound_NPC,
    Promoter_Bound_McClay,
    Promoter_Bound_Forrest,

    DTU_SCZ_Down,
    Minimum_DTU_SCZ_FDR,

    Expected_Overlap,
    Observed_Overlap,
    Excess_Overlap,
    Observed_Expected_Ratio
  )


# This should return exactly 100 TSS clusters.
nrow(NPC_SCZ_DTU_Down_overlap_TSS_clusters)


# Display the exact overlap clusters.
NPC_SCZ_DTU_Down_overlap_TSS_clusters
View(NPC_SCZ_DTU_Down_overlap_TSS_clusters)
write.csv(NPC_SCZ_DTU_Down_overlap_TSS_clusters,"NPC_SCZ_DTU_Down_overlap_TSS_clusters.csv")
