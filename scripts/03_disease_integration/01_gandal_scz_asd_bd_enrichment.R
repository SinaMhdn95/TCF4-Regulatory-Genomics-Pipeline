# TCF4 genes vs Gandal's gene dataset
#https://pmc.ncbi.nlm.nih.gov/articles/PMC6443102/#SD8

# Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(tidydr)
library(stringr)
library(readr)
library(tibble)

# Import Gandal's dataset
SCZ_AD_BPD_data <- read_excel("Transcriptome_SCZ_AD_BPD.xlsx", sheet = "DGE")

# Import my TCF4 datsets
# all peaks genes
gene_all <- read.csv("All_TCF4_peak_genes.csv")
gene_all <- gene_all$Gene
length(gene_all)

# E-box peaks genes
gene_motif <- read.csv("EBOX_peak_genes.csv")
gene_motif <- gene_motif$Gene
length(gene_motif)

# Promoter + E-box genes
promoter_ebox_genes <- read.csv("Promoter_EBOX_genes.csv")
promoter_ebox_genes <- promoter_ebox_genes$Gene
length(promoter_ebox_genes)

# We need a background gene list (reference) including all genes
background_gene_list <-  SCZ_AD_BPD_data$gene_name
n_distinct(background_gene_list)

# Make a significant gene list for each disease
SCZ_sig <- SCZ_AD_BPD_data %>%  # Significant SCZ genes based on FDR < 0.05
  filter(SCZ.fdr < 0.05)

SCZ_sig_gene_list <- SCZ_sig$gene_name # Only gene names
n_distinct(SCZ_sig_gene_list)

ASD_sig <- SCZ_AD_BPD_data %>%  # Significant ASD genes based on FDR < 0.05
  filter(ASD.fdr < 0.05)

ASD_sig_gene_list <- ASD_sig$gene_name # Only gene names
n_distinct(ASD_sig_gene_list)

BPD_sig <- SCZ_AD_BPD_data %>%  # Significant BPD genes based on FDR < 0.05
  filter(BD.fdr < 0.05)

BPD_sig_gene_list <- BPD_sig$gene_name # Only gene names
n_distinct(BPD_sig_gene_list)


# Create logical vectors for our TCF4 gene sets
is_all_gene <- background_gene_list %in% gene_all
sum(is_all_gene)

is_motif_gene <- background_gene_list %in% gene_motif
sum(is_motif_gene)

is_promoter_gene <- background_gene_list %in% promoter_ebox_genes
sum(is_promoter_gene)

# Convert significant gene lists into logical vectors aligned with background
is_SCZ_all <- background_gene_list %in% SCZ_sig_gene_list
sum(is_SCZ_all)

is_ASD_all <- background_gene_list %in% ASD_sig_gene_list
sum(is_ASD_all)

is_BPD_all <- background_gene_list %in% BPD_sig_gene_list
sum(is_BPD_all)

# Fisher test: Are TCF4_all genes enriched in significant genes?
# SCZ
fisher_SCZ <- fisher.test(table(is_all_gene, is_SCZ_all), alternative = "greater")
fisher_SCZ
table(is_all_gene, is_SCZ_all)

# ASD
fisher_ASD <- fisher.test(table(is_all_gene, is_ASD_all), alternative = "greater")
fisher_ASD
table(is_all_gene, is_ASD_all)

# BPD
fisher_BPD <- fisher.test(table(is_all_gene, is_BPD_all), alternative = "greater")
fisher_BPD
table(is_all_gene, is_BPD_all)

# Fisher test: Are TCF4_motif genes enriched?
# SCZ
fisher_SCZ_motif <- fisher.test(table(is_motif_gene, is_SCZ_all), alternative = "greater")
fisher_SCZ_motif
table(is_motif_gene, is_SCZ_all)

# ASD
fisher_ASD_motif <- fisher.test(table(is_motif_gene, is_ASD_all), alternative = "greater")
fisher_ASD_motif
table(is_motif_gene, is_ASD_all)

# BPD
fisher_BPD_motif <- fisher.test(table(is_motif_gene, is_BPD_all), alternative = "greater")
fisher_BPD_motif
table(is_motif_gene, is_BPD_all)

# Fisher test: Are TCF4_promoter_ebox genes enriched?
# SCZ
fisher_SCZ_promoter <- fisher.test(table(is_promoter_gene, is_SCZ_all), alternative = "greater")
fisher_SCZ_promoter
table(is_promoter_gene, is_SCZ_all)

# ASD
fisher_ASD_promoter <- fisher.test(table(is_promoter_gene, is_ASD_all), alternative = "greater")
fisher_ASD_promoter
table(is_promoter_gene, is_ASD_all)

# BPD
fisher_BPD_promoter <- fisher.test(table(is_promoter_gene, is_BPD_all), alternative = "greater")
fisher_BPD_promoter
table(is_promoter_gene, is_BPD_all)

# Separate downregulation and upregulation significat genes across diseases
#SCZ
SCZ_sig_up <- SCZ_sig %>%     # Upregulated significant genes
  filter(SCZ.log2FC > 0)

SCZ_sig_up_gene <- SCZ_sig_up$gene_name  # Extract only gene names

SCZ_sig_down <- SCZ_sig %>%   # Downregulated significant genes
  filter(SCZ.log2FC < 0)

SCZ_sig_down_gene <- SCZ_sig_down$gene_name # Extract only gene names

#ASD
ASD_sig_up <- ASD_sig %>%     # Upregulated significant genes
  filter(ASD.log2FC > 0)

ASD_sig_up_gene <- ASD_sig_up$gene_name   # Extract only gene names

ASD_sig_down <- ASD_sig %>%   # Downregulated significant genes
  filter(ASD.log2FC < 0)

ASD_sig_down_gene <- ASD_sig_down$gene_name  # Extract only gene names


#BPD
BPD_sig_up <- BPD_sig %>%    # Upregulated significant genes
  filter(BD.log2FC > 0)

BPD_sig_up_gene <- BPD_sig_up$gene_name # Extract only gene names

BPD_sig_down <- BPD_sig %>%  # Downregulated significant genes
  filter(BD.log2FC < 0)

BPD_sig_down_gene <- BPD_sig_down$gene_name # Extract only gene names

# Create logical vectors for upregulated and downregulated significant genes
is_SCZ_up <- background_gene_list %in% SCZ_sig_up_gene
is_SCZ_down <- background_gene_list %in% SCZ_sig_down_gene

is_ASD_up <- background_gene_list %in% ASD_sig_up_gene
is_ASD_down <- background_gene_list %in% ASD_sig_down_gene

is_BPD_up <- background_gene_list %in% BPD_sig_up_gene
is_BPD_down <- background_gene_list %in% BPD_sig_down_gene

# Fisher test: TCF4_all genes vs UP and DOWN genes
# SCZ up
fisher_SCZ_up <- fisher.test(table(is_all_gene, is_SCZ_up), alternative = "greater")
fisher_SCZ_up
table(is_all_gene, is_SCZ_up)

# SCZ down
fisher_SCZ_down <- fisher.test(table(is_all_gene, is_SCZ_down), alternative = "greater")
fisher_SCZ_down
table(is_all_gene, is_SCZ_down)


# ASD up
fisher_ASD_up <- fisher.test(table(is_all_gene, is_ASD_up), alternative = "greater")
fisher_ASD_up
table(is_all_gene, is_ASD_up)

# ASD down
fisher_ASD_down <- fisher.test(table(is_all_gene, is_ASD_down), alternative = "greater")
fisher_ASD_down
table(is_all_gene, is_ASD_down)


# BPD up
fisher_BPD_up <- fisher.test(table(is_all_gene, is_BPD_up), alternative = "greater")
fisher_BPD_up
table(is_all_gene, is_BPD_up)

# BPD down
fisher_BPD_down <- fisher.test(table(is_all_gene, is_BPD_down), alternative = "greater")
fisher_BPD_down
table(is_all_gene, is_BPD_down)


# Fisher test: TCF4_motif genes vs UP and DOWN genes
# SCZ up
fisher_SCZ_motif_up <- fisher.test(table(is_motif_gene, is_SCZ_up), alternative = "greater")
fisher_SCZ_motif_up
table(is_motif_gene, is_SCZ_up)

# SCZ down
fisher_SCZ_motif_down <- fisher.test(table(is_motif_gene, is_SCZ_down), alternative = "greater")
fisher_SCZ_motif_down
table(is_motif_gene, is_SCZ_down)


# ASD up
fisher_ASD_motif_up <- fisher.test(table(is_motif_gene, is_ASD_up), alternative = "greater")
fisher_ASD_motif_up
table(is_motif_gene, is_ASD_up)

# ASD down
fisher_ASD_motif_down <- fisher.test(table(is_motif_gene, is_ASD_down), alternative = "greater")
fisher_ASD_motif_down
table(is_motif_gene, is_ASD_down)


# BPD up
fisher_BPD_motif_up <- fisher.test(table(is_motif_gene, is_BPD_up), alternative = "greater")
fisher_BPD_motif_up
table(is_motif_gene, is_BPD_up)

# BPD down
fisher_BPD_motif_down <- fisher.test(table(is_motif_gene, is_BPD_down), alternative = "greater")
fisher_BPD_motif_down
table(is_motif_gene, is_BPD_down)


# Fisher test: TCF4_promoter_ebox genes vs UP and DOWN genes
# SCZ up
fisher_SCZ_promoter_up <- fisher.test(table(is_promoter_gene, is_SCZ_up), alternative = "greater")
fisher_SCZ_promoter_up
table(is_promoter_gene, is_SCZ_up)

# SCZ down
fisher_SCZ_promoter_down <- fisher.test(table(is_promoter_gene, is_SCZ_down), alternative = "greater")
fisher_SCZ_promoter_down
table(is_promoter_gene, is_SCZ_down)


# ASD up
fisher_ASD_promoter_up <- fisher.test(table(is_promoter_gene, is_ASD_up), alternative = "greater")
fisher_ASD_promoter_up
table(is_promoter_gene, is_ASD_up)

# ASD down
fisher_ASD_promoter_down <- fisher.test(table(is_promoter_gene, is_ASD_down), alternative = "greater")
fisher_ASD_promoter_down
table(is_promoter_gene, is_ASD_down)


# BPD up
fisher_BPD_promoter_up <- fisher.test(table(is_promoter_gene, is_BPD_up), alternative = "greater")
fisher_BPD_promoter_up
table(is_promoter_gene, is_BPD_up)

# BPD down
fisher_BPD_promoter_down <- fisher.test(table(is_promoter_gene, is_BPD_down), alternative = "greater")
fisher_BPD_promoter_down
table(is_promoter_gene, is_BPD_down)

# Find the overlaps for SCZ between TCF4 datset and Gandal's dataset
# TCF4 all peaks + SCZ upregulated genes
SCZ_all_up_overlap <- background_gene_list[
  is_all_gene & is_SCZ_up]

# TCF4 all peaks + SCZ downregulated genes
SCZ_all_down_overlap <- background_gene_list[
  is_all_gene & is_SCZ_down]


# TCF4 motif peaks + SCZ upregulated genes
SCZ_motif_up_overlap <- background_gene_list[
  is_motif_gene & is_SCZ_up]

# TCF4 motif peaks + SCZ downregulated genes
SCZ_motif_down_overlap <- background_gene_list[
  is_motif_gene & is_SCZ_down]


# TCF4 promoter + E-box peaks + SCZ upregulated genes
SCZ_promoter_up_overlap <- background_gene_list[
  is_promoter_gene & is_SCZ_up]

# TCF4 promoter + E-box peaks + SCZ downregulated genes
SCZ_promoter_down_overlap <- background_gene_list[
  is_promoter_gene & is_SCZ_down]

length(SCZ_all_up_overlap)
length(SCZ_all_down_overlap)

# Summary of overlaps
SCZ_overlap_summary <- data.frame(

  TCF4_Set = c(
    "All peaks",
    "Motif peaks",
    "Promoter E-box"
  ),

  Upregulated_Overlap = c(
    length(SCZ_all_up_overlap),
    length(SCZ_motif_up_overlap),
    length(SCZ_promoter_up_overlap)
  ),

  Downregulated_Overlap = c(
    length(SCZ_all_down_overlap),
    length(SCZ_motif_down_overlap),
    length(SCZ_promoter_down_overlap)
  )
)

SCZ_overlap_summary

# Make a comprehensive overlap dataset for further analysis
# Comprehensive SCZ overlap dataset
SCZ_overlap_all_up <- data.frame(
  gene = SCZ_all_up_overlap,
  TCF4_Set = "All peaks",
  SCZ_Direction = "Up"
)

SCZ_overlap_all_down <- data.frame(
  gene = SCZ_all_down_overlap,
  TCF4_Set = "All peaks",
  SCZ_Direction = "Down"
)

SCZ_overlap_motif_up <- data.frame(
  gene = SCZ_motif_up_overlap,
  TCF4_Set = "Motif peaks",
  SCZ_Direction = "Up"
)

SCZ_overlap_motif_down <- data.frame(
  gene = SCZ_motif_down_overlap,
  TCF4_Set = "Motif peaks",
  SCZ_Direction = "Down"
)

SCZ_overlap_promoter_up <- data.frame(
  gene = SCZ_promoter_up_overlap,
  TCF4_Set = "Promoter E-box",
  SCZ_Direction = "Up"
)

SCZ_overlap_promoter_down <- data.frame(
  gene = SCZ_promoter_down_overlap,
  TCF4_Set = "Promoter E-box",
  SCZ_Direction = "Down"
)

SCZ_TCF4_overlap_long <- bind_rows(
  SCZ_overlap_all_up,
  SCZ_overlap_all_down,
  SCZ_overlap_motif_up,
  SCZ_overlap_motif_down,
  SCZ_overlap_promoter_up,
  SCZ_overlap_promoter_down
)

SCZ_TCF4_overlap_long <- SCZ_TCF4_overlap_long %>%
  distinct(gene, TCF4_Set, SCZ_Direction) %>%
  arrange(TCF4_Set, SCZ_Direction, gene)

# Add SCZ log2FC and FDR
SCZ_TCF4_overlap_long <- SCZ_TCF4_overlap_long %>%
  left_join(
    SCZ_sig %>%
      select(gene = gene_name, SCZ.log2FC, SCZ.fdr),
    by = "gene"
  )

View(SCZ_TCF4_overlap_long)

# Save the overlap dataset
write.csv(
  SCZ_TCF4_overlap_long,
  "SCZ_TCF4_overlap_Ganda.csv",
  row.names = FALSE
)
