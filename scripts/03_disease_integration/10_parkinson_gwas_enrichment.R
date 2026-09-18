#Enrichment analysis using Parkinson's dataset
install.packages("tidyr")
# Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(tidydr)
library(stringr)
library(readr)
library(tibble)

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

# Import new reference
#New UCSC style reference, transcript level
final_ref_gene <- read.csv("final_ref_gene_2026_filtered.csv")

# Extract the unique gene symbols from the cleaned final reference gene list to create a
# comprehensive list of all genes in the genome for our analysis.
ref_genes <- final_ref_gene$gene
ref_genes <- unique(ref_genes[!is.na(ref_genes) & ref_genes != ""])
length(ref_genes)

# Create logical vectors for our TCF4 gene sets
# For each gene in the genome, check if it’s in our TCF4 list.
is_gene_all <- ref_genes %in% gene_all
is_gene_motif <- ref_genes %in% gene_motif  #Which genes have TCF4 + E-box
is_promoter_ebox <- ref_genes %in% promoter_ebox_genes  #Which genes have TCF4 + E-box + promoter

sum(is_gene_all)
sum(is_gene_motif)
sum(is_promoter_ebox)

# Load external datasets
# Parkinson trait from GWAS catalog (MONDO_0005180)
# Link: https://www.ebi.ac.uk/gwas/efotraits/MONDO_0005180
# Import dataset
PD_GWAScat_dataset <- read_tsv("Data integration/GWAS_Cat_Mono5180_PKD.tsv")

# Summarize daaset
PD_GWAScat_dataset %>%
  summarise(
    total_rows = n(),
    rows_with_genes = sum(!is.na(mappedGenes) & mappedGenes != "-"),
    rows_without_genes = sum(is.na(mappedGenes) | mappedGenes == "-")
  )

# Clean dataset (unique genes)
PD_GWAScat_genes <- PD_GWAScat_dataset %>%
  select(mappedGenes) %>%
  filter(!is.na(mappedGenes), mappedGenes != "-") %>%
  tidyr::separate_rows(mappedGenes, sep = ",") %>%
  mutate(gene = stringr::str_trim(mappedGenes)) %>%
  filter(gene != "") %>%
  select(gene) %>%
  distinct()

# Extract gene names
PD_GWAScat_gene_list <- PD_GWAScat_genes$gene

length(PD_GWAScat_gene_list)
head(PD_GWAScat_gene_list)

# Make sure my PD genes are restricted to your universe
PD_genes_clean <- intersect(PD_GWAScat_gene_list, ref_genes)
length(PD_genes_clean)

# Create logical vector for PD
is_PD <- ref_genes %in% PD_genes_clean
sum(is_PD)

# Run enrichment (Fisher tests)
# All TCF4 peaks vs PD
fisher_all_PD <- fisher.test(
  table(is_gene_all, is_PD),
  alternative = "greater"
)

fisher_all_PD

# Motif (E-box) vs PD
fisher_motif_PD <- fisher.test(
  table(is_gene_motif, is_PD),
  alternative = "greater"
)

fisher_motif_PD

# Promoter + E-box vs PD
fisher_promoter_PD <- fisher.test(
  table(is_promoter_ebox, is_PD),
  alternative = "greater"
)

fisher_promoter_PD

# Extract overlapping genes
overlap_all_PD <- ref_genes[is_gene_all & is_PD]
overlap_motif_PD <- ref_genes[is_gene_motif & is_PD]
overlap_promoter_PD <- ref_genes[is_promoter_ebox & is_PD]

length(overlap_all_PD)
length(overlap_motif_PD)
length(overlap_promoter_PD)

# Summary table
PD_summary <- data.frame(
  Category = c("All Peaks", "Motif", "Promoter + E-box"),
  Overlap = c(length(overlap_all_PD),
              length(overlap_motif_PD),
              length(overlap_promoter_PD)),
  P_value = c(fisher_all_PD$p.value,
              fisher_motif_PD$p.value,
              fisher_promoter_PD$p.value)
)

PD_summary
