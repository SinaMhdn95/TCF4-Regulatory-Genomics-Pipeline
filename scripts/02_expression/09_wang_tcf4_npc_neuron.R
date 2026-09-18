# Paper:
# The RNA-seq analysis identified:
#Day 3 NPCs
# 4891 DE genes
# 2330 upregulated
# 2561 downregulated
#Day 14 Glut_Ns
# 3152 DE genes
# 1862 upregulated
# 1290 downregulated

# TCF4 knockdown enrichment analysis
# Goal:
# 1. Test whether TCF4-bound genes are enriched among KD genes
# 2. Compare NPC Day 3 vs glutamatergic neuron Day 14
# 3. Compare UP vs DOWN genes
# 4. Run pathway enrichment on overlapping genes

# Paper link:https://pubmed.ncbi.nlm.nih.gov/35965434/

# Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(tidydr)
library(stringr)
library(readr)
library(tibble)
library(tidyverse)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

# Import data
NPC_KD <- readxl::read_xlsx("KD_genes_Decon_table_s7.xlsx" , sheet = "Day3")
Glut_KD <- readxl::read_xlsx("KD_genes_Decon_table_s7.xlsx" , sheet = "Day14")

# Import TCF4 gene sets
gene_all_df <- read.csv("All_TCF4_peak_genes.csv")
gene_motif_df <- read.csv("EBOX_peak_genes.csv")
promoter_ebox_df <- read.csv("Promoter_EBOX_genes.csv")

# Import background gene universe
background_gene <- read.csv("final_ref_gene_2026_filtered.csv")

# Clean background gene universe
# Use the same background for all Fisher tests.
# This should represent all genes that could have been detected.
allgenes <- background_gene %>%
  filter(!is.na(gene)) %>%
  filter(gene != "") %>%
  pull(gene) %>%
  unique()

length(allgenes)


# Clean TCF4 gene sets
gene_all <- gene_all_df %>%
  filter(!is.na(Gene)) %>%
  filter(Gene != "") %>%
  pull(Gene) %>%
  unique()

gene_motif <- gene_motif_df %>%
  filter(!is.na(Gene)) %>%
  filter(Gene != "") %>%
  pull(Gene) %>%
  unique()

promoter_ebox_genes <- promoter_ebox_df %>%
  filter(!is.na(Gene)) %>%
  filter(Gene != "") %>%
  pull(Gene) %>%
  unique()


# Restrict TCF4 gene sets to background universe
# Do not allow genes outside the background to enter Fisher tests.
gene_all <- intersect(gene_all, allgenes)
gene_motif <- intersect(gene_motif, allgenes)
promoter_ebox_genes <- intersect(promoter_ebox_genes, allgenes)


# 8. Clean knockdown datasets
# The column name is "-LogFC", but biologically it appears to be logFC.
# Positive values are treated as UP.
# Negative values are treated as DOWN.
NPC_KD_clean <- NPC_KD %>%
  filter(!is.na(Gene)) %>%
  filter(Gene != "") %>%
  filter(!is.na(`-LogFC`)) %>%
  distinct(Gene, .keep_all = TRUE)

Glut_KD_clean <- Glut_KD %>%
  filter(!is.na(Gene)) %>%
  filter(Gene != "") %>%
  filter(!is.na(`-LogFC`)) %>%
  distinct(Gene, .keep_all = TRUE)


# Create UP and DOWN gene lists
NPC_up <- NPC_KD_clean %>%
  filter(`-LogFC` > 0) %>%
  pull(Gene) %>%
  unique() %>%
  intersect(allgenes)

NPC_down <- NPC_KD_clean %>%
  filter(`-LogFC` < 0) %>%
  pull(Gene) %>%
  unique() %>%
  intersect(allgenes)

Glut_up <- Glut_KD_clean %>%
  filter(`-LogFC` > 0) %>%
  pull(Gene) %>%
  unique() %>%
  intersect(allgenes)

Glut_down <- Glut_KD_clean %>%
  filter(`-LogFC` < 0) %>%
  pull(Gene) %>%
  unique() %>%
  intersect(allgenes)


# Quick sanity check
length(gene_all)
length(gene_motif)
length(promoter_ebox_genes)

length(NPC_up)
length(NPC_down)
length(Glut_up)
length(Glut_down)

# Fisher test: All TCF4 peaks vs NPC UP
is_tcf4 <- allgenes %in% gene_all
is_query <- allgenes %in% NPC_up

tab_all_npc_up <- table(is_tcf4, is_query)

fisher_all_npc_up <- fisher.test(
  tab_all_npc_up,
  alternative = "greater"
)

res_all_npc_up <- data.frame(
  TCF4_Set = "All_TCF4_peaks",
  KD_Set = "NPC_UP",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_all_npc_up$estimate),
  P_Value = fisher_all_npc_up$p.value
)

# Fisher test: All TCF4 peaks vs NPC DOWN
is_tcf4 <- allgenes %in% gene_all
is_query <- allgenes %in% NPC_down

tab_all_npc_down <- table(is_tcf4, is_query)

fisher_all_npc_down <- fisher.test(
  tab_all_npc_down,
  alternative = "greater"
)

res_all_npc_down <- data.frame(
  TCF4_Set = "All_TCF4_peaks",
  KD_Set = "NPC_DOWN",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_all_npc_down$estimate),
  P_Value = fisher_all_npc_down$p.value
)

# Fisher test: All TCF4 peaks vs Glut UP
is_tcf4 <- allgenes %in% gene_all
is_query <- allgenes %in% Glut_up

tab_all_glut_up <- table(is_tcf4, is_query)

fisher_all_glut_up <- fisher.test(
  tab_all_glut_up,
  alternative = "greater"
)

res_all_glut_up <- data.frame(
  TCF4_Set = "All_TCF4_peaks",
  KD_Set = "Glut_UP",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_all_glut_up$estimate),
  P_Value = fisher_all_glut_up$p.value
)

# Fisher test: All TCF4 peaks vs Glut DOWN
is_tcf4 <- allgenes %in% gene_all
is_query <- allgenes %in% Glut_down

tab_all_glut_down <- table(is_tcf4, is_query)

fisher_all_glut_down <- fisher.test(
  tab_all_glut_down,
  alternative = "greater"
)

res_all_glut_down <- data.frame(
  TCF4_Set = "All_TCF4_peaks",
  KD_Set = "Glut_DOWN",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_all_glut_down$estimate),
  P_Value = fisher_all_glut_down$p.value
)


# Fisher tests: E-box TCF4 peaks
is_tcf4 <- allgenes %in% gene_motif
is_query <- allgenes %in% NPC_up

fisher_motif_npc_up <- fisher.test(
  table(is_tcf4, is_query),
  alternative = "greater"
)

res_motif_npc_up <- data.frame(
  TCF4_Set = "TCF4_EBOX_peaks",
  KD_Set = "NPC_UP",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_motif_npc_up$estimate),
  P_Value = fisher_motif_npc_up$p.value
)

is_query <- allgenes %in% NPC_down

fisher_motif_npc_down <- fisher.test(
  table(is_tcf4, is_query),
  alternative = "greater"
)

res_motif_npc_down <- data.frame(
  TCF4_Set = "TCF4_EBOX_peaks",
  KD_Set = "NPC_DOWN",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_motif_npc_down$estimate),
  P_Value = fisher_motif_npc_down$p.value
)

is_query <- allgenes %in% Glut_up

fisher_motif_glut_up <- fisher.test(
  table(is_tcf4, is_query),
  alternative = "greater"
)

res_motif_glut_up <- data.frame(
  TCF4_Set = "TCF4_EBOX_peaks",
  KD_Set = "Glut_UP",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_motif_glut_up$estimate),
  P_Value = fisher_motif_glut_up$p.value
)

is_query <- allgenes %in% Glut_down

fisher_motif_glut_down <- fisher.test(
  table(is_tcf4, is_query),
  alternative = "greater"
)

res_motif_glut_down <- data.frame(
  TCF4_Set = "TCF4_EBOX_peaks",
  KD_Set = "Glut_DOWN",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_motif_glut_down$estimate),
  P_Value = fisher_motif_glut_down$p.value
)


# Fisher tests: promoter E-box TCF4 peaks
is_tcf4 <- allgenes %in% promoter_ebox_genes
is_query <- allgenes %in% NPC_up

fisher_promoter_npc_up <- fisher.test(
  table(is_tcf4, is_query),
  alternative = "greater"
)

res_promoter_npc_up <- data.frame(
  TCF4_Set = "Promoter_EBOX_TCF4_peaks",
  KD_Set = "NPC_UP",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_promoter_npc_up$estimate),
  P_Value = fisher_promoter_npc_up$p.value
)

is_query <- allgenes %in% NPC_down

fisher_promoter_npc_down <- fisher.test(
  table(is_tcf4, is_query),
  alternative = "greater"
)

res_promoter_npc_down <- data.frame(
  TCF4_Set = "Promoter_EBOX_TCF4_peaks",
  KD_Set = "NPC_DOWN",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_promoter_npc_down$estimate),
  P_Value = fisher_promoter_npc_down$p.value
)

is_query <- allgenes %in% Glut_up

fisher_promoter_glut_up <- fisher.test(
  table(is_tcf4, is_query),
  alternative = "greater"
)

res_promoter_glut_up <- data.frame(
  TCF4_Set = "Promoter_EBOX_TCF4_peaks",
  KD_Set = "Glut_UP",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_promoter_glut_up$estimate),
  P_Value = fisher_promoter_glut_up$p.value
)

is_query <- allgenes %in% Glut_down

fisher_promoter_glut_down <- fisher.test(
  table(is_tcf4, is_query),
  alternative = "greater"
)

res_promoter_glut_down <- data.frame(
  TCF4_Set = "Promoter_EBOX_TCF4_peaks",
  KD_Set = "Glut_DOWN",
  TCF4_Genes = sum(is_tcf4),
  KD_Genes = sum(is_query),
  Overlap = sum(is_tcf4 & is_query),
  Odds_Ratio = as.numeric(fisher_promoter_glut_down$estimate),
  P_Value = fisher_promoter_glut_down$p.value
)


# 18. Combine Fisher results
fisher_results <- bind_rows(
  res_all_npc_up,
  res_all_npc_down,
  res_all_glut_up,
  res_all_glut_down,
  res_motif_npc_up,
  res_motif_npc_down,
  res_motif_glut_up,
  res_motif_glut_down,
  res_promoter_npc_up,
  res_promoter_npc_down,
  res_promoter_glut_up,
  res_promoter_glut_down
)

fisher_results <- fisher_results %>%
  mutate(
    FDR = p.adjust(P_Value, method = "BH")
  ) %>%
  arrange(FDR)


fisher_results


# Create overlap gene lists for pathway analysis
overlap_all_npc_up <- intersect(gene_all, NPC_up)
overlap_all_npc_down <- intersect(gene_all, NPC_down)
overlap_all_glut_up <- intersect(gene_all, Glut_up)
overlap_all_glut_down <- intersect(gene_all, Glut_down)

overlap_motif_npc_up <- intersect(gene_motif, NPC_up)
overlap_motif_npc_down <- intersect(gene_motif, NPC_down)

overlap_promoter_npc_up <- intersect(promoter_ebox_genes, NPC_up)
overlap_promoter_npc_down <- intersect(promoter_ebox_genes, NPC_down)