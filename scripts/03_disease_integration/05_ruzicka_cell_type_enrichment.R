# TCF4 genes vs Ruzika's gene (different cell types)
#https://pmc.ncbi.nlm.nih.gov/articles/PMC12772489/
# Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(tidydr)
library(stringr)
library(readr)
library(tibble)

# Import data from Ruzicka's dataset
file_path <- "Data integration/Ruzika_CellTypes.xlsx"
sheet_names <- excel_sheets(file_path)

# Import all tabs
data_list <- list()
for (s in sheet_names) {
  data_list[[s]] <- read_excel(file_path, sheet = s)
}

# Check tabs
names(data_list)
length(data_list)
head(data_list[[1]])

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

# For each cell type, we need a background gene list (reference) including all genes of that sheet + Sig genes of that sheet
# Make backround gene list for all 25 brain derived-cell types
background_list <- list()

for (cell_type in names(data_list)) {
  df <- data_list[[cell_type]]
  background_genes <- unique(df$gene)
  background_genes <- background_genes[!is.na(background_genes) & background_genes != ""]  ## Clean the gene list
  background_list[[cell_type]] <- background_genes  # Store the cleaned background gene list in the list
}

# Let's check if each cell type has it owns bckground
length(background_list[[1]])
length(background_list[[2]])
length(background_list[[3]])

sapply(background_list, length)

# Make a significant gene list (up and down) from each cell type gene list
up_list <- list()
down_list <- list()
no_direction_list <- list()

for (cell_type in names(data_list)) {
  df <- data_list[[cell_type]]
  no_direction_genes <- df$gene[df$Meta_adj.P.Val < 0.05]
  up_genes <- df$gene[df$Meta_adj.P.Val < 0.05 & df$Meta_logFC > 0]
  down_genes <- df$gene[df$Meta_adj.P.Val < 0.05 & df$Meta_logFC < 0]
  no_direction_genes <- unique(no_direction_genes[!is.na(no_direction_genes) & no_direction_genes !=""])
  up_genes <- unique(up_genes[!is.na(up_genes) & up_genes != ""])
  down_genes <- unique(down_genes[!is.na(down_genes) & down_genes != ""])
  no_direction_list [[cell_type]] <- no_direction_genes
  up_list[[cell_type]] <- up_genes
  down_list[[cell_type]] <- down_genes
}

# Check if the loop worked
sapply(no_direction_list, length)
sapply(up_list, length)
sapply(down_list, length)

# Create logical vectors for our TCF4 gene sets
# For each gene in each background genes, check if it’s in our TCF4 list or not
# We should do this inside a loop, because each cell type has its own background
tcf4_flags <- list()

for (cell_type in names(background_list)) {
  background_genes <- background_list[[cell_type]]
  is_gene_all <-  background_genes %in% gene_all
  is_gene_motif <- background_genes %in% gene_motif
  is_promoter_ebox <- background_genes %in% promoter_ebox_genes
  # Store results (save everything for one cell)
  tcf4_flags[[cell_type]] <- list(
    is_gene_all = is_gene_all,
    is_gene_motif = is_gene_motif,
    is_promoter_ebox = is_promoter_ebox
  )
}

# Check if the loop worked
tcf4_summary <- data.frame()

for (cell_type in names(tcf4_flags)) {
  tcf4_summary <- rbind(tcf4_summary, data.frame(
    Cell_Type = cell_type,
    TCF4_all = sum(tcf4_flags[[cell_type]]$is_gene_all),
    TCF4_motif = sum(tcf4_flags[[cell_type]]$is_gene_motif),
    TCF4_promoter = sum(tcf4_flags[[cell_type]]$is_promoter_ebox)
  ))
}

tcf4_summary

# Enrichment Analysis
# Step1: Create empty result table
enrichment_results <- data.frame()

# Step2: Loop through each cell type
for (cell_type in names(background_list)) {
  # Get the background gene universe for this cell type
  background_genes <- background_list[[cell_type]]

  # Get significant gene lists for this cell type
  sig_all <- no_direction_list[[cell_type]]   # all significant genes
  sig_up <- up_list[[cell_type]]              # upregulated genes
  sig_down <- down_list[[cell_type]]          # downregulated genes

  # Get TCF4 logical vectors
  # TRUE = gene belongs to TCF4 set and FALSE = gene does not belong to TCF4 set
  is_gene_all <- tcf4_flags[[cell_type]]$is_gene_all
  is_gene_motif <- tcf4_flags[[cell_type]]$is_gene_motif
  is_promoter_ebox <- tcf4_flags[[cell_type]]$is_promoter_ebox

  # Convert significant gene lists into logical vectors aligned with background (we should make everything align together for fisher test)
  # TRUE = gene is significant in this cell type
  # FALSE = gene is not significant
  is_sig_all <- background_genes %in% sig_all
  is_sig_up <- background_genes %in% sig_up
  is_sig_down <- background_genes %in% sig_down

  # Fisher test: Are TCF4_all genes enriched in significant genes?
  fisher_all_sig <- fisher.test(
    table(is_gene_all, is_sig_all),
    alternative = "greater"
  )

  # Fisher test: Are TCF4_motif genes enriched?
  fisher_motif_sig <- fisher.test(
    table(is_gene_motif, is_sig_all),
    alternative = "greater"
  )

  # Fisher test: Are TCF4_promoter_ebox genes enriched?
  fisher_promoter_sig <- fisher.test(
    table(is_promoter_ebox, is_sig_all),
    alternative = "greater"
  )
  # Save results for this cell type
    enrichment_results <- rbind(enrichment_results, data.frame(
    Cell_Type = cell_type,                                       # Cell type name
    Background_n = length(background_genes),                     # Total number of background genes
    Sig_all_n = length(sig_all),                                 # Number of significant genes
    TCF4_all_OR = as.numeric(fisher_all_sig$estimate),           # Odds ratio and p-value for TCF4_all
    TCF4_all_p = fisher_all_sig$p.value,
    TCF4_motif_OR = as.numeric(fisher_motif_sig$estimate),       # Odds ratio and p-value for TCF4_motif
    TCF4_motif_p = fisher_motif_sig$p.value,
    TCF4_promoter_OR = as.numeric(fisher_promoter_sig$estimate), # Odds ratio and p-value for TCF4_promoter
    TCF4_promoter_p = fisher_promoter_sig$p.value
  ))
}

# Step3: Adjust p-value (Apply BH correction to each [FDR])
enrichment_results$TCF4_all_FDR <- p.adjust(
  enrichment_results$TCF4_all_p,
  method = "BH"
)

enrichment_results$TCF4_motif_FDR <- p.adjust(
  enrichment_results$TCF4_motif_p,
  method = "BH"
)

enrichment_results$TCF4_promoter_FDR <- p.adjust(
  enrichment_results$TCF4_promoter_p,
  method = "BH"
)

View(enrichment_results)

# Get overlapped genes between TCF4 and significant genes in each cell type
# For TCF4_all
overlap_list <- list()

for (cell_type in names(background_list)) {
  sig_genes <- no_direction_list[[cell_type]]
  overlap <- intersect(sig_genes, gene_all)
  overlap_list[[cell_type]] <- overlap
}

# For motif and promoter separately
overlap_motif <- list()
overlap_promoter <- list()

for (cell_type in names(background_list)) {
  sig_genes <- no_direction_list[[cell_type]]
  overlap_motif[[cell_type]] <- intersect(sig_genes, gene_motif)
  overlap_promoter[[cell_type]] <- intersect(sig_genes, promoter_ebox_genes)
}

# Check how many genes overlap (overalp summary)
overlap_summary <- data.frame(
  Cell_Type = names(overlap_list),
  Overlap_all = sapply(overlap_list, length),
  Overlap_motif = sapply(overlap_motif, length),
  Overlap_promoter = sapply(overlap_promoter, length)
)

overlap_summary
View(overlap_list)
