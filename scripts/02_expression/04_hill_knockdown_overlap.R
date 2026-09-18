# Compare our TCF4 NPC genes with Hill's TCF4 KD genes
# Paper Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC5403663/#sec7
library(dplyr)
library(readr)

# Load Hill TCF4 KD data-
hill <- read_csv("NPC_KD_TCF4_Hill_2017.csv")

colnames(hill)
head(hill)

# Load our NPC gene files
npc_genes <- read_csv("All_TCF4_peak_genes.csv")
ebox_genes <- read_csv("EBOX_peak_genes.csv")
promoter_ebox_genes <- read_csv("Promoter_EBOX_genes.csv")

# Load background gene universe
background <- read_csv("final_ref_gene_2026_filtered.csv")

background_vec <- background %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  pull(gene)

length(background_vec)

# Check column names
colnames(npc_genes)
colnames(ebox_genes)
colnames(promoter_ebox_genes)

# Separate the gene column
npc_gene_vec <- npc_genes %>%
  filter(!is.na(Gene)) %>%
  distinct(Gene) %>%
  pull(Gene)

ebox_vec <- ebox_genes %>%
  filter(!is.na(Gene)) %>%
  distinct(Gene) %>%
  pull(Gene)

promoter_ebox_vec <- promoter_ebox_genes %>%
  filter(!is.na(Gene)) %>%
  distinct(Gene) %>%
  pull(Gene)

# Define significant Hill knockdown genes
hill_sig <- hill %>%
  filter(
    P_siRNA_1 < 0.05,
    P_siRNA_2 < 0.05,
    !is.na(Gene_Symbol)
  ) %>%
  filter(
    (Fold_change_siRNA_1 > 1 & Fold_change_siRNA_2 > 1) |
      (Fold_change_siRNA_1 < 1 & Fold_change_siRNA_2 < 1)
  ) %>%
  distinct(Gene_Symbol)

hill_sig_vec <- unique(hill_sig$Gene_Symbol)

length(hill_sig_vec)

# Separate upregulated and downregulated KD genes
hill_up <- hill %>%
  filter(
    P_siRNA_1 < 0.05,
    P_siRNA_2 < 0.05,
    Fold_change_siRNA_1 > 1,
    Fold_change_siRNA_2 > 1,
    !is.na(Gene_Symbol)
  ) %>%
  distinct(Gene_Symbol)

hill_down <- hill %>%
  filter(
    P_siRNA_1 < 0.05,
    P_siRNA_2 < 0.05,
    Fold_change_siRNA_1 < 1,
    Fold_change_siRNA_2 < 1,
    !is.na(Gene_Symbol)
  ) %>%
  distinct(Gene_Symbol)

hill_up_vec <- unique(hill_up$Gene_Symbol)
hill_down_vec <- unique(hill_down$Gene_Symbol)

length(hill_up_vec)
length(hill_down_vec)

# Fisher exact test function
fisher_overlap <- function(set1, set2, background){

  overlap <- length(intersect(set1, set2))

  set1_only <- length(setdiff(set1, set2))

  set2_only <- length(setdiff(set2, set1))

  neither <- length(background) -
    overlap -
    set1_only -
    set2_only

  fisher_matrix <- matrix(
    c(
      overlap,
      set1_only,
      set2_only,
      neither
    ),
    nrow = 2
  )

  fisher_result <- fisher.test(
    fisher_matrix,
    alternative = "greater"
  )

  return(
    data.frame(
      Overlap = overlap,
      Set1_Size = length(set1),
      Set2_Size = length(set2),
      Background_Size = length(background),
      Odds_Ratio = fisher_result$estimate,
      P_value = fisher_result$p.value
    )
  )
}

# Run fisher test
# NPC all genes vs Hill KD
npc_vs_hill <- fisher_overlap(
  npc_gene_vec,
  hill_sig_vec,
  background_vec
)

# NPC EBOX genes vs Hill KD
ebox_vs_hill <- fisher_overlap(
  ebox_vec,
  hill_sig_vec,
  background_vec
)

# NPC promoter EBOX genes vs Hill KD
promoter_vs_hill <- fisher_overlap(
  promoter_ebox_vec,
  hill_sig_vec,
  background_vec
)

# Do Upregulated and downregulated separately
npc_vs_hill_up <- fisher_overlap(
  npc_gene_vec,
  hill_up_vec,
  background_vec
)

npc_vs_hill_down <- fisher_overlap(
  npc_gene_vec,
  hill_down_vec,
  background_vec
)

# Combine results into a table
hill_results <- rbind(

  cbind(
    Comparison = "NPC all genes vs Hill KD",
    npc_vs_hill
  ),

  cbind(
    Comparison = "NPC EBOX genes vs Hill KD",
    ebox_vs_hill
  ),

  cbind(
    Comparison = "NPC promoter EBOX genes vs Hill KD",
    promoter_vs_hill
  ),

  cbind(
    Comparison = "NPC genes vs Hill UP genes",
    npc_vs_hill_up
  ),

  cbind(
    Comparison = "NPC genes vs Hill DOWN genes",
    npc_vs_hill_down
  )

)

hill_results
