# Since Blake’s dataset contains different genes from the Gandal and other datasets, I now want to compare the datasets to
# evaluate their overlap and determine whether the shared genes show the same direction of regulation across datasets
install.packages("writexl")
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(readxl)
library(ggplot2)
library(tidydr)
library(writexl)

# Import background genes across some datsets
blake_background_gene <- read.csv("Data integration/Blake final dataset/Blake_GSE48367_limma_all_results.csv")
colnames(blake_background_gene)

Gandal_background_gene <- read_excel("Transcriptome_SCZ_AD_BPD.xlsx", sheet = "DGE")
colnames(Gandal_background_gene)

# Find the overlaps between these two dataset and if they are in a same direction or not
# Clean Blake genes
blake_ref <- blake_background_gene %>%
  select(
    gene,
    Blake_logFC = logFC,
    Blake_FDR = adj.P.Val
  ) %>%
  mutate(
    gene = str_trim(gene),

    Blake_direction = case_when(
      Blake_logFC > 0 ~ "Up",
      Blake_logFC < 0 ~ "Down",
      TRUE ~ "No_change"
    )
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    gene != "-"
  ) %>%
  distinct(gene, .keep_all = TRUE)

# Clean Gandal SCZ genes
gandal_SCZ <- Gandal_background_gene %>%
  select(
    gene_name,
    SCZ.log2FC,
    SCZ.fdr
  ) %>%
  rename(gene = gene_name) %>%
  mutate(
    gene = str_trim(gene),

    SCZ_direction = case_when(
      SCZ.log2FC > 0 ~ "Up",
      SCZ.log2FC < 0 ~ "Down",
      TRUE ~ "No_change"
    )
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    gene != "-"
  ) %>%
  distinct(gene, .keep_all = TRUE)

# Find overlapping genes
overlap_SCZ <- inner_join(
  blake_ref,
  gandal_SCZ,
  by = "gene"
)

# Number of overlapping genes
nrow(overlap_SCZ)

head(overlap_SCZ)

# Compare direction between Blake and SCZ
overlap_SCZ <- overlap_SCZ %>%
  mutate(
    Same_direction = Blake_direction == SCZ_direction
  )

# Summary table
direction_summary_SCZ <- data.frame(

  Total_overlap = nrow(overlap_SCZ),

  Same_direction = sum(
    overlap_SCZ$Same_direction,
    na.rm = TRUE
  ),

  Opposite_direction = sum(
    !overlap_SCZ$Same_direction,
    na.rm = TRUE
  )
)

direction_summary_SCZ


#### Let's take a look to significant genes in both dataset and see if there is overlap
# Significant Blake genes only
blake_sig <- blake_background_gene %>%
  filter(adj.P.Val < 0.05) %>%
  select(
    gene,
    Blake_logFC = logFC,
    Blake_FDR = adj.P.Val
  ) %>%
  mutate(
    gene = str_trim(gene),

    Blake_direction = case_when(
      Blake_logFC > 0 ~ "Up",
      Blake_logFC < 0 ~ "Down",
      TRUE ~ "No_change"
    )
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    gene != "-"
  ) %>%
  distinct(gene, .keep_all = TRUE)


# Significant Gandal SCZ genes only
gandal_SCZ_sig <- Gandal_background_gene %>%
  filter(SCZ.fdr < 0.05) %>%
  select(
    gene_name,
    SCZ.log2FC,
    SCZ.fdr
  ) %>%
  rename(gene = gene_name) %>%
  mutate(
    gene = str_trim(gene),

    SCZ_direction = case_when(
      SCZ.log2FC > 0 ~ "Up",
      SCZ.log2FC < 0 ~ "Down",
      TRUE ~ "No_change"
    )
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    gene != "-"
  ) %>%
  distinct(gene, .keep_all = TRUE)


# Find overlapping significant genes
overlap_SCZ_sig <- inner_join(
  blake_sig,
  gandal_SCZ_sig,
  by = "gene"
)

# Number of overlapping significant genes
nrow(overlap_SCZ_sig)

head(overlap_SCZ_sig)


# Compare directions
overlap_SCZ_sig <- overlap_SCZ_sig %>%
  mutate(
    Same_direction = Blake_direction == SCZ_direction
  )


# Summary table
direction_summary_SCZ_sig <- data.frame(

  Total_overlap = nrow(overlap_SCZ_sig),

  Same_direction = sum(
    overlap_SCZ_sig$Same_direction,
    na.rm = TRUE
  ),

  Opposite_direction = sum(
    !overlap_SCZ_sig$Same_direction,
    na.rm = TRUE
  )
)

direction_summary_SCZ_sig

# Extract genes with opposite directions
opposite_direction_genes <- overlap_SCZ_sig %>%
  filter(Same_direction == FALSE)

write.csv(opposite_direction_genes, "opposite_blake&Gandal_genes")
# View first rows
View(opposite_direction_genes)
# Number of opposite-direction genes
nrow(opposite_direction_genes)


# ============================================================
# Compare Fromer SCZ genes with Blake and Gandal SCZ
# ============================================================
library(dplyr)
library(stringr)

# Import data
fromgenes = readLines("Data integration/Fromer_genes.txt")
fromgenesup = readLines("Data integration/Fromer-up.txt")
fromgenesdown = readLines("Data integration/Fromer-down.txt")

# Clean Fromer gene lists
fromer_all <- fromgenes %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

fromer_up <- fromgenesup %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

fromer_down <- fromgenesdown %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()


# Create Fromer dataframe with direction
fromer_ref <- data.frame(gene = fromer_all) %>%
  mutate(
    Fromer_direction = case_when(
      gene %in% fromer_up ~ "Up",
      gene %in% fromer_down ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Fromer_direction))


# ============================================================
# Blake significant DEG dataset
# ============================================================

blake_sig <- blake_results_clean %>%
  filter(adj.P.Val < 0.05) %>%
  select(
    gene,
    Blake_logFC = logFC,
    Blake_FDR = adj.P.Val
  ) %>%
  mutate(
    gene = str_trim(gene),
    Blake_direction = case_when(
      Blake_logFC > 0 ~ "Up",
      Blake_logFC < 0 ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(gene), gene != "", gene != "-") %>%
  distinct(gene, .keep_all = TRUE)


# ============================================================
# Gandal significant SCZ DEG dataset
# ============================================================

gandal_SCZ_sig <- Gandal_background_gene %>%
  filter(SCZ.fdr < 0.05) %>%
  select(
    gene = gene_name,
    Gandal_logFC = SCZ.log2FC,
    Gandal_FDR = SCZ.fdr
  ) %>%
  mutate(
    gene = str_trim(gene),
    Gandal_direction = case_when(
      Gandal_logFC > 0 ~ "Up",
      Gandal_logFC < 0 ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(gene), gene != "", gene != "-") %>%
  distinct(gene, .keep_all = TRUE)


# ============================================================
# Overlap: Fromer vs Blake
# ============================================================

Fromer_vs_Blake <- fromer_ref %>%
  inner_join(blake_sig, by = "gene") %>%
  mutate(
    Fromer_Blake_same_direction = Fromer_direction == Blake_direction
  )

nrow(Fromer_vs_Blake)

Fromer_Blake_direction_summary <- data.frame(
  Comparison = "Fromer vs Blake",
  Total_overlap = nrow(Fromer_vs_Blake),
  Same_direction = sum(Fromer_vs_Blake$Fromer_Blake_same_direction, na.rm = TRUE),
  Opposite_direction = sum(!Fromer_vs_Blake$Fromer_Blake_same_direction, na.rm = TRUE)
)

Fromer_Blake_direction_summary


# ============================================================
# Overlap: Fromer vs Gandal SCZ
# ============================================================

Fromer_vs_Gandal <- fromer_ref %>%
  inner_join(gandal_SCZ_sig, by = "gene") %>%
  mutate(
    Fromer_Gandal_same_direction = Fromer_direction == Gandal_direction
  )

nrow(Fromer_vs_Gandal)

Fromer_Gandal_direction_summary <- data.frame(
  Comparison = "Fromer vs Gandal SCZ",
  Total_overlap = nrow(Fromer_vs_Gandal),
  Same_direction = sum(Fromer_vs_Gandal$Fromer_Gandal_same_direction, na.rm = TRUE),
  Opposite_direction = sum(!Fromer_vs_Gandal$Fromer_Gandal_same_direction, na.rm = TRUE)
)

Fromer_Gandal_direction_summary


# ============================================================
# Overlap across all three: Fromer + Blake + Gandal SCZ
# ============================================================

Fromer_Blake_Gandal <- fromer_ref %>%
  inner_join(blake_sig, by = "gene") %>%
  inner_join(gandal_SCZ_sig, by = "gene") %>%
  mutate(
    Blake_Gandal_same_direction = Blake_direction == Gandal_direction,
    Fromer_Blake_same_direction = Fromer_direction == Blake_direction,
    Fromer_Gandal_same_direction = Fromer_direction == Gandal_direction,
    All_same_direction = Fromer_direction == Blake_direction &
      Fromer_direction == Gandal_direction
  )

nrow(Fromer_Blake_Gandal)


# Summary across all three
Fromer_Blake_Gandal_summary <- data.frame(
  Comparison = c(
    "Fromer vs Blake",
    "Fromer vs Gandal SCZ",
    "Blake vs Gandal SCZ",
    "All three same direction"
  ),
  Total_overlap = c(
    nrow(Fromer_Blake_Gandal),
    nrow(Fromer_Blake_Gandal),
    nrow(Fromer_Blake_Gandal),
    nrow(Fromer_Blake_Gandal)
  ),
  Same_direction = c(
    sum(Fromer_Blake_Gandal$Fromer_Blake_same_direction, na.rm = TRUE),
    sum(Fromer_Blake_Gandal$Fromer_Gandal_same_direction, na.rm = TRUE),
    sum(Fromer_Blake_Gandal$Blake_Gandal_same_direction, na.rm = TRUE),
    sum(Fromer_Blake_Gandal$All_same_direction, na.rm = TRUE)
  ),
  Opposite_or_not_all_same = c(
    sum(!Fromer_Blake_Gandal$Fromer_Blake_same_direction, na.rm = TRUE),
    sum(!Fromer_Blake_Gandal$Fromer_Gandal_same_direction, na.rm = TRUE),
    sum(!Fromer_Blake_Gandal$Blake_Gandal_same_direction, na.rm = TRUE),
    sum(!Fromer_Blake_Gandal$All_same_direction, na.rm = TRUE)
  )
)

Fromer_Blake_Gandal_summary

# Fromer vs Blake opposite genes
Fromer_Blake_opposite <- Fromer_vs_Blake %>%
  filter(Fromer_Blake_same_direction == FALSE)

# Fromer vs Gandal opposite genes
Fromer_Gandal_opposite <- Fromer_vs_Gandal %>%
  filter(Fromer_Gandal_same_direction == FALSE)

# Genes shared by all three but not all in same direction
Fromer_Blake_Gandal_not_all_same <- Fromer_Blake_Gandal %>%
  filter(All_same_direction == FALSE)

View(Fromer_Blake_opposite)
View(Fromer_Gandal_opposite)
View(Fromer_Blake_Gandal_not_all_same)

####################################################################
# Now let's compare the overlaps

overlap_Gandal <- read.csv("SCZ_TCF4_overlap_Ganda.csv")
overlap_Blake <- read.csv("Blake_TCF4_overlap_long_with_direction.csv")
overlap_Fromer <- read.csv("TCF4_Fromer_overlap.csv")

head(overlap_Gandal)
head(overlap_Blake)
head(overlap_Fromer)

# Clean each overlap dataset
Gandal_clean <- overlap_Gandal %>%
  select(
    gene,
    TCF4_Set,
    Gandal_direction = SCZ_Direction,
    SCZ.log2FC,
    SCZ.fdr
  ) %>%
  distinct(gene, TCF4_Set, .keep_all = TRUE)

Blake_clean <- overlap_Blake %>%
  select(
    gene,
    TCF4_Set,
    Blake_direction = Blake_Direction,
    Blake_logFC = Blake_logFC.x,
    Blake_FDR = Blake_FDR.x
  ) %>%
  distinct(gene, TCF4_Set, .keep_all = TRUE)

Fromer_clean <- overlap_Fromer %>%
  select(
    gene,
    TCF4_Set,
    Fromer_direction = Fromer_Direction
  ) %>%
  distinct(gene, TCF4_Set, .keep_all = TRUE)


# Pairwise overlaps and direction agreement
# Gandal vs Blake
Gandal_Blake_overlap <- inner_join(
  Gandal_clean,
  Blake_clean,
  by = c("gene", "TCF4_Set")
) %>%
  mutate(
    Same_direction = Gandal_direction == Blake_direction
  )

# Gandal vs Fromer
Gandal_Fromer_overlap <- inner_join(
  Gandal_clean,
  Fromer_clean,
  by = c("gene", "TCF4_Set")
) %>%
  mutate(
    Same_direction = Gandal_direction == Fromer_direction
  )

# Blake vs Fromer
Blake_Fromer_overlap <- inner_join(
  Blake_clean,
  Fromer_clean,
  by = c("gene", "TCF4_Set")
) %>%
  mutate(
    Same_direction = Blake_direction == Fromer_direction
  )

# Three-way overlap and direction agreement
Gandal_Blake_Fromer_overlap <- Gandal_clean %>%
  inner_join(Blake_clean, by = c("gene", "TCF4_Set")) %>%
  inner_join(Fromer_clean, by = c("gene", "TCF4_Set")) %>%
  mutate(
    Gandal_Blake_same = Gandal_direction == Blake_direction,
    Gandal_Fromer_same = Gandal_direction == Fromer_direction,
    Blake_Fromer_same = Blake_direction == Fromer_direction,
    All_three_same = Gandal_direction == Blake_direction &
      Gandal_direction == Fromer_direction
  )

# Clean table showing direction for every shared gene
Gandal_Blake_Fromer_direction_table <- Gandal_Blake_Fromer_overlap %>%
  select(
    gene,
    TCF4_Set,
    Gandal_direction,
    Blake_direction,
    Fromer_direction,
    Gandal_Blake_same,
    Gandal_Fromer_same,
    Blake_Fromer_same,
    All_three_same,
    SCZ.log2FC,
    SCZ.fdr,
    Blake_logFC,
    Blake_FDR
  ) %>%
  arrange(TCF4_Set, gene)

View(Gandal_Blake_Fromer_direction_table)

# 5. Summary tables
pairwise_direction_summary <- data.frame(
  Comparison = c(
    "Gandal vs Blake",
    "Gandal vs Fromer",
    "Blake vs Fromer"
  ),
  Total_overlap = c(
    nrow(Gandal_Blake_overlap),
    nrow(Gandal_Fromer_overlap),
    nrow(Blake_Fromer_overlap)
  ),
  Same_direction = c(
    sum(Gandal_Blake_overlap$Same_direction, na.rm = TRUE),
    sum(Gandal_Fromer_overlap$Same_direction, na.rm = TRUE),
    sum(Blake_Fromer_overlap$Same_direction, na.rm = TRUE)
  ),
  Opposite_direction = c(
    sum(!Gandal_Blake_overlap$Same_direction, na.rm = TRUE),
    sum(!Gandal_Fromer_overlap$Same_direction, na.rm = TRUE),
    sum(!Blake_Fromer_overlap$Same_direction, na.rm = TRUE)
  )
)

pairwise_direction_summary

# Genes with the SAME direction between Gandal and Blake
Gandal_Blake_same_direction_genes <- Gandal_Blake_overlap %>%
  filter(Same_direction == TRUE) %>%
  select(
    gene,
    TCF4_Set,
    Gandal_direction,
    Blake_direction,
    SCZ.log2FC,
    SCZ.fdr,
    Blake_logFC,
    Blake_FDR
  ) %>%
  arrange(TCF4_Set, gene)

View(Gandal_Blake_same_direction_genes)
nrow(Gandal_Blake_same_direction_genes)


# Save all overlap/direction tables into one Excel file
write_xlsx(
  list(
    "Gandal_overlaps" = Gandal_clean,
    "Blake_overlaps" = Blake_clean,
    "Fromer_overlaps" = Fromer_clean,
    "Gandal_vs_Blake" = Gandal_Blake_overlap,
    "Gandal_vs_Fromer" = Gandal_Fromer_overlap,
    "Blake_vs_Fromer" = Blake_Fromer_overlap,
    "All_three_overlap" = Gandal_Blake_Fromer_overlap,
    "All_three_direction_table" = Gandal_Blake_Fromer_direction_table
  ),
  path = "TCF4_Gandal_Blake_Fromer_overlap_direction_tables.xlsx"
)



## Genes with OPPOSITE direction between Gandal and Blake
Gandal_Blake_opposite_direction_genes <- Gandal_Blake_overlap %>%
  filter(Same_direction == FALSE) %>%
  select(
    gene,
    TCF4_Set,
    Gandal_direction,
    Blake_direction,
    SCZ.log2FC,
    SCZ.fdr,
    Blake_logFC,
    Blake_FDR
  ) %>%
  arrange(TCF4_Set, gene)

View(Gandal_Blake_opposite_direction_genes)
nrow(Gandal_Blake_opposite_direction_genes)

Gandal_Blake_direction_summary <- data.frame(
  Total_overlap = nrow(Gandal_Blake_overlap),
  Same_direction = nrow(Gandal_Blake_same_direction_genes),
  Opposite_direction = nrow(Gandal_Blake_opposite_direction_genes)
)

Gandal_Blake_direction_summary
