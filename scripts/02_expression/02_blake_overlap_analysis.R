####################################################################
# Overlap Analysis: TCF4 Gene Sets vs Old Blake DEG Lists
####################################################################
library(dplyr)
library(stringr)

## Import TCF4 peak annotation files and blake's background gene

npc_annotation <- read.csv("npc_hg19_peak_annotation.csv")
mcclay_annotation <- read.csv("mcclay_hg19_peak_annotation.csv")
forrest_annotation <- read.csv("forrest_hg19_peak_annotation.csv")

## Extract clean TCF4 gene lists
## Change SYMBOL to the correct column name if needed

colnames(npc_annotation)
colnames(mcclay_annotation)
colnames(forrest_annotation)

npc_tcf4_genes <- npc_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA"] %>%
  unique()

mcclay_tcf4_genes <- mcclay_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA"] %>%
  unique()

forrest_tcf4_genes <- forrest_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA"] %>%
  unique()


## Check TCF4 gene counts
length(npc_tcf4_genes)
length(mcclay_tcf4_genes)
length(forrest_tcf4_genes)



# Blake Background and Old Blake Lists
# load Blake gene-level result
blake_gene_level_results <- read.csv(
  "Blake_GSE48367_limma_all_results_probe_level.csv",
  stringsAsFactors = FALSE
)

# Check column names
colnames(blake_gene_level_results)

# Create correct Blake background
blake_background_genes <- blake_gene_level_results$gene %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

length(blake_background_genes)
head(blake_background_genes)

####################################################################
# Import Old Blake DEG Lists
####################################################################
old_blake_all <- readLines("Data integration/TCF4 KD - 1205 DEG.txt") %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA"] %>%
  unique()

old_blake_up <- readLines("Data integration/TCF4 blake up n470.txt") %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA"] %>%
  unique()

old_blake_down <- readLines("Data integration/TCF4 blake down n665.txt") %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA"] %>%
  unique()


## Check old Blake gene counts
length(old_blake_all)
length(old_blake_up)
length(old_blake_down)

####################################################################
# Calculate Simple Gene Overlaps
####################################################################
## NPC TCF4 overlaps
npc_old_blake_all_overlap <- intersect(npc_tcf4_genes, old_blake_all)
npc_old_blake_up_overlap <- intersect(npc_tcf4_genes, old_blake_up)
npc_old_blake_down_overlap <- intersect(npc_tcf4_genes, old_blake_down)


## McClay TCF4 overlaps
mcclay_old_blake_all_overlap <- intersect(mcclay_tcf4_genes, old_blake_all)
mcclay_old_blake_up_overlap <- intersect(mcclay_tcf4_genes, old_blake_up)
mcclay_old_blake_down_overlap <- intersect(mcclay_tcf4_genes, old_blake_down)


## Forrest TCF4 overlaps
forrest_old_blake_all_overlap <- intersect(forrest_tcf4_genes, old_blake_all)
forrest_old_blake_up_overlap <- intersect(forrest_tcf4_genes, old_blake_up)
forrest_old_blake_down_overlap <- intersect(forrest_tcf4_genes, old_blake_down)

####################################################################
# Make Simple Overlap Summary Table
####################################################################
old_blake_overlap_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC_TCF4", "NPC_TCF4", "NPC_TCF4",
    "McClay_TCF4", "McClay_TCF4", "McClay_TCF4",
    "Forrest_TCF4", "Forrest_TCF4", "Forrest_TCF4"
  ),

  Blake_Set = c(
    "Old_Blake_all_DEGs",
    "Old_Blake_upregulated",
    "Old_Blake_downregulated",
    "Old_Blake_all_DEGs",
    "Old_Blake_upregulated",
    "Old_Blake_downregulated",
    "Old_Blake_all_DEGs",
    "Old_Blake_upregulated",
    "Old_Blake_downregulated"
  ),

  TCF4_Gene_Count = c(
    length(npc_tcf4_genes),
    length(npc_tcf4_genes),
    length(npc_tcf4_genes),
    length(mcclay_tcf4_genes),
    length(mcclay_tcf4_genes),
    length(mcclay_tcf4_genes),
    length(forrest_tcf4_genes),
    length(forrest_tcf4_genes),
    length(forrest_tcf4_genes)
  ),

  Blake_Gene_Count = c(
    length(old_blake_all),
    length(old_blake_up),
    length(old_blake_down),
    length(old_blake_all),
    length(old_blake_up),
    length(old_blake_down),
    length(old_blake_all),
    length(old_blake_up),
    length(old_blake_down)
  ),

  Overlap_Genes = c(
    length(npc_old_blake_all_overlap),
    length(npc_old_blake_up_overlap),
    length(npc_old_blake_down_overlap),
    length(mcclay_old_blake_all_overlap),
    length(mcclay_old_blake_up_overlap),
    length(mcclay_old_blake_down_overlap),
    length(forrest_old_blake_all_overlap),
    length(forrest_old_blake_up_overlap),
    length(forrest_old_blake_down_overlap)
  )
)

old_blake_overlap_summary <- old_blake_overlap_summary %>%
  mutate(
    Percent_of_TCF4_Genes = round(
      Overlap_Genes / TCF4_Gene_Count * 100,
      3
    ),
    Percent_of_Blake_Genes = round(
      Overlap_Genes / Blake_Gene_Count * 100,
      3
    )
  )

old_blake_overlap_summary

####################################################################
# Simple Fisher Exact Tests: TCF4 Gene Sets vs Old Blake DEG Lists
####################################################################
## Restrict TCF4 gene sets to Blake background
npc_in_background <- intersect(npc_tcf4_genes, blake_background_genes)
mcclay_in_background <- intersect(mcclay_tcf4_genes, blake_background_genes)
forrest_in_background <- intersect(forrest_tcf4_genes, blake_background_genes)


## Restrict old Blake DEG lists to Blake background
old_blake_all_in_background <- intersect(old_blake_all, blake_background_genes)
old_blake_up_in_background <- intersect(old_blake_up, blake_background_genes)
old_blake_down_in_background <- intersect(old_blake_down, blake_background_genes)

length(npc_in_background)
length(old_blake_all_in_background)
length(blake_background_genes)

####################################################################
# NPC TCF4 vs Old Blake
####################################################################
## NPC vs all old Blake DEGs
npc_all_table <- table(
  blake_background_genes %in% npc_in_background,
  blake_background_genes %in% old_blake_all_in_background
)

npc_all_fisher <- fisher.test(
  npc_all_table,
  alternative = "greater"
)


## NPC vs old Blake upregulated genes
npc_up_table <- table(
  blake_background_genes %in% npc_in_background,
  blake_background_genes %in% old_blake_up_in_background
)

npc_up_fisher <- fisher.test(
  npc_up_table,
  alternative = "greater"
)


## NPC vs old Blake downregulated genes
npc_down_table <- table(
  blake_background_genes %in% npc_in_background,
  blake_background_genes %in% old_blake_down_in_background
)

npc_down_fisher <- fisher.test(
  npc_down_table,
  alternative = "greater"
)

####################################################################
# McClay TCF4 vs Old Blake
####################################################################

## McClay vs all old Blake DEGs
mcclay_all_table <- table(
  blake_background_genes %in% mcclay_in_background,
  blake_background_genes %in% old_blake_all_in_background
)

mcclay_all_fisher <- fisher.test(
  mcclay_all_table,
  alternative = "greater"
)


## McClay vs old Blake upregulated genes
mcclay_up_table <- table(
  blake_background_genes %in% mcclay_in_background,
  blake_background_genes %in% old_blake_up_in_background
)

mcclay_up_fisher <- fisher.test(
  mcclay_up_table,
  alternative = "greater"
)


## McClay vs old Blake downregulated genes
mcclay_down_table <- table(
  blake_background_genes %in% mcclay_in_background,
  blake_background_genes %in% old_blake_down_in_background
)

mcclay_down_fisher <- fisher.test(
  mcclay_down_table,
  alternative = "greater"
)

####################################################################
# Forrest TCF4 vs Old Blake
####################################################################

## Forrest vs all old Blake DEGs
forrest_all_table <- table(
  blake_background_genes %in% forrest_in_background,
  blake_background_genes %in% old_blake_all_in_background
)

forrest_all_fisher <- fisher.test(
  forrest_all_table,
  alternative = "greater"
)


## Forrest vs old Blake upregulated genes
forrest_up_table <- table(
  blake_background_genes %in% forrest_in_background,
  blake_background_genes %in% old_blake_up_in_background
)

forrest_up_fisher <- fisher.test(
  forrest_up_table,
  alternative = "greater"
)


## Forrest vs old Blake downregulated genes
forrest_down_table <- table(
  blake_background_genes %in% forrest_in_background,
  blake_background_genes %in% old_blake_down_in_background
)

forrest_down_fisher <- fisher.test(
  forrest_down_table,
  alternative = "greater"
)


####################################################################
# Create Simple Fisher Summary Table
####################################################################
old_blake_fisher_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC_TCF4", "NPC_TCF4", "NPC_TCF4",
    "McClay_TCF4", "McClay_TCF4", "McClay_TCF4",
    "Forrest_TCF4", "Forrest_TCF4", "Forrest_TCF4"
  ),

  Blake_Set = c(
    "Old_Blake_all_DEGs",
    "Old_Blake_upregulated",
    "Old_Blake_downregulated",
    "Old_Blake_all_DEGs",
    "Old_Blake_upregulated",
    "Old_Blake_downregulated",
    "Old_Blake_all_DEGs",
    "Old_Blake_upregulated",
    "Old_Blake_downregulated"
  ),

  Background_Genes = rep(length(blake_background_genes), 9),

  TCF4_Genes_in_Background = c(
    length(npc_in_background),
    length(npc_in_background),
    length(npc_in_background),
    length(mcclay_in_background),
    length(mcclay_in_background),
    length(mcclay_in_background),
    length(forrest_in_background),
    length(forrest_in_background),
    length(forrest_in_background)
  ),

  Blake_Genes_in_Background = c(
    length(old_blake_all_in_background),
    length(old_blake_up_in_background),
    length(old_blake_down_in_background),
    length(old_blake_all_in_background),
    length(old_blake_up_in_background),
    length(old_blake_down_in_background),
    length(old_blake_all_in_background),
    length(old_blake_up_in_background),
    length(old_blake_down_in_background)
  ),

  Overlap_Genes = c(
    length(intersect(npc_in_background, old_blake_all_in_background)),
    length(intersect(npc_in_background, old_blake_up_in_background)),
    length(intersect(npc_in_background, old_blake_down_in_background)),
    length(intersect(mcclay_in_background, old_blake_all_in_background)),
    length(intersect(mcclay_in_background, old_blake_up_in_background)),
    length(intersect(mcclay_in_background, old_blake_down_in_background)),
    length(intersect(forrest_in_background, old_blake_all_in_background)),
    length(intersect(forrest_in_background, old_blake_up_in_background)),
    length(intersect(forrest_in_background, old_blake_down_in_background))
  ),

  Odds_Ratio = c(
    npc_all_fisher$estimate,
    npc_up_fisher$estimate,
    npc_down_fisher$estimate,
    mcclay_all_fisher$estimate,
    mcclay_up_fisher$estimate,
    mcclay_down_fisher$estimate,
    forrest_all_fisher$estimate,
    forrest_up_fisher$estimate,
    forrest_down_fisher$estimate
  ),

  P_Value = c(
    npc_all_fisher$p.value,
    npc_up_fisher$p.value,
    npc_down_fisher$p.value,
    mcclay_all_fisher$p.value,
    mcclay_up_fisher$p.value,
    mcclay_down_fisher$p.value,
    forrest_all_fisher$p.value,
    forrest_up_fisher$p.value,
    forrest_down_fisher$p.value
  )
)

old_blake_fisher_summary <- old_blake_fisher_summary %>%
  mutate(
    Fisher_FDR = p.adjust(P_Value, method = "BH"),

    Percent_TCF4_Overlap = round(
      Overlap_Genes / TCF4_Genes_in_Background * 100,
      3
    ),

    Percent_Blake_Overlap = round(
      Overlap_Genes / Blake_Genes_in_Background * 100,
      3
    )
  )

old_blake_fisher_summary

####################################################################
# Simple Enrichment Analysis: New Blake DEGs vs TCF4 Gene Sets
####################################################################
# 1. Import TCF4 Annotation Files
npc_annotation <- read.csv("npc_hg19_peak_annotation.csv")
mcclay_annotation <- read.csv("mcclay_hg19_peak_annotation.csv")
forrest_annotation <- read.csv("forrest_hg19_peak_annotation.csv")

## Check gene-symbol column names
colnames(npc_annotation)
colnames(mcclay_annotation)
colnames(forrest_annotation)

## Extract TCF4 genes
## If your gene column is not SYMBOL, replace SYMBOL with the correct column name
npc_tcf4_genes <- npc_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

mcclay_tcf4_genes <- mcclay_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

forrest_tcf4_genes <- forrest_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

length(npc_tcf4_genes)
length(mcclay_tcf4_genes)
length(forrest_tcf4_genes)

####################################################################
# 2. Import New Blake Gene-Level Results
####################################################################
blake_gene_level_results <- read.csv(
  "Blake_GSE48367_limma_all_results_probe_level.csv"
)

colnames(blake_gene_level_results)

## Blake background = all genes tested in the new Blake analysis

blake_background_genes <- blake_gene_level_results$gene %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

length(blake_background_genes)


# 3. Create New Blake DEG Lists
## New Blake DEGs at FDR < 0.05
new_blake_all_FDR005 <- blake_gene_level_results %>%
  filter(adj.P.Val < 0.05) %>%
  pull(gene) %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

new_blake_up_FDR005 <- blake_gene_level_results %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  pull(gene) %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

new_blake_down_FDR005 <- blake_gene_level_results %>%
  filter(adj.P.Val < 0.05, logFC < 0) %>%
  pull(gene) %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()


## New Blake DEGs at FDR < 0.01
new_blake_all_FDR001 <- blake_gene_level_results %>%
  filter(adj.P.Val < 0.01) %>%
  pull(gene) %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

new_blake_up_FDR001 <- blake_gene_level_results %>%
  filter(adj.P.Val < 0.01, logFC > 0) %>%
  pull(gene) %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

new_blake_down_FDR001 <- blake_gene_level_results %>%
  filter(adj.P.Val < 0.01, logFC < 0) %>%
  pull(gene) %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()


## Check new Blake DEG counts
new_blake_count_summary <- data.frame(
  Blake_Set = c(
    "New_Blake_all_FDR005",
    "New_Blake_up_FDR005",
    "New_Blake_down_FDR005",
    "New_Blake_all_FDR001",
    "New_Blake_up_FDR001",
    "New_Blake_down_FDR001"
  ),
  Gene_Count = c(
    length(new_blake_all_FDR005),
    length(new_blake_up_FDR005),
    length(new_blake_down_FDR005),
    length(new_blake_all_FDR001),
    length(new_blake_up_FDR001),
    length(new_blake_down_FDR001)
  )
)

new_blake_count_summary


# 4. Restrict TCF4 and Blake DEG Lists to Blake Background
npc_in_background <- intersect(npc_tcf4_genes, blake_background_genes)
mcclay_in_background <- intersect(mcclay_tcf4_genes, blake_background_genes)
forrest_in_background <- intersect(forrest_tcf4_genes, blake_background_genes)

new_blake_all_FDR005_in_background <- intersect(new_blake_all_FDR005, blake_background_genes)
new_blake_up_FDR005_in_background <- intersect(new_blake_up_FDR005, blake_background_genes)
new_blake_down_FDR005_in_background <- intersect(new_blake_down_FDR005, blake_background_genes)

new_blake_all_FDR001_in_background <- intersect(new_blake_all_FDR001, blake_background_genes)
new_blake_up_FDR001_in_background <- intersect(new_blake_up_FDR001, blake_background_genes)
new_blake_down_FDR001_in_background <- intersect(new_blake_down_FDR001, blake_background_genes)

length(npc_in_background)
length(mcclay_in_background)
length(forrest_in_background)

length(new_blake_all_FDR005_in_background)
length(new_blake_up_FDR005_in_background)
length(new_blake_down_FDR005_in_background)

length(new_blake_all_FDR001_in_background)
length(new_blake_up_FDR001_in_background)
length(new_blake_down_FDR001_in_background)



# 5. Simple Fisher Tests: NPC TCF4
## NPC vs new Blake all DEGs FDR < 0.05
a <- length(intersect(npc_in_background, new_blake_all_FDR005_in_background))
b <- length(npc_in_background) - a
c <- length(new_blake_all_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

npc_all_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
npc_all_FDR005_fisher <- fisher.test(npc_all_FDR005_table, alternative = "greater")


## NPC vs new Blake upregulated FDR < 0.05
a <- length(intersect(npc_in_background, new_blake_up_FDR005_in_background))
b <- length(npc_in_background) - a
c <- length(new_blake_up_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

npc_up_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
npc_up_FDR005_fisher <- fisher.test(npc_up_FDR005_table, alternative = "greater")


## NPC vs new Blake downregulated FDR < 0.05
a <- length(intersect(npc_in_background, new_blake_down_FDR005_in_background))
b <- length(npc_in_background) - a
c <- length(new_blake_down_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

npc_down_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
npc_down_FDR005_fisher <- fisher.test(npc_down_FDR005_table, alternative = "greater")


## NPC vs new Blake all DEGs FDR < 0.01
a <- length(intersect(npc_in_background, new_blake_all_FDR001_in_background))
b <- length(npc_in_background) - a
c <- length(new_blake_all_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

npc_all_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
npc_all_FDR001_fisher <- fisher.test(npc_all_FDR001_table, alternative = "greater")


## NPC vs new Blake upregulated FDR < 0.01
a <- length(intersect(npc_in_background, new_blake_up_FDR001_in_background))
b <- length(npc_in_background) - a
c <- length(new_blake_up_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

npc_up_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
npc_up_FDR001_fisher <- fisher.test(npc_up_FDR001_table, alternative = "greater")


## NPC vs new Blake downregulated FDR < 0.01
a <- length(intersect(npc_in_background, new_blake_down_FDR001_in_background))
b <- length(npc_in_background) - a
c <- length(new_blake_down_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

npc_down_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
npc_down_FDR001_fisher <- fisher.test(npc_down_FDR001_table, alternative = "greater")


# 6. Simple Fisher Tests: McClay TCF4
## McClay vs FDR < 0.05
a <- length(intersect(mcclay_in_background, new_blake_all_FDR005_in_background))
b <- length(mcclay_in_background) - a
c <- length(new_blake_all_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

mcclay_all_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
mcclay_all_FDR005_fisher <- fisher.test(mcclay_all_FDR005_table, alternative = "greater")


a <- length(intersect(mcclay_in_background, new_blake_up_FDR005_in_background))
b <- length(mcclay_in_background) - a
c <- length(new_blake_up_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

mcclay_up_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
mcclay_up_FDR005_fisher <- fisher.test(mcclay_up_FDR005_table, alternative = "greater")


a <- length(intersect(mcclay_in_background, new_blake_down_FDR005_in_background))
b <- length(mcclay_in_background) - a
c <- length(new_blake_down_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

mcclay_down_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
mcclay_down_FDR005_fisher <- fisher.test(mcclay_down_FDR005_table, alternative = "greater")


## McClay vs FDR < 0.01
a <- length(intersect(mcclay_in_background, new_blake_all_FDR001_in_background))
b <- length(mcclay_in_background) - a
c <- length(new_blake_all_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

mcclay_all_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
mcclay_all_FDR001_fisher <- fisher.test(mcclay_all_FDR001_table, alternative = "greater")


a <- length(intersect(mcclay_in_background, new_blake_up_FDR001_in_background))
b <- length(mcclay_in_background) - a
c <- length(new_blake_up_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

mcclay_up_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
mcclay_up_FDR001_fisher <- fisher.test(mcclay_up_FDR001_table, alternative = "greater")


a <- length(intersect(mcclay_in_background, new_blake_down_FDR001_in_background))
b <- length(mcclay_in_background) - a
c <- length(new_blake_down_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

mcclay_down_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
mcclay_down_FDR001_fisher <- fisher.test(mcclay_down_FDR001_table, alternative = "greater")

####################################################################
# 7. Simple Fisher Tests: Forrest TCF4
####################################################################

## Forrest vs FDR < 0.05
a <- length(intersect(forrest_in_background, new_blake_all_FDR005_in_background))
b <- length(forrest_in_background) - a
c <- length(new_blake_all_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

forrest_all_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
forrest_all_FDR005_fisher <- fisher.test(forrest_all_FDR005_table, alternative = "greater")


a <- length(intersect(forrest_in_background, new_blake_up_FDR005_in_background))
b <- length(forrest_in_background) - a
c <- length(new_blake_up_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

forrest_up_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
forrest_up_FDR005_fisher <- fisher.test(forrest_up_FDR005_table, alternative = "greater")


a <- length(intersect(forrest_in_background, new_blake_down_FDR005_in_background))
b <- length(forrest_in_background) - a
c <- length(new_blake_down_FDR005_in_background) - a
d <- length(blake_background_genes) - a - b - c

forrest_down_FDR005_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
forrest_down_FDR005_fisher <- fisher.test(forrest_down_FDR005_table, alternative = "greater")


## Forrest vs FDR < 0.01
a <- length(intersect(forrest_in_background, new_blake_all_FDR001_in_background))
b <- length(forrest_in_background) - a
c <- length(new_blake_all_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

forrest_all_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
forrest_all_FDR001_fisher <- fisher.test(forrest_all_FDR001_table, alternative = "greater")


a <- length(intersect(forrest_in_background, new_blake_up_FDR001_in_background))
b <- length(forrest_in_background) - a
c <- length(new_blake_up_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

forrest_up_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
forrest_up_FDR001_fisher <- fisher.test(forrest_up_FDR001_table, alternative = "greater")


a <- length(intersect(forrest_in_background, new_blake_down_FDR001_in_background))
b <- length(forrest_in_background) - a
c <- length(new_blake_down_FDR001_in_background) - a
d <- length(blake_background_genes) - a - b - c

forrest_down_FDR001_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
forrest_down_FDR001_fisher <- fisher.test(forrest_down_FDR001_table, alternative = "greater")


# 8. Make Fisher Summary Table
new_blake_fisher_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC_TCF4", "NPC_TCF4", "NPC_TCF4",
    "NPC_TCF4", "NPC_TCF4", "NPC_TCF4",
    "McClay_TCF4", "McClay_TCF4", "McClay_TCF4",
    "McClay_TCF4", "McClay_TCF4", "McClay_TCF4",
    "Forrest_TCF4", "Forrest_TCF4", "Forrest_TCF4",
    "Forrest_TCF4", "Forrest_TCF4", "Forrest_TCF4"
  ),

  Blake_Set = c(
    "New_Blake_all_FDR005",
    "New_Blake_up_FDR005",
    "New_Blake_down_FDR005",
    "New_Blake_all_FDR001",
    "New_Blake_up_FDR001",
    "New_Blake_down_FDR001",

    "New_Blake_all_FDR005",
    "New_Blake_up_FDR005",
    "New_Blake_down_FDR005",
    "New_Blake_all_FDR001",
    "New_Blake_up_FDR001",
    "New_Blake_down_FDR001",

    "New_Blake_all_FDR005",
    "New_Blake_up_FDR005",
    "New_Blake_down_FDR005",
    "New_Blake_all_FDR001",
    "New_Blake_up_FDR001",
    "New_Blake_down_FDR001"
  ),

  Background_Genes = rep(length(blake_background_genes), 18),

  TCF4_Genes_in_Background = c(
    rep(length(npc_in_background), 6),
    rep(length(mcclay_in_background), 6),
    rep(length(forrest_in_background), 6)
  ),

  Blake_Genes_in_Background = c(
    length(new_blake_all_FDR005_in_background),
    length(new_blake_up_FDR005_in_background),
    length(new_blake_down_FDR005_in_background),
    length(new_blake_all_FDR001_in_background),
    length(new_blake_up_FDR001_in_background),
    length(new_blake_down_FDR001_in_background),

    length(new_blake_all_FDR005_in_background),
    length(new_blake_up_FDR005_in_background),
    length(new_blake_down_FDR005_in_background),
    length(new_blake_all_FDR001_in_background),
    length(new_blake_up_FDR001_in_background),
    length(new_blake_down_FDR001_in_background),

    length(new_blake_all_FDR005_in_background),
    length(new_blake_up_FDR005_in_background),
    length(new_blake_down_FDR005_in_background),
    length(new_blake_all_FDR001_in_background),
    length(new_blake_up_FDR001_in_background),
    length(new_blake_down_FDR001_in_background)
  ),

  Overlap_Genes = c(
    length(intersect(npc_in_background, new_blake_all_FDR005_in_background)),
    length(intersect(npc_in_background, new_blake_up_FDR005_in_background)),
    length(intersect(npc_in_background, new_blake_down_FDR005_in_background)),
    length(intersect(npc_in_background, new_blake_all_FDR001_in_background)),
    length(intersect(npc_in_background, new_blake_up_FDR001_in_background)),
    length(intersect(npc_in_background, new_blake_down_FDR001_in_background)),

    length(intersect(mcclay_in_background, new_blake_all_FDR005_in_background)),
    length(intersect(mcclay_in_background, new_blake_up_FDR005_in_background)),
    length(intersect(mcclay_in_background, new_blake_down_FDR005_in_background)),
    length(intersect(mcclay_in_background, new_blake_all_FDR001_in_background)),
    length(intersect(mcclay_in_background, new_blake_up_FDR001_in_background)),
    length(intersect(mcclay_in_background, new_blake_down_FDR001_in_background)),

    length(intersect(forrest_in_background, new_blake_all_FDR005_in_background)),
    length(intersect(forrest_in_background, new_blake_up_FDR005_in_background)),
    length(intersect(forrest_in_background, new_blake_down_FDR005_in_background)),
    length(intersect(forrest_in_background, new_blake_all_FDR001_in_background)),
    length(intersect(forrest_in_background, new_blake_up_FDR001_in_background)),
    length(intersect(forrest_in_background, new_blake_down_FDR001_in_background))
  ),

  Odds_Ratio = c(
    npc_all_FDR005_fisher$estimate,
    npc_up_FDR005_fisher$estimate,
    npc_down_FDR005_fisher$estimate,
    npc_all_FDR001_fisher$estimate,
    npc_up_FDR001_fisher$estimate,
    npc_down_FDR001_fisher$estimate,

    mcclay_all_FDR005_fisher$estimate,
    mcclay_up_FDR005_fisher$estimate,
    mcclay_down_FDR005_fisher$estimate,
    mcclay_all_FDR001_fisher$estimate,
    mcclay_up_FDR001_fisher$estimate,
    mcclay_down_FDR001_fisher$estimate,

    forrest_all_FDR005_fisher$estimate,
    forrest_up_FDR005_fisher$estimate,
    forrest_down_FDR005_fisher$estimate,
    forrest_all_FDR001_fisher$estimate,
    forrest_up_FDR001_fisher$estimate,
    forrest_down_FDR001_fisher$estimate
  ),

  P_Value = c(
    npc_all_FDR005_fisher$p.value,
    npc_up_FDR005_fisher$p.value,
    npc_down_FDR005_fisher$p.value,
    npc_all_FDR001_fisher$p.value,
    npc_up_FDR001_fisher$p.value,
    npc_down_FDR001_fisher$p.value,

    mcclay_all_FDR005_fisher$p.value,
    mcclay_up_FDR005_fisher$p.value,
    mcclay_down_FDR005_fisher$p.value,
    mcclay_all_FDR001_fisher$p.value,
    mcclay_up_FDR001_fisher$p.value,
    mcclay_down_FDR001_fisher$p.value,

    forrest_all_FDR005_fisher$p.value,
    forrest_up_FDR005_fisher$p.value,
    forrest_down_FDR005_fisher$p.value,
    forrest_all_FDR001_fisher$p.value,
    forrest_up_FDR001_fisher$p.value,
    forrest_down_FDR001_fisher$p.value
  )
)

new_blake_fisher_summary <- new_blake_fisher_summary %>%
  mutate(
    Fisher_FDR = p.adjust(P_Value, method = "BH"),

    Percent_TCF4_Overlap = round(
      Overlap_Genes / TCF4_Genes_in_Background * 100,
      3
    ),

    Percent_Blake_Overlap = round(
      Overlap_Genes / Blake_Genes_in_Background * 100,
      3
    )
  )

new_blake_fisher_summary
