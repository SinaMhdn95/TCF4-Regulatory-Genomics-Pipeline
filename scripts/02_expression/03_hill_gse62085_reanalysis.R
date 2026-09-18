####################################################################
# Reproduce GSE62085 TCF4 Knockdown Microarray Analysis
# Dataset:GSE62085_Non-normalized_data.txt
# Hill paper
# Goal:
# 1. Import non-normalized Illumina HT-12 v4 microarray data
# 2. Apply VST + RSN normalization using lumi
# 3. Filter probes detected in all 12 samples
# 4. Compare:
#    Control vs TCF4 siRNA #1
#    Control vs TCF4 siRNA #2
# 5. Keep shared significant probes changing in the same direction
# 6. Create DEG lists and background genes
####################################################################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("lumi")

# Load libraries
library(dplyr)
library(stringr)
library(tibble)
library(readr)
library(limma)
library(lumi)


## File name
expression_file <- "GSE62085_Non-normalized_data.txt"

## Read file manually
gse62085_raw <- read.delim(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

dim(gse62085_raw)
colnames(gse62085_raw)[1:10]

## Probe IDs
probe_ids <- gse62085_raw$ID_REF

## Intensity columns are every even column after ID_REF:
## 2, 4, 6, ..., 24
intensity_column_positions <- seq(
  from = 2,
  to = ncol(gse62085_raw),
  by = 2
)

## Detection p-value columns are every odd column after ID_REF:
## 3, 5, 7, ..., 25
detection_column_positions <- seq(
  from = 3,
  to = ncol(gse62085_raw),
  by = 2
)

intensity_column_positions
detection_column_positions

## Build matrices
expression_matrix <- as.matrix(
  gse62085_raw[, intensity_column_positions]
)

detection_p_matrix <- as.matrix(
  gse62085_raw[, detection_column_positions]
)

## Add row names
rownames(expression_matrix) <- probe_ids
rownames(detection_p_matrix) <- probe_ids

## Rename columns cleanly
sample_names <- paste0("Sample_", 1:12)

colnames(expression_matrix) <- sample_names
colnames(detection_p_matrix) <- sample_names

## Convert to numeric
mode(expression_matrix) <- "numeric"
mode(detection_p_matrix) <- "numeric"

## Check matrices
dim(expression_matrix)
dim(detection_p_matrix)

expression_matrix[1:5, 1:5]
detection_p_matrix[1:5, 1:5]


# Check Raw Expression and Detection Distributions
## Raw expression intensity distribution
raw_expression_quantiles <- quantile(
  expression_matrix,
  probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
  na.rm = TRUE
)

raw_expression_quantiles


## Detection p-value distribution
detection_pvalue_quantiles <- quantile(
  detection_p_matrix,
  probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
  na.rm = TRUE
)

detection_pvalue_quantiles


# Detection Filtering
# Keep probes detected in all 12 samples
detected_in_all_samples <- rowSums(
  detection_p_matrix < 0.05,
  na.rm = TRUE
) == ncol(detection_p_matrix)

table(detected_in_all_samples)

filtered_expression_matrix <- expression_matrix[
  detected_in_all_samples,
]

filtered_detection_matrix <- detection_p_matrix[
  detected_in_all_samples,
]

dim(filtered_expression_matrix)

probe_background <- rownames(filtered_expression_matrix)

length(probe_background)


table(detected_in_all_samples)
dim(filtered_expression_matrix)
length(probe_background)

####################################################################
# Normalize Filtered Expression Matrix
# Goal:
# 1. Log2 transform raw intensities
# 2. Normalize arrays before differential expression
####################################################################
## Log2 transform
filtered_expression_matrix[
  filtered_expression_matrix <= 0
] <- NA

log2_filtered_expression_matrix <- log2(
  filtered_expression_matrix
)

## Quantile normalization
## This is a practical replacement because lumiR() did not read this file.
normalized_expression_matrix <- normalizeBetweenArrays(
  log2_filtered_expression_matrix,
  method = "quantile"
)

## Check normalized expression distribution
normalized_expression_quantiles <- quantile(
  normalized_expression_matrix,
  probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
  na.rm = TRUE
)

normalized_expression_quantiles

dim(normalized_expression_matrix)

####################################################################
# Assign Sample Groups
# Based on GEO description:
# Sample 1-4 = Control siRNA
# Sample 5-8 = TCF4 siRNA #1
# Sample 9-12 = TCF4 siRNA #2
####################################################################
sample_info <- data.frame(
  sample = colnames(normalized_expression_matrix),
  group = c(
    rep("Control", 4),
    rep("TCF4_siRNA1", 4),
    rep("TCF4_siRNA2", 4)
  )
)

sample_info

table(sample_info$group)

sample_order_check <- data.frame(
  Matrix_Column = colnames(normalized_expression_matrix),
  Assigned_Group = sample_info$group
)

sample_order_check

####################################################################
# Two-Sample t-tests
# Control vs TCF4 siRNA #1
# Control vs TCF4 siRNA #2
####################################################################
control_samples <- sample_info$sample[
  sample_info$group == "Control"
]

sirna1_samples <- sample_info$sample[
  sample_info$group == "TCF4_siRNA1"
]

sirna2_samples <- sample_info$sample[
  sample_info$group == "TCF4_siRNA2"
]

control_matrix <- normalized_expression_matrix[, control_samples]
sirna1_matrix <- normalized_expression_matrix[, sirna1_samples]
sirna2_matrix <- normalized_expression_matrix[, sirna2_samples]

control_mean <- rowMeans(control_matrix, na.rm = TRUE)
sirna1_mean <- rowMeans(sirna1_matrix, na.rm = TRUE)
sirna2_mean <- rowMeans(sirna2_matrix, na.rm = TRUE)

sirna1_logFC <- sirna1_mean - control_mean
sirna2_logFC <- sirna2_mean - control_mean

sirna1_pvalue <- apply(
  normalized_expression_matrix,
  1,
  function(x) {
    t.test(
      x[sample_info$group == "TCF4_siRNA1"],
      x[sample_info$group == "Control"],
      alternative = "two.sided"
    )$p.value
  }
)

sirna2_pvalue <- apply(
  normalized_expression_matrix,
  1,
  function(x) {
    t.test(
      x[sample_info$group == "TCF4_siRNA2"],
      x[sample_info$group == "Control"],
      alternative = "two.sided"
    )$p.value
  }
)


# Create Probe-Level Differential Expression Results
gse62085_probe_results <- data.frame(
  probe_id = rownames(normalized_expression_matrix),

  Control_mean = control_mean,
  siRNA1_mean = sirna1_mean,
  siRNA2_mean = sirna2_mean,

  siRNA1_logFC = sirna1_logFC,
  siRNA2_logFC = sirna2_logFC,

  siRNA1_pvalue = sirna1_pvalue,
  siRNA2_pvalue = sirna2_pvalue
) %>%
  mutate(
    siRNA1_FDR_BH = p.adjust(siRNA1_pvalue, method = "BH"),
    siRNA2_FDR_BH = p.adjust(siRNA2_pvalue, method = "BH"),

    siRNA1_direction = case_when(
      siRNA1_logFC > 0 ~ "Up",
      siRNA1_logFC < 0 ~ "Down",
      TRUE ~ "No_change"
    ),

    siRNA2_direction = case_when(
      siRNA2_logFC > 0 ~ "Up",
      siRNA2_logFC < 0 ~ "Down",
      TRUE ~ "No_change"
    )
  )

head(gse62085_probe_results)


# Check DEG Counts Against Paper
gse62085_sirna1_nominal <- gse62085_probe_results %>%
  filter(siRNA1_pvalue < 0.05)

gse62085_sirna2_nominal <- gse62085_probe_results %>%
  filter(siRNA2_pvalue < 0.05)

gse62085_shared_nominal_probes <- gse62085_probe_results %>%
  filter(
    siRNA1_pvalue < 0.05,
    siRNA2_pvalue < 0.05,
    siRNA1_direction == siRNA2_direction
  ) %>%
  mutate(
    shared_direction = siRNA1_direction
  )

nrow(gse62085_sirna1_nominal)
nrow(gse62085_sirna2_nominal)
nrow(gse62085_shared_nominal_probes)

table(gse62085_shared_nominal_probes$shared_direction)


# Save Current Approximation Results
write.csv(
  gse62085_probe_results,
  "GSE62085_probe_level_results_log2_quantile.csv",
  row.names = FALSE
)

write.csv(
  gse62085_shared_nominal_probes,
  "GSE62085_shared_nominal_same_direction_log2_quantile.csv",
  row.names = FALSE
)

gse62085_diagnostic_summary_log2_quantile <- data.frame(
  Metric = c(
    "Total probes in raw file",
    "Probes detected in all 12 samples",
    "siRNA1 nominal p < 0.05 probes",
    "siRNA2 nominal p < 0.05 probes",
    "Shared same-direction probes",
    "Shared upregulated probes",
    "Shared downregulated probes"
  ),
  Value = c(
    nrow(expression_matrix),
    nrow(filtered_expression_matrix),
    nrow(gse62085_sirna1_nominal),
    nrow(gse62085_sirna2_nominal),
    nrow(gse62085_shared_nominal_probes),
    sum(gse62085_shared_nominal_probes$shared_direction == "Up"),
    sum(gse62085_shared_nominal_probes$shared_direction == "Down")
  )
)

gse62085_diagnostic_summary_log2_quantile

write.csv(
  gse62085_diagnostic_summary_log2_quantile,
  "GSE62085_diagnostic_summary_log2_quantile.csv",
  row.names = FALSE
)


####################################################################
# 13. Read GPL10558 Annotation File
####################################################################
####################################################################
# Parse GPL10558 Illumina Manifest Annotation File
####################################################################
library(dplyr)
library(stringr)
library(readr)

annotation_file_R2 <- "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"

## Read the whole annotation file as text
gpl_lines <- readLines(
  gzfile(annotation_file_R2)
)

## Check the first lines
head(gpl_lines, 30)

####################################################################
# Find Probe Annotation Section
####################################################################
probes_start <- grep(
  "^\\[Probes\\]",
  gpl_lines
)

controls_start <- grep(
  "^\\[Controls\\]",
  gpl_lines
)

probes_start
controls_start

####################################################################
# Extract Probe Annotation Lines
####################################################################
probe_lines <- gpl_lines[
  (probes_start + 1):(controls_start - 1)
]

## Check first few probe-section lines
head(probe_lines, 10)

####################################################################
# Read Probe Annotation Table Correctly
####################################################################
probe_temp_file <- tempfile(fileext = ".txt")

writeLines(
  probe_lines,
  con = probe_temp_file
)

gpl10558_probes <- read.delim(
  probe_temp_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  sep = "\t"
)

dim(gpl10558_probes)
colnames(gpl10558_probes)
head(gpl10558_probes)

####################################################################
# Clean GPL10558 Probe Annotation Column Names
####################################################################
gpl10558_clean <- gpl10558_probes %>%
  rename_with(~ str_replace_all(., " ", "_")) %>%
  rename_with(~ str_replace_all(., "-", "_")) %>%
  rename_with(~ str_replace_all(., "\\.", "_"))

colnames(gpl10558_clean)

####################################################################
# Create Probe-to-Gene Annotation Table
####################################################################
probe_annotation <- gpl10558_clean %>%
  transmute(
    probe_id = Probe_Id,
    gene_symbol = Symbol
  ) %>%
  mutate(
    probe_id = as.character(probe_id),
    gene_symbol = as.character(gene_symbol),
    probe_id = str_trim(probe_id),
    gene_symbol = str_trim(gene_symbol)
  ) %>%
  filter(
    !is.na(probe_id),
    probe_id != "",
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol != "-",
    gene_symbol != "NA",
    gene_symbol != "NULL"
  ) %>%
  distinct(probe_id, .keep_all = TRUE)

dim(probe_annotation)
head(probe_annotation)

####################################################################
# Check Probe ID Matching
####################################################################
length(
  intersect(
    gse62085_shared_nominal_probes$probe_id,
    probe_annotation$probe_id
  )
)

head(gse62085_shared_nominal_probes$probe_id)
head(probe_annotation$probe_id)

####################################################################
# Annotate Reanalyzed Shared Same-Direction Probes
####################################################################
gse62085_shared_nominal_annotated <- gse62085_shared_nominal_probes %>%
  left_join(
    probe_annotation,
    by = "probe_id"
  )

dim(gse62085_shared_nominal_annotated)

head(
  gse62085_shared_nominal_annotated %>%
    select(
      probe_id,
      gene_symbol,
      siRNA1_logFC,
      siRNA2_logFC,
      shared_direction
    )
)

sum(is.na(gse62085_shared_nominal_annotated$gene_symbol))

####################################################################
# Annotate Reanalyzed Shared Same-Direction Probes
####################################################################
gse62085_shared_nominal_annotated <- gse62085_shared_nominal_probes %>%
  left_join(
    probe_annotation,
    by = "probe_id"
  )

dim(gse62085_shared_nominal_annotated)

head(
  gse62085_shared_nominal_annotated %>%
    select(
      probe_id,
      gene_symbol,
      siRNA1_logFC,
      siRNA2_logFC,
      shared_direction
    )
)

sum(is.na(gse62085_shared_nominal_annotated$gene_symbol))

####################################################################
# Create Reanalyzed Gene Lists
####################################################################
reanalyzed_genes_all <- gse62085_shared_nominal_annotated %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol != "-",
    gene_symbol != "NA",
    gene_symbol != "NULL"
  ) %>%
  pull(gene_symbol) %>%
  as.character() %>%
  str_trim() %>%
  unique()

reanalyzed_genes_up <- gse62085_shared_nominal_annotated %>%
  filter(
    shared_direction == "Up",
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol != "-",
    gene_symbol != "NA",
    gene_symbol != "NULL"
  ) %>%
  pull(gene_symbol) %>%
  as.character() %>%
  str_trim() %>%
  unique()

reanalyzed_genes_down <- gse62085_shared_nominal_annotated %>%
  filter(
    shared_direction == "Down",
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol != "-",
    gene_symbol != "NA",
    gene_symbol != "NULL"
  ) %>%
  pull(gene_symbol) %>%
  as.character() %>%
  str_trim() %>%
  unique()

length(reanalyzed_genes_all)
length(reanalyzed_genes_up)
length(reanalyzed_genes_down)

head(reanalyzed_genes_all)

####################################################################
# Compare Reanalyzed GSE62085 Genes with Published Hill 2017 Genes
####################################################################
library(dplyr)
library(stringr)

## Import published Hill 2017 gene list
hill_published <- read.csv(
  "NPC_KD_TCF4_Hill_2017.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(hill_published)
colnames(hill_published)
head(hill_published)

####################################################################
# Check Why Hill Published List Has 624 Rows but 599 Unique Genes
####################################################################
## Total rows in the file
nrow(hill_published)

## Check gene column
colnames(hill_published)

## Number of raw entries in Gene_Symbol column
length(hill_published$Gene_Symbol)

## Number of non-empty gene symbols before unique()
hill_gene_symbols_clean_all_rows <- hill_published$Gene_Symbol %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"]

length(hill_gene_symbols_clean_all_rows)

## Number of unique gene symbols
length(unique(hill_gene_symbols_clean_all_rows))

## Number of duplicated gene-symbol entries
sum(duplicated(hill_gene_symbols_clean_all_rows))

## Which genes are duplicated?
duplicated_hill_genes <- hill_gene_symbols_clean_all_rows[
  duplicated(hill_gene_symbols_clean_all_rows) |
    duplicated(hill_gene_symbols_clean_all_rows, fromLast = TRUE)
]

sort(table(duplicated_hill_genes), decreasing = TRUE)

####################################################################
# Clean Published Hill Gene List
####################################################################
hill_published_genes <- hill_published$Gene_Symbol %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

length(hill_published_genes)
head(hill_published_genes)

####################################################################
# Compare Published Hill Genes with Reanalyzed Genes
####################################################################
shared_genes <- intersect(
  hill_published_genes,
  reanalyzed_genes_all
)

published_only_genes <- setdiff(
  hill_published_genes,
  reanalyzed_genes_all
)

reanalyzed_only_genes <- setdiff(
  reanalyzed_genes_all,
  hill_published_genes
)

length(shared_genes)
length(published_only_genes)
length(reanalyzed_only_genes)

####################################################################
# Create Hill Gene-Level Background and Non-DEG Genes
####################################################################
## Annotate all detected background probes
hill_background_probe_annotation <- data.frame(
  probe_id = probe_background
) %>%
  left_join(
    probe_annotation,
    by = "probe_id"
  )

dim(hill_background_probe_annotation)
head(hill_background_probe_annotation)

## Create Hill gene-level background
hill_background_genes <- hill_background_probe_annotation %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol != "-",
    gene_symbol != "NA",
    gene_symbol != "NULL"
  ) %>%
  pull(gene_symbol) %>%
  as.character() %>%
  str_trim() %>%
  unique()

length(hill_background_genes)
head(hill_background_genes)

####################################################################
# Restrict Published Hill DEG Genes to Hill Background
####################################################################
hill_published_genes_in_background <- intersect(
  hill_published_genes,
  hill_background_genes
)

hill_nonDGE_genes <- setdiff(
  hill_background_genes,
  hill_published_genes_in_background
)

length(hill_background_genes)
length(hill_published_genes)
length(hill_published_genes_in_background)
length(hill_nonDGE_genes)

####################################################################
# Fisher Enrichment: TCF4 Gene Sets vs Published Hill DEG Genes
####################################################################
# Reload TCF4 Gene Sets Correctly
library(dplyr)
library(stringr)

## Import TCF4 peak annotation files
npc_annotation <- read.csv(
  "npc_hg19_peak_annotation.csv",
  stringsAsFactors = FALSE
)

mcclay_annotation <- read.csv(
  "mcclay_hg19_peak_annotation.csv",
  stringsAsFactors = FALSE
)

forrest_annotation <- read.csv(
  "forrest_hg19_peak_annotation.csv",
  stringsAsFactors = FALSE
)

## Check columns
colnames(npc_annotation)
colnames(mcclay_annotation)
colnames(forrest_annotation)

####################################################################
# Extract Clean TCF4 Gene Symbols
####################################################################
npc_tcf4_genes <- npc_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  str_to_upper() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

mcclay_tcf4_genes <- mcclay_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  str_to_upper() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

forrest_tcf4_genes <- forrest_annotation$SYMBOL %>%
  as.character() %>%
  str_trim() %>%
  str_to_upper() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

length(npc_tcf4_genes)
length(mcclay_tcf4_genes)
length(forrest_tcf4_genes)

head(npc_tcf4_genes)
head(mcclay_tcf4_genes)
head(forrest_tcf4_genes)

####################################################################
# Clean Hill Background and Published Hill Genes
####################################################################
hill_background_genes_clean <- hill_background_genes %>%
  as.character() %>%
  str_trim() %>%
  str_to_upper() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

hill_published_genes_clean <- hill_published_genes %>%
  as.character() %>%
  str_trim() %>%
  str_to_upper() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

length(hill_background_genes_clean)
length(hill_published_genes_clean)

head(hill_background_genes_clean)
head(hill_published_genes_clean)

####################################################################
# Restrict Hill Published DEGs and TCF4 Gene Sets to Hill Background
####################################################################
hill_published_genes_in_background <- base::intersect(
  hill_published_genes_clean,
  hill_background_genes_clean
)

hill_nonDGE_genes <- base::setdiff(
  hill_background_genes_clean,
  hill_published_genes_in_background
)

npc_tcf4_in_hill_background <- base::intersect(
  npc_tcf4_genes,
  hill_background_genes_clean
)

mcclay_tcf4_in_hill_background <- base::intersect(
  mcclay_tcf4_genes,
  hill_background_genes_clean
)

forrest_tcf4_in_hill_background <- base::intersect(
  forrest_tcf4_genes,
  hill_background_genes_clean
)

length(hill_background_genes_clean)
length(hill_published_genes_in_background)
length(hill_nonDGE_genes)

length(npc_tcf4_in_hill_background)
length(mcclay_tcf4_in_hill_background)
length(forrest_tcf4_in_hill_background)

####################################################################
# Fisher Test: NPC TCF4 vs Published Hill DEGs
####################################################################
a <- length(base::intersect(
  npc_tcf4_in_hill_background,
  hill_published_genes_in_background
))

b <- length(npc_tcf4_in_hill_background) - a
c <- length(hill_published_genes_in_background) - a
d <- length(hill_background_genes_clean) - a - b - c

npc_hill_table <- matrix(
  c(a, b, c, d),
  nrow = 2,
  byrow = TRUE
)

rownames(npc_hill_table) <- c("NPC_TCF4_gene", "Not_NPC_TCF4_gene")
colnames(npc_hill_table) <- c("Hill_DEG", "Hill_nonDEG")

npc_hill_table

npc_hill_fisher <- fisher.test(
  npc_hill_table,
  alternative = "greater"
)

npc_hill_fisher

####################################################################
# Fisher Test: McClay TCF4 vs Published Hill DEGs
####################################################################
a <- length(base::intersect(
  mcclay_tcf4_in_hill_background,
  hill_published_genes_in_background
))

b <- length(mcclay_tcf4_in_hill_background) - a
c <- length(hill_published_genes_in_background) - a
d <- length(hill_background_genes_clean) - a - b - c

mcclay_hill_table <- matrix(
  c(a, b, c, d),
  nrow = 2,
  byrow = TRUE
)

rownames(mcclay_hill_table) <- c("McClay_TCF4_gene", "Not_McClay_TCF4_gene")
colnames(mcclay_hill_table) <- c("Hill_DEG", "Hill_nonDEG")

mcclay_hill_table

mcclay_hill_fisher <- fisher.test(
  mcclay_hill_table,
  alternative = "greater"
)

mcclay_hill_fisher

####################################################################
# Fisher Test: Forrest TCF4 vs Published Hill DEGs
####################################################################
a <- length(base::intersect(
  forrest_tcf4_in_hill_background,
  hill_published_genes_in_background
))

b <- length(forrest_tcf4_in_hill_background) - a
c <- length(hill_published_genes_in_background) - a
d <- length(hill_background_genes_clean) - a - b - c

forrest_hill_table <- matrix(
  c(a, b, c, d),
  nrow = 2,
  byrow = TRUE
)

rownames(forrest_hill_table) <- c("Forrest_TCF4_gene", "Not_Forrest_TCF4_gene")
colnames(forrest_hill_table) <- c("Hill_DEG", "Hill_nonDEG")

forrest_hill_table

forrest_hill_fisher <- fisher.test(
  forrest_hill_table,
  alternative = "greater"
)

forrest_hill_fisher

####################################################################
# Create Hill Published DEG vs TCF4 Fisher Summary Table
####################################################################
hill_tcf4_fisher_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),

  Background_Genes = rep(
    length(hill_background_genes_clean),
    3
  ),

  TCF4_Genes_in_Background = c(
    length(npc_tcf4_in_hill_background),
    length(mcclay_tcf4_in_hill_background),
    length(forrest_tcf4_in_hill_background)
  ),

  Hill_DEG_Genes_in_Background = rep(
    length(hill_published_genes_in_background),
    3
  ),

  Overlap_Genes = c(
    length(base::intersect(npc_tcf4_in_hill_background, hill_published_genes_in_background)),
    length(base::intersect(mcclay_tcf4_in_hill_background, hill_published_genes_in_background)),
    length(base::intersect(forrest_tcf4_in_hill_background, hill_published_genes_in_background))
  ),

  Odds_Ratio = c(
    npc_hill_fisher$estimate,
    mcclay_hill_fisher$estimate,
    forrest_hill_fisher$estimate
  ),

  P_Value = c(
    npc_hill_fisher$p.value,
    mcclay_hill_fisher$p.value,
    forrest_hill_fisher$p.value
  )
)

hill_tcf4_fisher_summary <- hill_tcf4_fisher_summary %>%
  mutate(
    Fisher_FDR = p.adjust(P_Value, method = "BH"),

    Percent_TCF4_Overlap = round(
      Overlap_Genes / TCF4_Genes_in_Background * 100,
      3
    ),

    Percent_Hill_DEG_Overlap = round(
      Overlap_Genes / Hill_DEG_Genes_in_Background * 100,
      3
    )
  )

hill_tcf4_fisher_summary

####################################################################
# Directional Enrichment: TCF4 Gene Sets vs Hill Published DEGs
####################################################################

library(dplyr)
library(stringr)

####################################################################
# 1. Inspect Hill Published File for Direction Columns
####################################################################

colnames(hill_published)
head(hill_published)

####################################################################
# 2A. Create Hill Up and Down Gene Lists Using Direction Column
####################################################################
hill_published_up_genes <- hill_published %>%
  filter(
    Direction == "Up"
  ) %>%
  pull(Gene_Symbol) %>%
  as.character() %>%
  str_trim() %>%
  str_to_upper() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

hill_published_down_genes <- hill_published %>%
  filter(
    Direction == "Down"
  ) %>%
  pull(Gene_Symbol) %>%
  as.character() %>%
  str_trim() %>%
  str_to_upper() %>%
  .[!is.na(.) & . != "" & . != "-" & . != "NA" & . != "NULL"] %>%
  unique()

length(hill_published_up_genes)
length(hill_published_down_genes)