# Blake TCF4 knockdown dataset
# https://pmc.ncbi.nlm.nih.gov/articles/PMC3751932/
# GEO dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48367
# Goal:Extract all detected gene symbols from the raw Blake microarray files.
# These are all genes represented/detected in the dataset, not DE genes.

# extract  all detected genes by Blake
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# Import tar file
blake_gz <- "Data integration/BlakeFiles"

# List all txt.gz files
files <- list.files(
  blake_gz,
  pattern = "\\.txt\\.gz$",
  full.names = TRUE
)

# Inspect files
# Read one file first
test_file <- read.delim(
  files[1],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(test_file) # View column names
head(test_file)     # View first few rows

# Now, I want to extract all genes
# In this dataset, the gene symbol column is named "symbol"
# Some rows may contain missing values, empty strings, or "-"
# Create an empty dataframe
all_ref_genes <- data.frame(gene = character())

# Loop through all files
for (i in seq_along(files)) {     # seq_along(): create numbers for each file, so we say go through all files one by one
  # Read one sample file
  sample_data <- read.delim(      # read.delim(): Reads a tab-separated text file because GEO files are tab-delimited microarray files
    files[i],
    stringsAsFactors = FALSE,     # keep text as character, not factor
    check.names = FALSE           # keep original column names exactly
  )

  # Extract gene symbols from this sample
  sample_genes <- sample_data %>%
    select(symbol) %>%
    rename(gene = symbol) %>%
    mutate(gene = str_trim(gene)) %>%   # str_trim(): Removes extra spaces
    filter(
      !is.na(gene),
      gene != "",
      gene != "-"
    ) %>%
    distinct()

  # Add this sample's genes to the full list
  all_ref_genes <- bind_rows(all_ref_genes, sample_genes)
}

# Create final unique Blake reference gene list
blake_ref_genes <- all_ref_genes %>%
  mutate(gene = str_trim(gene)) %>%
  filter(
    !is.na(gene),
    gene != "",
    gene != "-"
  ) %>%
  distinct(gene) %>%
  arrange(gene)

# Check number of unique genes
nrow(blake_ref_genes)

# View first genes
head(blake_ref_genes)

# Create a vector for Fisher tests
blake_ref_gene_list <- blake_ref_genes$gene

# Now I want to find the DEGs and then upregulated and downregulated genes in blake KD samples
# so this is our new mission:
# Blake GSE48367: Differential expression analysis
# Goal: TCF4 knockdown vs controls
# Controls: mock + GAPDH KD
# TCF4 knockdown: KD1 + KD2
# Gene symbol column: symbol
# Each .txt.gz file contains: 1) probes 2) gene names 3) expression intensity value

# The pipeline overview:
# 1. Read all files
# 2. Build expression matrix: we transform the data into matrix structure. this is what limma needs
# 3. Log2 transform: Microarray values are usually huge numbers. So, statistical models work better with log2(values)
# 4. Assign groups: define controls and TCF4 knockdown
# 5. limma model: For EVERY probe/gene: limma asks: “Is expression different between TCF4 KD and controls?”
#    limma calculates fold change, t-statistics, p-value, and FDR
# 6. Multiple testing correction: create adj.P.Val
# 7. Define DEGs: adj.P.Val < 0.05 -> blake_DEG
# 8. Split UP vs DOWN genes

library(GEOquery)
library(tibble)
library(limma)

# Step 1: Build expression matrix
# Store expression values from all samples in a table (dataframe)
expr_list <- data.frame()

# This table will store probe annotation information
annotation_table <- data.frame()

for (i in seq_along(files)) {       # seq_along(files) creates 1, 2, 3, ... for each file

  # Read one raw sample file
  sample_data <- read.delim(
    files[i],
    stringsAsFactors = FALSE,     # keep text as character, not factor
    check.names = FALSE           # keep original column names exactly
  )

  # Create a clean sample name
  sample_name <- basename(files[i])     # basename() removes the folder path
  sample_name <- str_remove(sample_name, "\\.txt\\.gz$")  # str_remove() removes .txt.gz from the sample name

  # Extract expression values
  sample_expr <- sample_data %>%
    select(ID, symbol, `global normalization`) %>%
    rename(
      probe_id = ID,
      gene = symbol,
      expression = `global normalization`
    ) %>%
    mutate(
      gene = str_trim(gene),
      sample = sample_name
    )

  # Add this sample's expression data to the full expression table
  expr_list <- bind_rows(expr_list, sample_expr)

  # Extract annotation information
  sample_annot <- sample_data %>%
    select(ID, symbol, description) %>%
    rename(
      probe_id = ID,
      gene = symbol
    ) %>%
    mutate(gene = str_trim(gene)) %>%
    distinct()    # Remove duplicated rows

  # Add this sample's annotation data to the full annotation table
  annotation_table <- bind_rows(annotation_table, sample_annot)
}

# Keep only one annotation row per probe
annotation_table <- annotation_table %>%
  distinct(probe_id, .keep_all = TRUE)

# Step 2: Convert long expression table to matrix format
# limma needs an expression matrix:
# rows = probes
# columns = samples
# Make sure expression column is numeric
expr_list <- expr_list %>%
  mutate(expression = as.numeric(expression))

# Convert long expression table to wide format
expr_wide <- expr_list %>%
  select(probe_id, sample, expression) %>%
  group_by(probe_id, sample) %>%
  summarise(
    expression = mean(expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = expression
  )

# Convert to dataframe
expr_wide <- as.data.frame(expr_wide)

# Set probe IDs as row names
rownames(expr_wide) <- expr_wide$probe_id

# Remove probe_id column
expr_wide$probe_id <- NULL

# Convert to numeric matrix
expr_matrix <- as.matrix(expr_wide)

# Check matrix
dim(expr_matrix)

expr_matrix[1:5, 1:5]

# Step 3: Log2 transform if needed
# Microarray intensities are sometimes already log2 transformed.
# Checks the range of values to decide whether log2 transformation is needed.
qx <- as.numeric(
  quantile(
    expr_matrix,
    c(0, 0.25, 0.5, 0.75, 0.99, 1),
    na.rm = TRUE
  )
)

# If values are very large, they are likely not log2 transformed
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0)


# Apply log2 transformation only if needed
if (LogC) {

  # log2 cannot handle zero or negative values
  expr_matrix[expr_matrix <= 0] <- NA

  # log2 transform expression values
  expr_matrix <- log2(expr_matrix)
}


# Remove probes with missing values
expr_matrix <- na.omit(expr_matrix)

# Check final matrix size
dim(expr_matrix)

# Step 4: Create sample group labels
# Download GEO metadata
gset <- getGEO("GSE48367", GSEMatrix = TRUE)

# If multiple platforms exist
if (length(gset) > 1) {
  idx <- grep("GPL13915", attr(gset, "names"))
} else {
  idx <- 1
}

gset <- gset[[idx]]

# Sample metadata
pheno <- pData(gset)

# Check columns
colnames(pheno)

# Look at important metadata
pheno[, c("geo_accession", "title", "source_name_ch1")]

# Create sample group labels from GEO metadata
sample_info <- pheno %>%
  select(geo_accession, title, source_name_ch1) %>%
  mutate(
    group = case_when(
      str_detect(title, regex("Mock", ignore_case = TRUE)) ~ "Control",
      str_detect(title, regex("GAPDH", ignore_case = TRUE)) ~ "Control",
      str_detect(title, regex("TCF4 KD1|TCF4 KD2", ignore_case = TRUE)) ~ "TCF4_KD",
      TRUE ~ NA_character_
    )
  )

sample_info
table(sample_info$group)

# Step 5: Match GEO sample metadata to expression matrix columns
# The expression matrix column names are longer than the GEO IDs
# The GEO metadata uses only:GSM1176499
# So we extract the GSM ID from each expression matrix column
sample_info_final <- data.frame(
  sample = colnames(expr_matrix)
) %>%
  mutate(
    geo_accession = str_extract(sample, "GSM[0-9]+")
  ) %>%
  left_join(sample_info, by = "geo_accession")


# Check that each expression matrix sample received the correct group label
sample_info_final
table(sample_info_final$group)

# Check if any samples failed to match
sample_info_final %>%
  filter(is.na(group))

# Step 6: Create design matrix for limma
# Convert group column to a factor
group <- factor(sample_info_final$group)

# Set Control as the reference group
group <- relevel(group, ref = "Control")

# Create design matrix
# This tells limma which samples belong to each group.
design <- model.matrix(~ 0 + group)

# Clean column names
colnames(design) <- levels(group)

# View the design matrix
design

# Step 7: Fit limma model
# lmFit() fits a linear model for every probe in the expression matrix.
# Rows = probes
# Columns = samples
fit <- lmFit(expr_matrix, design)

# Define the contrast/comparison:
# TCF4 knockdown minus Control
# Positive logFC means:
# higher expression in TCF4_KD compared with Control
# Negative logFC means:
# lower expression in TCF4_KD compared with Control

contrast_matrix <- makeContrasts(
  TCF4_KD_vs_Control = TCF4_KD - Control,
  levels = design
)

contrast_matrix

# Apply the contrast to the fitted model
fit2 <- contrasts.fit(fit, contrast_matrix)

# eBayes() applies empirical Bayes moderation.
# This improves variance estimation, especially with small sample size.
fit2 <- eBayes(fit2)


# Step 8: Extract differential expression results
# topTable() extracts the limma results.
# number = Inf means return all probes.
# adjust.method = "BH" means Benjamini-Hochberg FDR correction.

blake_results <- topTable(
  fit2,
  coef = "TCF4_KD_vs_Control",
  number = Inf,
  adjust.method = "BH"
)

# View first rows
head(blake_results)

# Check number of probes tested
nrow(blake_results)


# Step 9: Add gene symbols and descriptions
# blake_results has probe IDs as row names.
# We convert row names into a column called probe_id,
# then join annotation_table to add gene symbols and descriptions.
blake_results_clean <- blake_results %>%
  rownames_to_column("probe_id") %>%
  left_join(annotation_table, by = "probe_id") %>%
  filter(
    !is.na(gene),
    gene != "",
    gene != "-"
  )

# View cleaned results
head(blake_results_clean)

# Check number of annotated probes
nrow(blake_results_clean)

# Step 10: Identify significant DEGs
# adj.P.Val is the FDR-adjusted p-value.
# FDR < 0.05 is the standard cutoff for significant DEGs.
blake_DEG <- blake_results_clean %>%
  filter(adj.P.Val < 0.05)


# Upregulated genes after TCF4 knockdown
# logFC > 0 means higher in TCF4_KD compared with Control.
blake_DEG_up <- blake_DEG %>%
  filter(logFC > 0)


# Downregulated genes after TCF4 knockdown
# logFC < 0 means lower in TCF4_KD compared with Control.
blake_DEG_down <- blake_DEG %>%
  filter(logFC < 0)


# Count unique genes
n_distinct(blake_DEG$gene)
n_distinct(blake_DEG_up$gene)
n_distinct(blake_DEG_down$gene)

# Step 11: Create clean gene lists
# All significant DEG genes
blake_DEG_gene_list <- blake_DEG %>%
  distinct(gene) %>%
  arrange(gene)

# Upregulated significant genes
blake_DEG_up_gene_list <- blake_DEG_up %>%
  distinct(gene) %>%
  arrange(gene)

# Downregulated significant genes
blake_DEG_down_gene_list <- blake_DEG_down %>%
  distinct(gene) %>%
  arrange(gene)


# Preview lists
head(blake_DEG_gene_list)
head(blake_DEG_up_gene_list)
head(blake_DEG_down_gene_list)

# Step 12: Save outputs
write.csv(
  blake_results_clean,
  "Blake_GSE48367_limma_all_results.csv",
  row.names = FALSE
)

write.csv(
  blake_DEG,
  "Blake_DEG.csv",
  row.names = FALSE
)

write.csv(
  blake_DEG_gene_list,
  "Blake_DEG_gene_list.csv",
  row.names = FALSE
)

write.csv(
  blake_DEG_up_gene_list,
  "Blake_up_gene_list.csv",
  row.names = FALSE
)

write.csv(
  blake_DEG_down_gene_list,
  "Blake_DEG_down_gene_list.csv",
  row.names = FALSE
)

# Make stricter FDR (FDR < 0.01)
blake_DEG_FDR001 <- blake_results_clean %>%
  filter(adj.P.Val < 0.01)

blake_DEG_up_FDR001 <- blake_DEG_FDR001 %>%
  filter(logFC > 0)

blake_DEG_down_FDR001 <- blake_DEG_FDR001 %>%
  filter(logFC < 0)

n_distinct(blake_DEG_FDR001$gene)
n_distinct(blake_DEG_up_FDR001$gene)
n_distinct(blake_DEG_down_FDR001$gene)

# Clean gene list with FDR < 0.01
blake_DEG_gene_list_FDR001 <- blake_DEG_FDR001 %>%
  distinct(gene) %>%
  arrange(gene)

blake_DEG_up_gene_list_FDR001 <- blake_DEG_up_FDR001 %>%
  distinct(gene) %>%
  arrange(gene)

blake_DEG_down_gene_list_FDR001 <- blake_DEG_down_FDR001 %>%
  distinct(gene) %>%
  arrange(gene)

# I have an older Blake's dataset and I want to compare it with new one with FDR < 0.05
#Load genes differentially expressed in Blake TCF4 siRNA paper in SH-SY5Y
blake_all <- readLines("Data integration/TCF4 KD - 1205 DEG.txt")
blake_up <- readLines("Data integration/TCF4 blake up n470.txt")
blake_down <- readLines("Data integration/TCF4 blake down n665.txt")

# Clean old Blake gene lists
old_blake_all <- blake_all %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

old_blake_up <- blake_up %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

old_blake_down <- blake_down %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

# Clean new Blake gene lists
new_blake_all <- blake_DEG_gene_list$gene %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

new_blake_up <- blake_DEG_up_gene_list$gene %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

new_blake_down <- blake_DEG_down_gene_list$gene %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

# Compare DEG counts
comparison_counts <- data.frame(
    Category = c("All DEGs", "Upregulated", "Downregulated"),
  Old_Blake = c(
    length(old_blake_all),
    length(old_blake_up),
    length(old_blake_down)
  ),
  New_Blake = c(
    length(new_blake_all),
    length(new_blake_up),
    length(new_blake_down)
  )
)

comparison_counts

# Find overlaps
overlap_all <- intersect(old_blake_all, new_blake_all)
overlap_up <- intersect(old_blake_up, new_blake_up)
overlap_down <- intersect(old_blake_down, new_blake_down)

only_old_all <- setdiff(old_blake_all, new_blake_all)
only_new_all <- setdiff(new_blake_all, old_blake_all)

only_old_up <- setdiff(old_blake_up, new_blake_up)
only_new_up <- setdiff(new_blake_up, old_blake_up)

only_old_down <- setdiff(old_blake_down, new_blake_down)
only_new_down <- setdiff(new_blake_down, old_blake_down)

# Create overlap summary table
overlap_summary <- data.frame(
  Category = c("All DEGs", "Upregulated", "Downregulated"),
  Old_Count = c(
    length(old_blake_all),
    length(old_blake_up),
    length(old_blake_down)
  ),
  New_Count = c(
    length(new_blake_all),
    length(new_blake_up),
    length(new_blake_down)
  ),
  Overlap = c(
    length(overlap_all),
    length(overlap_up),
    length(overlap_down)
  ),
  Old_Only = c(
    length(only_old_all),
    length(only_old_up),
    length(only_old_down)
  ),
  New_Only = c(
    length(only_new_all),
    length(only_new_up),
    length(only_new_down)
  )
)

overlap_summary

# Direction consistency check
old_up_new_down <- intersect(old_blake_up, new_blake_down)
old_down_new_up <- intersect(old_blake_down, new_blake_up)

direction_summary <- data.frame(
  Comparison = c(
    "Old up and new down",
    "Old down and new up"
  ),
  Count = c(
    length(old_up_new_down),
    length(old_down_new_up)
  )
)

direction_summary



# Compare old Blake DEG lists with new DEG lists (FDR < 0.01)
library(dplyr)
library(stringr)

# Clean old Blake DEG lists
old_blake_all <- blake_all %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

old_blake_up <- blake_up %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

old_blake_down <- blake_down %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()


# Clean new FDR < 0.01 DEG lists
new_blake_all_FDR001 <- blake_DEG_gene_list_FDR001$gene %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

new_blake_up_FDR001 <- blake_DEG_up_gene_list_FDR001$gene %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

new_blake_down_FDR001 <- blake_DEG_down_gene_list_FDR001$gene %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()


# Compare DEG counts
comparison_counts_FDR001 <- data.frame(
  Category = c("All DEGs", "Upregulated", "Downregulated"),

  Old_Blake = c(
    length(old_blake_all),
    length(old_blake_up),
    length(old_blake_down)
  ),

  New_Blake_FDR001 = c(
    length(new_blake_all_FDR001),
    length(new_blake_up_FDR001),
    length(new_blake_down_FDR001)
  )
)

comparison_counts_FDR001

# Find overlaps
overlap_all_FDR001 <- intersect(
  old_blake_all,
  new_blake_all_FDR001
)

overlap_up_FDR001 <- intersect(
  old_blake_up,
  new_blake_up_FDR001
)

overlap_down_FDR001 <- intersect(
  old_blake_down,
  new_blake_down_FDR001
)


# Genes found only in old dataset
only_old_all_FDR001 <- setdiff(
  old_blake_all,
  new_blake_all_FDR001
)

only_old_up_FDR001 <- setdiff(
  old_blake_up,
  new_blake_up_FDR001
)

only_old_down_FDR001 <- setdiff(
  old_blake_down,
  new_blake_down_FDR001
)


# Genes found only in new dataset
only_new_all_FDR001 <- setdiff(
  new_blake_all_FDR001,
  old_blake_all
)

only_new_up_FDR001 <- setdiff(
  new_blake_up_FDR001,
  old_blake_up
)

only_new_down_FDR001 <- setdiff(
  new_blake_down_FDR001,
  old_blake_down
)

# Create overlap summary table
overlap_summary_FDR001 <- data.frame(

  Category = c(
    "All DEGs",
    "Upregulated",
    "Downregulated"
  ),

  Old_Count = c(
    length(old_blake_all),
    length(old_blake_up),
    length(old_blake_down)
  ),

  New_Count_FDR001 = c(
    length(new_blake_all_FDR001),
    length(new_blake_up_FDR001),
    length(new_blake_down_FDR001)
  ),

  Overlap = c(
    length(overlap_all_FDR001),
    length(overlap_up_FDR001),
    length(overlap_down_FDR001)
  ),

  Old_Only = c(
    length(only_old_all_FDR001),
    length(only_old_up_FDR001),
    length(only_old_down_FDR001)
  ),

  New_Only = c(
    length(only_new_all_FDR001),
    length(only_new_up_FDR001),
    length(only_new_down_FDR001)
  )
)

overlap_summary_FDR001


# Check direction consistency
# Old up but new down
old_up_new_down_FDR001 <- intersect(
  old_blake_up,
  new_blake_down_FDR001
)

# Old down but new up
old_down_new_up_FDR001 <- intersect(
  old_blake_down,
  new_blake_up_FDR001
)

direction_summary_FDR001 <- data.frame(

  Comparison = c(
    "Old up / New down",
    "Old down / New up"
  ),

  Count = c(
    length(old_up_new_down_FDR001),
    length(old_down_new_up_FDR001)
  )
)

direction_summary_FDR001

# Do enrichment analysis with FDR 5% + new background
# Blake experiment background genes

blake_background <- blake_results_clean$gene %>%
  str_trim() %>%
  .[!is.na(.) & . != "" & . != "-"] %>%
  unique()

length(blake_background)

# Create logical vectors for TCF4 ChIP gene sets
is_all_gene_new_blake <- blake_background %in% gene_all
is_motif_gene_new_blake <- blake_background %in% gene_motif
is_promoter_gene_new_blake <- blake_background %in% promoter_ebox_genes

sum(is_all_gene_new_blake)
sum(is_motif_gene_new_blake)
sum(is_promoter_gene_new_blake)

# Create logical vectors for Blake's DEG sets
# Blake DEG genes
is_blake_all <- blake_background %in% blake_DEG_gene_list$gene
is_blake_up <- blake_background %in% blake_DEG_up_gene_list$gene
is_blake_down <- blake_background %in% blake_DEG_down_gene_list$gene

sum(is_blake_all)
sum(is_blake_up)
sum(is_blake_down)

# Fisher tests for TCF4 all-peak genes
fisher_allpeaks_allgenes <- fisher.test(
  table(is_all_gene_new_blake, is_blake_all),
  alternative = "greater")

fisher_allpeaks_UPgenes <- fisher.test(
  table(is_all_gene_new_blake, is_blake_up),
  alternative = "greater")

fisher_allpeaks_DOWNgenes <- fisher.test(
  table(is_all_gene_new_blake, is_blake_down),
  alternative = "greater")

# Fisher tests for TCF4 motif genes
fisher_motifpeaks_allgenes <- fisher.test(
  table(is_motif_gene_new_blake, is_blake_all),
  alternative = "greater")

fisher_motifpeaks_UPgenes <- fisher.test(
  table(is_motif_gene_new_blake, is_blake_up),
  alternative = "greater")

fisher_motifpeaks_DOWNgenes <- fisher.test(
  table(is_motif_gene_new_blake, is_blake_down),
  alternative = "greater")

# Fisher tests for TCF4 promoter + E-box genes
fisher_promoterpeaks_allgenes <- fisher.test(
  table(is_promoter_gene_new_blake, is_blake_all),
  alternative = "greater")

fisher_promoterpeaks_UPgenes <- fisher.test(
  table(is_promoter_gene_new_blake, is_blake_up),
  alternative = "greater")

fisher_promoterpeaks_DOWNgenes <- fisher.test(
  table(is_promoter_gene_new_blake, is_blake_down),
  alternative = "greater")

# Create clean enrichment summary table
blake_enrichment_summary <- data.frame(

  Gene_Set = c(
    "TCF4_all_peaks",
    "TCF4_all_peaks",
    "TCF4_all_peaks",

    "TCF4_motif_peaks",
    "TCF4_motif_peaks",
    "TCF4_motif_peaks",

    "TCF4_promoter_Ebox",
    "TCF4_promoter_Ebox",
    "TCF4_promoter_Ebox"
  ),

  Blake_Gene_Set = c(
    "All_DEGs",
    "Upregulated",
    "Downregulated",

    "All_DEGs",
    "Upregulated",
    "Downregulated",

    "All_DEGs",
    "Upregulated",
    "Downregulated"
  ),

  Overlap = c(
    sum(is_all_gene_new_blake & is_blake_all),
    sum(is_all_gene_new_blake & is_blake_up),
    sum(is_all_gene_new_blake & is_blake_down),

    sum(is_motif_gene_new_blake & is_blake_all),
    sum(is_motif_gene_new_blake & is_blake_up),
    sum(is_motif_gene_new_blake & is_blake_down),

    sum(is_promoter_gene_new_blake & is_blake_all),
    sum(is_promoter_gene_new_blake & is_blake_up),
    sum(is_promoter_gene_new_blake & is_blake_down)
  ),

  Odds_Ratio = c(
    fisher_allpeaks_allgenes$estimate,
    fisher_allpeaks_UPgenes$estimate,
    fisher_allpeaks_DOWNgenes$estimate,

    fisher_motifpeaks_allgenes$estimate,
    fisher_motifpeaks_UPgenes$estimate,
    fisher_motifpeaks_DOWNgenes$estimate,

    fisher_promoterpeaks_allgenes$estimate,
    fisher_promoterpeaks_UPgenes$estimate,
    fisher_promoterpeaks_DOWNgenes$estimate
  ),

  P_Value = c(
    fisher_allpeaks_allgenes$p.value,
    fisher_allpeaks_UPgenes$p.value,
    fisher_allpeaks_DOWNgenes$p.value,

    fisher_motifpeaks_allgenes$p.value,
    fisher_motifpeaks_UPgenes$p.value,
    fisher_motifpeaks_DOWNgenes$p.value,

    fisher_promoterpeaks_allgenes$p.value,
    fisher_promoterpeaks_UPgenes$p.value,
    fisher_promoterpeaks_DOWNgenes$p.value
  )
)

# View summary table
blake_enrichment_summary

# ============================================================
# Enrichment analysis with Blake FDR < 0.01 DEGs
# Background = all Blake genes tested in the experiment
# ============================================================
# Create logical vectors for Blake FDR < 0.01 DEG sets
is_blake_all_FDR001 <- blake_background %in% new_blake_all_FDR001
is_blake_up_FDR001 <- blake_background %in% new_blake_up_FDR001
is_blake_down_FDR001 <- blake_background %in% new_blake_down_FDR001

sum(is_blake_all_FDR001)
sum(is_blake_up_FDR001)
sum(is_blake_down_FDR001)


# Fisher tests for TCF4 all-peak genes
fisher_allpeaks_allgenes_FDR001 <- fisher.test(
  table(is_all_gene_new_blake, is_blake_all_FDR001),
  alternative = "greater"
)

fisher_allpeaks_UPgenes_FDR001 <- fisher.test(
  table(is_all_gene_new_blake, is_blake_up_FDR001),
  alternative = "greater"
)

fisher_allpeaks_DOWNgenes_FDR001 <- fisher.test(
  table(is_all_gene_new_blake, is_blake_down_FDR001),
  alternative = "greater"
)


# Fisher tests for TCF4 motif genes
fisher_motifpeaks_allgenes_FDR001 <- fisher.test(
  table(is_motif_gene_new_blake, is_blake_all_FDR001),
  alternative = "greater"
)

fisher_motifpeaks_UPgenes_FDR001 <- fisher.test(
  table(is_motif_gene_new_blake, is_blake_up_FDR001),
  alternative = "greater"
)

fisher_motifpeaks_DOWNgenes_FDR001 <- fisher.test(
  table(is_motif_gene_new_blake, is_blake_down_FDR001),
  alternative = "greater"
)


# Fisher tests for TCF4 promoter + E-box genes
fisher_promoterpeaks_allgenes_FDR001 <- fisher.test(
  table(is_promoter_gene_new_blake, is_blake_all_FDR001),
  alternative = "greater"
)

fisher_promoterpeaks_UPgenes_FDR001 <- fisher.test(
  table(is_promoter_gene_new_blake, is_blake_up_FDR001),
  alternative = "greater"
)

fisher_promoterpeaks_DOWNgenes_FDR001 <- fisher.test(
  table(is_promoter_gene_new_blake, is_blake_down_FDR001),
  alternative = "greater"
)


# Create clean enrichment summary table
blake_enrichment_summary_FDR001 <- data.frame(

  Gene_Set = c(
    "TCF4_all_peaks",
    "TCF4_all_peaks",
    "TCF4_all_peaks",
    "TCF4_motif_peaks",
    "TCF4_motif_peaks",
    "TCF4_motif_peaks",
    "TCF4_promoter_Ebox",
    "TCF4_promoter_Ebox",
    "TCF4_promoter_Ebox"
  ),

  Blake_Gene_Set = c(
    "All_DEGs_FDR001",
    "Upregulated_FDR001",
    "Downregulated_FDR001",
    "All_DEGs_FDR001",
    "Upregulated_FDR001",
    "Downregulated_FDR001",
    "All_DEGs_FDR001",
    "Upregulated_FDR001",
    "Downregulated_FDR001"
  ),

  Overlap = c(
    sum(is_all_gene_new_blake & is_blake_all_FDR001),
    sum(is_all_gene_new_blake & is_blake_up_FDR001),
    sum(is_all_gene_new_blake & is_blake_down_FDR001),
    sum(is_motif_gene_new_blake & is_blake_all_FDR001),
    sum(is_motif_gene_new_blake & is_blake_up_FDR001),
    sum(is_motif_gene_new_blake & is_blake_down_FDR001),
    sum(is_promoter_gene_new_blake & is_blake_all_FDR001),
    sum(is_promoter_gene_new_blake & is_blake_up_FDR001),
    sum(is_promoter_gene_new_blake & is_blake_down_FDR001)
  ),

  Odds_Ratio = c(
    fisher_allpeaks_allgenes_FDR001$estimate,
    fisher_allpeaks_UPgenes_FDR001$estimate,
    fisher_allpeaks_DOWNgenes_FDR001$estimate,
    fisher_motifpeaks_allgenes_FDR001$estimate,
    fisher_motifpeaks_UPgenes_FDR001$estimate,
    fisher_motifpeaks_DOWNgenes_FDR001$estimate,
    fisher_promoterpeaks_allgenes_FDR001$estimate,
    fisher_promoterpeaks_UPgenes_FDR001$estimate,
    fisher_promoterpeaks_DOWNgenes_FDR001$estimate
  ),

  P_Value = c(
    fisher_allpeaks_allgenes_FDR001$p.value,
    fisher_allpeaks_UPgenes_FDR001$p.value,
    fisher_allpeaks_DOWNgenes_FDR001$p.value,
    fisher_motifpeaks_allgenes_FDR001$p.value,
    fisher_motifpeaks_UPgenes_FDR001$p.value,
    fisher_motifpeaks_DOWNgenes_FDR001$p.value,
    fisher_promoterpeaks_allgenes_FDR001$p.value,
    fisher_promoterpeaks_UPgenes_FDR001$p.value,
    fisher_promoterpeaks_DOWNgenes_FDR001$p.value
  )
)

# Add FDR correction across the 9 Fisher tests
blake_enrichment_summary_FDR001 <- blake_enrichment_summary_FDR001 %>%
  mutate(
    Fisher_FDR = p.adjust(P_Value, method = "BH")
  )

# View summary table
blake_enrichment_summary_FDR001

# ============================================================
# Extract overlapping genes between TCF4 sets and Blake DEGs
# FDR < 0.05
# ============================================================
# TCF4 all peaks + Blake upregulated genes
Blake_all_up_overlap <- blake_background[
  is_all_gene_new_blake & is_blake_up]

# TCF4 all peaks + Blake downregulated genes
Blake_all_down_overlap <- blake_background[
  is_all_gene_new_blake & is_blake_down]

# TCF4 motif peaks + Blake upregulated genes
Blake_motif_up_overlap <- blake_background[
  is_motif_gene_new_blake & is_blake_up]

# TCF4 motif peaks + Blake downregulated genes
Blake_motif_down_overlap <- blake_background[
  is_motif_gene_new_blake & is_blake_down]

# TCF4 promoter + E-box genes + Blake upregulated genes
Blake_promoter_up_overlap <- blake_background[
  is_promoter_gene_new_blake & is_blake_up]

# TCF4 promoter + E-box genes + Blake downregulated genes
Blake_promoter_down_overlap <- blake_background[
  is_promoter_gene_new_blake & is_blake_down]

# Create comprehensive Blake overlap dataset
Blake_overlap_all_up <- data.frame(
  gene = Blake_all_up_overlap,
  TCF4_Set = "All peaks",
  Blake_Direction = "Up"
)

Blake_overlap_all_down <- data.frame(
  gene = Blake_all_down_overlap,
  TCF4_Set = "All peaks",
  Blake_Direction = "Down"
)

Blake_overlap_motif_up <- data.frame(
  gene = Blake_motif_up_overlap,
  TCF4_Set = "Motif peaks",
  Blake_Direction = "Up"
)

Blake_overlap_motif_down <- data.frame(
  gene = Blake_motif_down_overlap,
  TCF4_Set = "Motif peaks",
  Blake_Direction = "Down"
)

Blake_overlap_promoter_up <- data.frame(
  gene = Blake_promoter_up_overlap,
  TCF4_Set = "Promoter E-box",
  Blake_Direction = "Up"
)

Blake_overlap_promoter_down <- data.frame(
  gene = Blake_promoter_down_overlap,
  TCF4_Set = "Promoter E-box",
  Blake_Direction = "Down"
)

Blake_TCF4_overlap_long <- bind_rows(
  Blake_overlap_all_up,
  Blake_overlap_all_down,
  Blake_overlap_motif_up,
  Blake_overlap_motif_down,
  Blake_overlap_promoter_up,
  Blake_overlap_promoter_down
) %>%
  distinct(gene, TCF4_Set, Blake_Direction) %>%
  arrange(TCF4_Set, Blake_Direction, gene)

View(Blake_TCF4_overlap_long)

# Add Blake statistics
Blake_TCF4_overlap_long <- Blake_TCF4_overlap_long %>%
  left_join(
    blake_DEG %>%
      select(
        gene,
        Blake_logFC = logFC,
        Blake_FDR = adj.P.Val
      ) %>%
      distinct(gene, .keep_all = TRUE),
    by = "gene"
  )

View(Blake_TCF4_overlap_long)

# Save the table
write.csv(
  Blake_TCF4_overlap_long,
  "Blake_TCF4_overlap_long_with_direction.csv",
  row.names = FALSE)
