# TCF4 x developing-brain cross-ancestry atlas integration
# Paper: Cross-ancestry atlas of gene, isoform, and splicing regulation in the
# developing human brain (Science 2024; eadh0829)
# https://pubmed.ncbi.nlm.nih.gov/38781368/
# The paper analyzed RNA-seq on GRCh37 with GENCODE v29lift37.
# This paper has 7 supplementary tables
# THIS CODE IS FOR TABLES1

# Load libraries
library(readxl)  # imports Excel workbooks into R
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)


# Import hg19 data annotated peaks
npc_peak_anno <- read.csv("npc_hg19_peak_annotation.csv")
mcclay_peak_anno <- read.csv("mcclay_hg19_peak_annotation.csv")
forrest_peak_anno <- read.csv("forrest_hg19_peak_annotation.csv")

# dim() reports:
#   - number of rows, representing peaks
#   - number of columns
dim(npc_peak_anno)
dim(mcclay_peak_anno)
dim(forrest_peak_anno)

# Check columns names
colnames(npc_peak_anno)
colnames(mcclay_peak_anno)
colnames(forrest_peak_anno)

# First start witrh

# Count total peaks in each TCF4 dataset
total_peak_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_Peaks = c(
    nrow(npc_peak_anno),
    nrow(mcclay_peak_anno),
    nrow(forrest_peak_anno)
  )
)

total_peak_summary

# Count promoter peaks and promoter-bound genes
# The symbol ^ means that the annotation must begin with "Promoter"
promoter_peak_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"),

  Total_Peaks = c(
    nrow(npc_peak_anno),
    nrow(mcclay_peak_anno),
    nrow(forrest_peak_anno)),

  Promoter_Peaks = c(
    sum(grepl("^Promoter",npc_peak_anno$annotation),na.rm = TRUE),

    sum(grepl("^Promoter",mcclay_peak_anno$annotation),na.rm = TRUE),

    sum(grepl("^Promoter",forrest_peak_anno$annotation),na.rm = TRUE)
  )
)

promoter_peak_summary

# Create gene sets using all annotated TCF4 peaks
# Create the NPC gene set
npc_all_genes <- npc_peak_anno %>%

  # Keep the Ensembl gene identifier and gene symbol
  transmute(

    # Remove a possible Ensembl version suffix
    # Example: ENSG000001234.5 becomes ENSG000001234
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      as.character(ENSEMBL)),

    gene_symbol = as.character(SYMBOL)
  ) %>%

  # Remove peaks without an assigned Ensembl gene
  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%

  # Multiple peaks can be assigned to the same gene
  # Keep only one row per Ensembl gene
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  ) %>%

  # Arrange genes alphabetically by Ensembl ID
  arrange(
    ensembl_gene_id
  )


# Create the McClay gene set
mcclay_all_genes <- mcclay_peak_anno %>%
  transmute(
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      as.character(ENSEMBL)
    ),

    gene_symbol = as.character(SYMBOL)
  ) %>%

  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%

  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  ) %>%

  arrange(
    ensembl_gene_id
  )


# Create the Forrest gene set
forrest_all_genes <- forrest_peak_anno %>%
  transmute(
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      as.character(ENSEMBL)
    ),

    gene_symbol = as.character(SYMBOL)
  ) %>%

  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%

  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  ) %>%

  arrange(
    ensembl_gene_id
  )

# Summarize the all-peak gene sets
all_gene_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_Peaks = c(
    nrow(npc_peak_anno),
    nrow(mcclay_peak_anno),
    nrow(forrest_peak_anno)
  ),

  Peaks_With_Ensembl_Gene = c(

    sum(
      !is.na(npc_peak_anno$ENSEMBL) &
        npc_peak_anno$ENSEMBL != ""
    ),

    sum(
      !is.na(mcclay_peak_anno$ENSEMBL) &
        mcclay_peak_anno$ENSEMBL != ""
    ),

    sum(
      !is.na(forrest_peak_anno$ENSEMBL) &
        forrest_peak_anno$ENSEMBL != ""
    )
  ),

  Unique_Peak_Associated_Genes = c(
    nrow(npc_all_genes),
    nrow(mcclay_all_genes),
    nrow(forrest_all_genes)
  )
)

all_gene_summary

# Inspect results
head(npc_all_genes)
head(mcclay_all_genes)
head(forrest_all_genes)

# Count genes with missing gene symbols
missing_symbol_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),
  Genes_With_Missing_Symbol = c(
    sum(is.na(npc_all_genes$gene_symbol) |
        npc_all_genes$gene_symbol == ""),
    sum(is.na(mcclay_all_genes$gene_symbol) |
        mcclay_all_genes$gene_symbol == ""),
    sum(is.na(forrest_all_genes$gene_symbol) |
        forrest_all_genes$gene_symbol == "")
  )
)

missing_symbol_summary

# Let's start with Table S1
# eQTL = a genetic variant associated with gene expression
# eGene = the gene affected by that eQTL
# Import the Table S1 eGene results

# Load readxl
# Select the official Table S1 Excel workbook
# A file-selection window will open
table_s1_file <- file.choose("Cross_Ancestry/science.adh0829_table_s1.xlsx")

# Import significant fetal-brain eGenes
# An eGene is a gene with evidence that genetic variation
# is associated with its total gene-expression level
s1_eGene <- read_excel(
  path = table_s1_file,
  sheet = "ST1-2-eGene")

# Display every sheet contained in Table S1
excel_sheets(table_s1_file)

dim(s1_eGene)       # Check the number of rows and columns
colnames(s1_eGene)  # Display all column names
head(s1_eGene)      # Display the first six rows

#Column	Meaning#
# pid =	Ensembl gene identifier
# gene_name	= Gene symbol
# sid	= Primary associated genetic variant
# dist = Variant distance from the TSS
# slope	= Direction and size of the eQTL effect
#npval = Nominal P value
# ppval = Permutation P value
# qval = Multiple-testing-adjusted value
# fetal_only = Whether the eQTL was reported as fetal-only
# gene_type = Gene biotype


# Integrate TCF4 peak-associated genes with Table S1 eGenes
# Create a clean Table S1 eGene table
s1_eGene_clean <- s1_eGene %>%
  transmute(

    # Remove a possible Ensembl version suffix
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      as.character(pid)
    ),

    # Gene information
    s1_gene_symbol = as.character(gene_name),
    gene_type = as.character(gene_type),

    # Primary eQTL variant
    primary_variant = as.character(sid),

    # Number of variants tested for the gene
    number_of_variants_tested = as.numeric(nvar),

    # Distance between the primary variant and the gene TSS
    distance_to_TSS_bp = as.numeric(dist),
    distance_to_TSS_kb = as.numeric(dist_kb),

    # Direction and size of the eQTL association
    eQTL_slope = as.numeric(slope),

    # Statistical results
    nominal_p_value = as.numeric(npval),
    permutation_p_value = as.numeric(ppval),
    beta_approximated_p_value = as.numeric(bpval),
    q_value = as.numeric(qval),

    # Gene-specific nominal significance threshold
    nominal_p_value_threshold = as.numeric(
      pval_nominal_threshold
    ),

    # Whether the authors classified the eQTL as fetal-only
    fetal_only = as.character(fetal_only)
  ) %>%

  # Remove rows without a usable Ensembl gene identifier
  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%

  # Keep one row per Ensembl gene
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  ) %>%

  arrange(
    ensembl_gene_id
  )

# Summarize the cleaned Table S1 eGene data
s1_eGene_QC_summary <- data.frame(
  Original_Rows = nrow(
    s1_eGene),
  Clean_Rows = nrow(
    s1_eGene_clean),
  Unique_Ensembl_Genes = dplyr::n_distinct(
    s1_eGene_clean$ensembl_gene_id),
  Missing_Ensembl_IDs = sum(
    is.na(s1_eGene_clean$ensembl_gene_id) |
      s1_eGene_clean$ensembl_gene_id == ""),
  Duplicated_Ensembl_IDs = sum(
    duplicated(
      s1_eGene_clean$ensembl_gene_id)
  )
)

s1_eGene_QC_summary

# Examine the range of adjusted q-values
summary(s1_eGene_clean$q_value)

# Count rows with q-value below 0.05
sum(s1_eGene_clean$q_value < 0.05,na.rm = TRUE)

# Now we identify the overlaps
# Find NPC TCF4 genes in the significant fetal-brain eGene list
npc_s1_eGene_overlap <- npc_all_genes %>%
  inner_join(
    s1_eGene_clean,
    by = "ensembl_gene_id"
  ) %>%
  arrange(
    q_value
  )

# Find McClay TCF4 genes in the significant fetal-brain eGene list
mcclay_s1_eGene_overlap <- mcclay_all_genes %>%
  inner_join(
    s1_eGene_clean,
    by = "ensembl_gene_id"
  ) %>%
  arrange(
    q_value
  )


# Find Forrest TCF4 genes in the significant fetal-brain eGene list
forrest_s1_eGene_overlap <- forrest_all_genes %>%
  inner_join(
    s1_eGene_clean,
    by = "ensembl_gene_id"
  ) %>%
  arrange(
    q_value
  )

# Summarize TCF4 and fetal-brain eGene overlaps
s1_eGene_overlap_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_TCF4_Peak_Associated_Genes = c(
    nrow(npc_all_genes),
    nrow(mcclay_all_genes),
    nrow(forrest_all_genes)
  ),

  Overlap_With_Significant_eGenes = c(
    nrow(npc_s1_eGene_overlap),
    nrow(mcclay_s1_eGene_overlap),
    nrow(forrest_s1_eGene_overlap)
  )
)


# Calculate the descriptive overlap percentage
s1_eGene_overlap_summary <- s1_eGene_overlap_summary %>%
  mutate(

    Percent_of_TCF4_Set_in_eGene_List = round(
      100 *
        Overlap_With_Significant_eGenes /
        Total_TCF4_Peak_Associated_Genes,
      digits = 2
    )
  )

s1_eGene_overlap_summary

# Display the first 10 overlapping genes ordered by q-value
head(npc_s1_eGene_overlap,10)
head(mcclay_s1_eGene_overlap,10)
head(forrest_s1_eGene_overlap,10)


# Compare unique and shared TCF4 gene sets
# Prepare NPC membership information
npc_membership <- npc_all_genes %>%
  transmute(
    ensembl_gene_id,
    npc_gene_symbol = gene_symbol,
    NPC = TRUE)


# Prepare McClay membership information
mcclay_membership <- mcclay_all_genes %>%
  transmute(
    ensembl_gene_id,
    mcclay_gene_symbol = gene_symbol,
    McClay = TRUE)


# Prepare Forrest membership information
forrest_membership <- forrest_all_genes %>%
  transmute(
    ensembl_gene_id,
    forrest_gene_symbol = gene_symbol,
    Forrest = TRUE)


# Combine all three TCF4 datasets
tcf4_gene_membership <- npc_membership %>%
  full_join(
    mcclay_membership,
    by = "ensembl_gene_id"
  ) %>%

  full_join(
    forrest_membership,
    by = "ensembl_gene_id"
  ) %>%

  mutate(
    # Genes missing from a dataset receive FALSE
    NPC = replace_na(
      NPC,
      FALSE
    ),
    McClay = replace_na(
      McClay,
      FALSE
    ),
    Forrest = replace_na(
      Forrest,
      FALSE
    ),

    # Select an available gene symbol
    gene_symbol = coalesce(      # selects the first non-missing value from several columns.
      npc_gene_symbol,
      mcclay_gene_symbol,
      forrest_gene_symbol)
  )

# Create the seven TCF4 membership groups
# assign every gene to one mutually exclusive membership group
tcf4_gene_membership <- tcf4_gene_membership %>%
  mutate(
    Membership_Group = case_when(
       NPC & !McClay & !Forrest ~ "NPC only",
      !NPC & McClay & !Forrest ~ "McClay only",
      !NPC & !McClay & Forrest ~ "Forrest only",
      NPC & McClay & !Forrest ~ "NPC + McClay",
      NPC & !McClay & Forrest ~ "NPC + Forrest",
      !NPC & McClay & Forrest ~ "McClay + Forrest",
      NPC & McClay & Forrest ~ "All three")
  ) %>%

  dplyr::select(
    ensembl_gene_id,
    gene_symbol,
    NPC,
    McClay,
    Forrest,
    Membership_Group
  )

# Define the order of the membership groups
membership_group_order <- c(
  "NPC only",
  "McClay only",
  "Forrest only",
  "NPC + McClay",
  "NPC + Forrest",
  "McClay + Forrest",
  "All three"
)


# Convert the membership column to an ordered factor
tcf4_gene_membership <- tcf4_gene_membership %>%
  mutate(
    Membership_Group = factor(
      Membership_Group,
      levels = membership_group_order)
  ) %>%
  arrange(
    Membership_Group,
    ensembl_gene_id)

# Count genes in each mutually exclusive group
tcf4_membership_summary <- tcf4_gene_membership %>%
  dplyr::count(
    Membership_Group,
    name = "Number_of_Genes")

tcf4_membership_summary

# Add fetal-brain eGene information
# Prepare significant fetal-brain eGene status from Table S1
s1_eGene_status <- s1_eGene_clean %>%
  transmute(
    ensembl_gene_id,
    s1_gene_symbol,

    # Every row in this table represents a significant eGene
    Significant_eGene = TRUE,
    fetal_only = tolower(
      trimws(
        as.character(fetal_only))),
    q_value,
    eQTL_slope,
    primary_variant)

# Join the eGene information to the TCF4 membership table
# Add eGene and fetal-only status to every TCF4 gene
tcf4_membership_eGene <- tcf4_gene_membership %>%
  left_join(
    s1_eGene_status,
    by = "ensembl_gene_id"
  ) %>%

  mutate(
    # Genes absent from ST1-2 are marked FALSE
    Significant_eGene = replace_na(
      Significant_eGene,
      FALSE
    ),

    # TRUE only when the gene is an eGene and fetal_only is yes
    Fetal_Only_eGene = replace_na(
      fetal_only == "yes",
      FALSE)
  )

# Now calculate eGene and fetal-only percentages
# Summarize fetal-brain regulatory evidence by membership group
tcf4_membership_eGene_summary <- tcf4_membership_eGene %>%
  group_by(
    Membership_Group) %>%

  summarise(
    Total_TCF4_Genes = n(),
    Significant_eGenes = sum(
      Significant_eGene),

    Percent_Significant_eGenes = round(
      100 *
        Significant_eGenes /
        Total_TCF4_Genes,
      digits = 2),

    Fetal_Only_eGenes = sum(
      Fetal_Only_eGene),

    Percent_Fetal_Only_of_All_Genes = round(
      100 *
        Fetal_Only_eGenes /
        Total_TCF4_Genes,
      digits = 2),

    Percent_Fetal_Only_Among_eGenes = round(
      100 *
        Fetal_Only_eGenes /
        Significant_eGenes,
      digits = 2),

    .groups = "drop")

tcf4_membership_eGene_summary

# Permutation test for significant eGenes
# Step 1: Starting table
# Step 2: Convert group names into zeros and ones
# Step 3: Give the columns readable names
# Step 4: Convert eGene status to zero and one
# Step 5: Calculate the real observed overlaps
# Step 6: Set the random seed
# Step 7: Choose 10,000 permutations (More permutations provide a more precise empirical P value but require more computation.)
# Step 8: Shuffle the eGene labels (sample() rearranges the 10,459 zero/one labels)
# Step 9: Count the randomized overlaps
# Step 10: Repeat the randomization 10,000 times
# Step 11: Calculate the mean random overlap
# Step 12: Calculate the random standard deviation
# Step 13: Calculate the empirical P value
# Step 14: Calculate observed/random ratio
# Step 15: Correct for seven tests

# Create a membership indicator matrix
# Each row represents one TCF4 gene
# Each column represents one membership group
# A value of 1 means that the gene belongs to that group
membership_matrix <- model.matrix(    # model.matrix() converts the categorical Membership_Group column
                                      # into a numeric indicator matrix containing zeros and ones
  ~ Membership_Group - 1,             # removes the intercept and tell R Do not include an intercept column
  data = tcf4_membership_eGene)


# Give the matrix columns clear names
colnames(membership_matrix) <- membership_group_order


# Convert eGene status from TRUE/FALSE to 1/0
eGene_status <- as.integer(tcf4_membership_eGene$Significant_eGene)


# Calculate the observed number of eGenes in every group
observed_eGene_overlap <- colSums(membership_matrix * eGene_status)

observed_eGene_overlap

# Set the random seed so the result can be reproduced
set.seed(12345)


# Choose the number of permutations
number_of_permutations <- 10000


# Randomly shuffle eGene status across all TCF4 genes
# replicate() repeats the enclosed analysis 10,000 times.
# sample() rearranges the eGene labels without changing
# the total number of significant eGenes.
permuted_eGene_overlaps <- replicate(
  number_of_permutations,
  {
    shuffled_eGene_status <- sample(
      eGene_status,
      size = length(eGene_status),
      replace = FALSE)

    colSums(
      membership_matrix *
        shuffled_eGene_status)
  }
)

# Summarize the eGene permutation results
eGene_membership_permutation_results <- data.frame(
  Membership_Group = membership_group_order,
  Total_Genes_in_Group = as.numeric(
    colSums(membership_matrix)
  ),

  Observed_eGenes = as.numeric(
    observed_eGene_overlap
  ),

  Mean_Random_eGenes = as.numeric(
    rowMeans(permuted_eGene_overlaps)
  ),

  Random_Standard_Deviation = as.numeric(
    apply(
      permuted_eGene_overlaps,
      MARGIN = 1,
      FUN = sd)
  ),

  Empirical_P_Value = as.numeric(
    (rowSums(permuted_eGene_overlaps >= observed_eGene_overlap) + 1)/(number_of_permutations + 1)
  )
) %>%
  mutate(
    Observed_to_Random_Ratio = round(
      Observed_eGenes /
        Mean_Random_eGenes,
      digits = 3),

    Empirical_FDR = p.adjust(
      Empirical_P_Value,
      method = "BH")
  ) %>%

  arrange(
    Empirical_FDR,
    Empirical_P_Value
  )

eGene_membership_permutation_results


# Import the Table S1 isoGene results
# Steps:
# 1. Imported significant isoform results
# 2. Cleaned the identifiers (separate the identifiers into Transcript + Gene)
# 3. Handled duplicated features (Two isoform features appeared twice with nearly identical results.we kept one with smaller qvalue)
# 4. Created a gene-level isoGene table (Because the TCF4 peak tables are currently gene-level, we summarized the isoform results by parent gene.)
# 5. Joined isoGenes to the TCF4 gene sets (comparison and overlaps)
# Import significant fetal-brain isoGenes
# An isoGene has at least one transcript or isoform whose
# abundance is significantly associated with a genetic variant.
s1_isoGene <- readxl::read_excel(path = table_s1_file,sheet = "ST1-3-isoGene")

# Check the number of rows and columns
dim(s1_isoGene)

# Display all column names
colnames(s1_isoGene)

# Display the first six rows
head(s1_isoGene)

# Inspect the Table S1 isoGene identifiers
# Display the most important columns without truncating their values
s1_isoGene %>%
  dplyr::select(
    gene_id,
    group_id,
    gene_name,
    variant_id,
    tss_distance,
    slope,
    slope_se,
    pval_nominal,
    pval_perm,
    pval_beta,
    qval,
    group_size,
    gene_type) %>%

  dplyr::slice_head(n = 10) %>%
  print(width = Inf)

# Summarize the isoGene identifiers
s1_isoGene_identifier_summary <- data.frame(

  Total_Rows = nrow(s1_isoGene),
  Unique_Transcript_Features = dplyr::n_distinct(
    s1_isoGene$gene_id,
    na.rm = TRUE),
  Unique_Ensembl_Gene_Groups = dplyr::n_distinct(
    s1_isoGene$group_id,
    na.rm = TRUE),
  Unique_Gene_Symbols = dplyr::n_distinct(
    s1_isoGene$gene_name,
    na.rm = TRUE),
  Missing_Transcript_IDs = sum(
    is.na(s1_isoGene$gene_id) |
      s1_isoGene$gene_id == ""),
  Missing_Gene_IDs = sum(
    is.na(s1_isoGene$group_id) |
      s1_isoGene$group_id == ""),
  Missing_Gene_Symbols = sum(
    is.na(s1_isoGene$gene_name) |
      s1_isoGene$gene_name == "")
)

s1_isoGene_identifier_summary

# Count significant isoform rows for each parent gene
s1_isoGene_rows_per_gene <- s1_isoGene %>%
  dplyr::count(group_id,name = "Number_of_Significant_Isoform_Rows") %>%
  dplyr::arrange(dplyr::desc(Number_of_Significant_Isoform_Rows))

head(s1_isoGene_rows_per_gene,10)

# Examine the adjusted q-value distribution
summary(s1_isoGene$qval)

# Count isoform rows with q-value below 0.05
sum(s1_isoGene$qval < 0.05,na.rm = TRUE)


# Clean the Table S1 isoGene results
# Inspect isoform feature identifiers that appear more than once
duplicated_isoGene_features <- s1_isoGene %>%
  dplyr::group_by(gene_id) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(gene_id,qval)

duplicated_isoGene_features %>%
  dplyr::select(
    gene_id,
    group_id,
    gene_name,
    variant_id,
    slope,
    pval_nominal,
    pval_perm,
    pval_beta,
    qval) %>%
  print(width = Inf)

# Create a clean transcript-level isoQTL table
s1_isoGene_clean <- s1_isoGene %>%
  dplyr::transmute(

    # Preserve the complete original feature identifier
    original_isoform_feature_id = as.character(
      gene_id),

    # Extract only the stable ENST transcript identifier
    # Example:
    # ENST00000456328.2_1:ENSG00000223972
    # becomes:
    # ENST00000456328
    ensembl_transcript_id = sub(
      "^(ENST[0-9]+).*$",
      "\\1",
      as.character(gene_id)),

    # Extract the stable parent Ensembl gene identifier
    ensembl_gene_id = sub(
      "^(ENSG[0-9]+).*$",
      "\\1",
      as.character(group_id)),

    # Gene information
    s1_gene_symbol = as.character(gene_name),

    gene_type = as.character(gene_type),

    # Number of transcript features included in the tested group
    group_size = as.numeric(group_size),

    # Primary isoQTL variant
    primary_variant = as.character(variant_id),

    # Number of genetic variants tested
    number_of_variants_tested = as.numeric(num_var),

    # Variant distance from the transcript TSS
    distance_to_transcript_TSS_bp = as.numeric(tss_distance),

    # Minor-allele information
    minor_allele_samples = as.numeric(ma_samples),

    minor_allele_count = as.numeric(ma_count),

    minor_allele_frequency = as.numeric(maf),

    # Direction and size of the variant–isoform relationship
    isoQTL_slope = as.numeric(slope),

    isoQTL_slope_standard_error = as.numeric(slope_se),

    # Statistical results
    nominal_p_value = as.numeric(pval_nominal),

    permutation_p_value = as.numeric(pval_perm),

    beta_approximated_p_value = as.numeric(pval_beta),

    q_value = as.numeric(qval),

    nominal_p_value_threshold = as.numeric(pval_nominal_threshold)
    ) %>%

  # Remove rows without valid gene or transcript identifiers
  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    !is.na(ensembl_transcript_id),
    ensembl_transcript_id != ""
  ) %>%

  # Put the strongest associations first
  dplyr::arrange(
    q_value,
    nominal_p_value
  )

# Inspect the cleaned transcript-level data
dim(s1_isoGene_clean)

head(s1_isoGene_clean,10)

# Summarize significant isoform evidence at the gene level
# Create one summary row per gene
# We need two versions:
# s1_isoGene_clean: preserves all significant transcript features;
# s1_isoGene_gene_summary: one row per gene for comparison with the TCF4 gene sets.

s1_isoGene_gene_summary <- s1_isoGene_clean %>%

  # The table is already ordered from smallest to largest q-value
  dplyr::group_by(
    ensembl_gene_id
  ) %>%
  dplyr::summarise(

    # Use the gene symbol from the strongest isoform row
    s1_gene_symbol = dplyr::first(s1_gene_symbol),

    # Count distinct significant isoform features for this gene
    Number_of_Significant_Isoform_Features =
      dplyr::n_distinct(original_isoform_feature_id),

    # Record the strongest transcript-level result
    Strongest_Transcript_ID = dplyr::first(ensembl_transcript_id),

    Strongest_Isoform_Feature_ID = dplyr::first(original_isoform_feature_id),

    Strongest_Primary_Variant = dplyr::first(primary_variant),

    Strongest_IsoQTL_Slope = dplyr::first(isoQTL_slope),

    Minimum_Isoform_Q_Value = min(q_value, na.rm = TRUE),

    # Every gene in this table has at least one significant isoform
    Significant_isoGene = TRUE, .groups = "drop") %>%

  dplyr::arrange(Minimum_Isoform_Q_Value)

# Confirm one row per parent gene
s1_isoGene_gene_QC <- data.frame(
  Total_Gene_Rows = nrow(
    s1_isoGene_gene_summary),

  Unique_Ensembl_Genes = dplyr::n_distinct(
    s1_isoGene_gene_summary$ensembl_gene_id),

  Duplicated_Ensembl_Genes = sum(
    duplicated(s1_isoGene_gene_summary$ensembl_gene_id)
  )
)

s1_isoGene_gene_QC


# Complete the observed isoGene overlap analysis
# Create one row per significant isoform feature
# The original s1_isoGene_clean object remains unchanged.
# For duplicated features, retain the row with the smallest q-value.
s1_isoGene_unique_features <- s1_isoGene_clean %>%

  dplyr::arrange(
    q_value,
    nominal_p_value
  ) %>%

  dplyr::distinct(
    original_isoform_feature_id,
    .keep_all = TRUE
  )

# Confirm that each isoform feature appears only once
s1_isoGene_unique_feature_QC <- data.frame(
  Original_Clean_Rows = nrow(s1_isoGene_clean),

  Unique_Feature_Rows = nrow(s1_isoGene_unique_features),

  Unique_Feature_IDs = dplyr::n_distinct(
    s1_isoGene_unique_features$original_isoform_feature_id),

  Duplicated_Feature_IDs = sum(
    duplicated(s1_isoGene_unique_features$original_isoform_feature_id)
  )
)

s1_isoGene_unique_feature_QC

# rebuild the gene-level summary from the unique features
# Create one summary row per gene
s1_isoGene_gene_summary <- s1_isoGene_unique_features %>%

  # The table is already ordered by q-value
  dplyr::group_by(
    ensembl_gene_id
  ) %>%

  dplyr::summarise(

    # Information from the strongest significant isoform
    s1_gene_symbol = dplyr::first(
      s1_gene_symbol
    ),

    Number_of_Significant_Isoform_Features =
      dplyr::n_distinct(
        original_isoform_feature_id
      ),

    Strongest_Transcript_ID = dplyr::first(
      ensembl_transcript_id
    ),

    Strongest_Isoform_Feature_ID = dplyr::first(
      original_isoform_feature_id
    ),

    Strongest_Primary_Variant = dplyr::first(
      primary_variant
    ),

    Strongest_IsoQTL_Slope = dplyr::first(
      isoQTL_slope
    ),

    Minimum_Isoform_Q_Value = min(
      q_value,
      na.rm = TRUE
    ),

    # Every gene in this table has at least one significant isoQTL
    Significant_isoGene = TRUE,

    .groups = "drop"
  ) %>%

  dplyr::arrange(
    Minimum_Isoform_Q_Value
  )

# Confirm one row per Ensembl gene
s1_isoGene_gene_QC <- data.frame(

  Total_Gene_Rows = nrow(
    s1_isoGene_gene_summary),

  Unique_Ensembl_Genes = dplyr::n_distinct(
    s1_isoGene_gene_summary$ensembl_gene_id),

  Duplicated_Ensembl_Genes = sum(
    duplicated(s1_isoGene_gene_summary$ensembl_gene_id)
  )
)

s1_isoGene_gene_QC

# Calculate the overall overlap for each TCF4 dataset
# Find NPC TCF4 genes with significant isoQTL evidence
npc_s1_isoGene_overlap <- npc_all_genes %>%
  dplyr::inner_join(
    s1_isoGene_gene_summary,
    by = "ensembl_gene_id"
  ) %>%

  dplyr::arrange(
    Minimum_Isoform_Q_Value
  )


# Find McClay TCF4 genes with significant isoQTL evidence
mcclay_s1_isoGene_overlap <- mcclay_all_genes %>%
  dplyr::inner_join(
    s1_isoGene_gene_summary,
    by = "ensembl_gene_id"
  ) %>%

  dplyr::arrange(
    Minimum_Isoform_Q_Value
  )


# Find Forrest TCF4 genes with significant isoQTL evidence
forrest_s1_isoGene_overlap <- forrest_all_genes %>%
  dplyr::inner_join(
    s1_isoGene_gene_summary,
    by = "ensembl_gene_id"
  ) %>%

  dplyr::arrange(
    Minimum_Isoform_Q_Value
  )

# Summarize significant isoGene overlap for each complete dataset
s1_isoGene_overlap_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_TCF4_Peak_Associated_Genes = c(
    nrow(npc_all_genes),
    nrow(mcclay_all_genes),
    nrow(forrest_all_genes)
  ),

  Overlap_With_Significant_isoGenes = c(
    nrow(npc_s1_isoGene_overlap),
    nrow(mcclay_s1_isoGene_overlap),
    nrow(forrest_s1_isoGene_overlap)
  )
) %>%

  dplyr::mutate(

    Percent_of_TCF4_Set_in_isoGene_List = round(
      100 *
        Overlap_With_Significant_isoGenes /
        Total_TCF4_Peak_Associated_Genes,
      digits = 2
    )
  )

s1_isoGene_overlap_summary

# Calculate isoGene evidence in the seven membership groups
# Add significant isoGene status to all TCF4-associated genes
tcf4_membership_isoGene <- tcf4_gene_membership %>%
  dplyr::left_join(
    s1_isoGene_gene_summary,
    by = "ensembl_gene_id"
  ) %>%

  dplyr::mutate(

    # Genes absent from the significant isoGene list receive FALSE
    Significant_isoGene = tidyr::replace_na(
      Significant_isoGene,
      FALSE
    )
  )

# Summarize isoGene evidence in the seven mutually exclusive groups
tcf4_membership_isoGene_summary <- tcf4_membership_isoGene %>%

  dplyr::group_by(
    Membership_Group
  ) %>%

  dplyr::summarise(

    Total_TCF4_Genes = dplyr::n(),

    Significant_isoGenes = sum(
      Significant_isoGene
    ),

    Percent_Significant_isoGenes = round(
      100 *
        Significant_isoGenes /
        Total_TCF4_Genes,
      digits = 2
    ),

    .groups = "drop"
  )

tcf4_membership_isoGene_summary



# Run the within-TCF4 isoGene permutation analysis
# Create a membership indicator matrix
# Each row represents one TCF4-associated gene.
# Each column represents one of the seven membership groups.
# A value of 1 means that the gene belongs to that group.
isoGene_membership_matrix <- model.matrix(
  ~ Membership_Group - 1,
  data = tcf4_membership_isoGene
)

# Apply the previously defined membership-group names
colnames(isoGene_membership_matrix) <- membership_group_order

# Check its dimensions
# 7 membership-group columns
dim(isoGene_membership_matrix)

# Display the first six rows
head(isoGene_membership_matrix)

# Convert isoGene status to numeric values
# TRUE  becomes 1
# FALSE becomes 0
isoGene_status <- as.integer(
  tcf4_membership_isoGene$Significant_isoGene
)

# Total number of TCF4-associated genes
length(isoGene_status)

# Total number of significant isoGenes in the TCF4 union
# Expected from the membership summary:
# 985 + 958 + 879 + 277 + 235 + 507 + 163 = 4,004
sum(isoGene_status)

# Calculate the observed isoGene overlaps
# Count the real significant isoGenes in each membership group
observed_isoGene_overlap <- colSums(
  isoGene_membership_matrix *
    isoGene_status)

observed_isoGene_overlap

# Set a random seed so the analysis is reproducible
set.seed(12345)


# Define the number of random permutations
number_of_isoGene_permutations <- 10000

# Randomly redistribute isoGene labels among all TCF4 genes
#
# The membership groups remain unchanged.
# The total number of significant isoGenes remains 4,004.
# Only the genes carrying the isoGene label are shuffled.
permuted_isoGene_overlaps <- replicate(
  number_of_isoGene_permutations,
  {
    # Randomly rearrange the isoGene labels
    shuffled_isoGene_status <- sample(
      isoGene_status,
      size = length(
        isoGene_status
      ),
      replace = FALSE
    )

    # Count randomly assigned isoGenes in all seven groups
    colSums(
      isoGene_membership_matrix *
        shuffled_isoGene_status
    )
  }
)

# Expected:
#   7 rows, representing membership groups
#   10,000 columns, representing permutations
dim(permuted_isoGene_overlaps)

# Summarize observed and random isoGene overlaps
isoGene_membership_permutation_results <- data.frame(
  Membership_Group = membership_group_order,
  Total_Genes_in_Group = as.numeric(
    colSums(isoGene_membership_matrix)
  ),

  Observed_isoGenes = as.numeric(observed_isoGene_overlap),

  Mean_Random_isoGenes = as.numeric(
    rowMeans(permuted_isoGene_overlaps)
  ),

  Random_Standard_Deviation = as.numeric(
    apply(
      permuted_isoGene_overlaps,
      MARGIN = 1,
      FUN = sd
    )
  ),

  Empirical_P_Value = as.numeric(
    (
      rowSums(permuted_isoGene_overlaps >=observed_isoGene_overlap) + 1) /
      (number_of_isoGene_permutations + 1)
  )
) %>%
  dplyr::mutate(

    # Difference between observed and mean random counts
    Observed_Minus_Random = round(
      Observed_isoGenes -
        Mean_Random_isoGenes,
      digits = 2
    ),

    # Fold difference relative to random expectation
    Observed_to_Random_Ratio = round(
      Observed_isoGenes /
        Mean_Random_isoGenes,
      digits = 3
    ),

    # Correct the seven empirical P values
    Empirical_FDR = p.adjust(
      Empirical_P_Value,
      method = "BH"
    ),

    # TRUE means the enrichment remains significant after correction
    Significant_After_FDR = (
      Empirical_FDR < 0.05
    )
  )

# Display the complete permutation result table
isoGene_membership_permutation_results

# Display results ordered by empirical FDR
isoGene_membership_permutation_results %>%
  dplyr::arrange(
    Empirical_FDR,
    Empirical_P_Value
  )


# Import the Table S1 sGene results
# Steps:
# 1. Imported the fetal-brain sGene table
# 2. Cleaned the identifiers and column names (This step allowes us to match the
# sGene results to our TCF4 peak-associated genes using stable Ensembl gene IDs.)
# 3. Investigated repeated splicing-feature rows
# 4. Combined splicing features by gene (Different splicing features can belong to the same gene.
# Therefore, we collaps the 7,658 features to their parent genes.)
# 5. Added sGene status to the TCF4 genes
# 6. Calculated the percentage in each group

# Import significant fetal-brain splicing QTL results
s1_sGene <- readxl::read_excel(
  path = table_s1_file,
  sheet = "ST1-4-sGene")

# Check the number of rows and columns
dim(s1_sGene)

# Display all column names
colnames(s1_sGene)

# Display the first six rows
head(s1_sGene)

# Inspect the important sGene identifiers and statistics
s1_sGene %>%
  dplyr::select(
    gene_id,
    group_id,
    variant_id,
    tss_distance,
    slope,
    slope_se,
    pval_nominal,
    pval_perm,
    pval_beta,
    qval,
    group_size) %>%
  dplyr::slice_head(n = 10) %>%
  print(width = Inf)

# Summarize identifiers and missing values
s1_sGene_identifier_summary <- data.frame(
  Total_Rows = nrow(s1_sGene),

  Unique_Splicing_Features = dplyr::n_distinct(
    s1_sGene$gene_id,
    na.rm = TRUE),

  Unique_Gene_Groups = dplyr::n_distinct(
    s1_sGene$group_id,
    na.rm = TRUE),

  Missing_Splicing_Feature_IDs = sum(
    is.na(s1_sGene$gene_id) |
      s1_sGene$gene_id == ""),

  Missing_Gene_Group_IDs = sum(
    is.na(s1_sGene$group_id) |
      s1_sGene$group_id == "")
)

s1_sGene_identifier_summary

# Examine the adjusted q-value distribution
summary(s1_sGene$qval)

# Count rows with q-value below 0.05
sum(s1_sGene$qval < 0.05,na.rm = TRUE)


# Inspect repeated Table S1 splicing features
# Count how many times each complete splicing-feature identifier occurs
s1_sGene_feature_counts <- s1_sGene %>%
  dplyr::count(
    gene_id,
    name = "Number_of_Rows") %>%    # gives the count column a clear name
  dplyr::arrange(
    dplyr::desc(
      Number_of_Rows)
  )

# Summarize the extent of feature duplication
s1_sGene_duplicate_summary <- s1_sGene_feature_counts %>%
  dplyr::filter(Number_of_Rows > 1) %>%

  dplyr::summarise(
    Number_of_Repeated_Feature_IDs = dplyr::n(),
    Total_Rows_From_Repeated_Features = sum(Number_of_Rows),

    Extra_Repeated_Rows = sum(Number_of_Rows - 1),

    Maximum_Rows_for_One_Feature = max(Number_of_Rows)
  )

s1_sGene_duplicate_summary

# Examine whether repeated feature IDs have different associations
s1_sGene_repeated_feature_details <- s1_sGene %>%
  dplyr::group_by(gene_id) %>%

  dplyr::filter(dplyr::n() > 1) %>%

  dplyr::summarise(
    Number_of_Rows = dplyr::n(),

    Number_of_Gene_Groups = dplyr::n_distinct(
      group_id),

    Number_of_Primary_Variants = dplyr::n_distinct(
      variant_id),

    Number_of_Distinct_Slopes = dplyr::n_distinct(
      slope),

    Minimum_Q_Value = min(
      qval,
      na.rm = TRUE),

    Maximum_Q_Value = max(
      qval,
      na.rm = TRUE),

    .groups = "drop") %>%

  dplyr::arrange(dplyr::desc(
      Number_of_Primary_Variants),
    dplyr::desc(Number_of_Rows))

head(s1_sGene_repeated_feature_details,20)

# Display repeated rows and their statistical results
s1_sGene %>%
  dplyr::group_by(gene_id) %>%

  dplyr::filter(dplyr::n() > 1) %>%

  dplyr::ungroup() %>%

  dplyr::arrange(
    gene_id,
    qval
  ) %>%

  dplyr::select(
    gene_id,
    group_id,
    variant_id,
    slope,
    pval_nominal,
    pval_perm,
    pval_beta,
    qval
  ) %>%

  dplyr::slice_head(n = 30) %>%
  print(width = Inf)

# Extract stable gene and splicing-feature information
s1_sGene_clean_all_rows <- s1_sGene %>%

  dplyr::transmute(

    # Preserve the complete original identifier
    original_splicing_feature_id = as.character(
      gene_id
    ),

    # Retain the coordinate and cluster portion
    #
    # Example:
    # 1:744003:744111:clu_7663_NA
    splicing_feature_coordinates = sub(
      ":ENSG.*$",
      "",
      as.character(gene_id)
    ),

    # Extract the stable Ensembl gene identifier
    # Example:
    # ENSG00000228327.3_7 becomes ENSG00000228327
    ensembl_gene_id = sub(
      "^(ENSG[0-9]+).*$",
      "\\1",
      as.character(group_id)
    ),

    # Extract the final gene-symbol portion from gene_id
    # Example:
    # ...:AL669831.1 becomes AL669831.1
    s1_gene_symbol = sub(
      "^.*:",
      "",
      as.character(gene_id)
    ),

    # Number of splicing features in the tested group
    group_size = as.numeric(
      group_size
    ),

    # Primary sQTL variant
    primary_variant = as.character(
      variant_id
    ),

    # Number of genetic variants tested
    number_of_variants_tested = as.numeric(
      num_var
    ),

    # Variant distance from the TSS
    distance_to_TSS_bp = as.numeric(
      tss_distance
    ),

    # Minor-allele information
    minor_allele_samples = as.numeric(
      ma_samples
    ),

    minor_allele_count = as.numeric(
      ma_count
    ),

    minor_allele_frequency = as.numeric(
      maf
    ),

    # Direction and size of the variant–splicing association
    sQTL_slope = as.numeric(
      slope
    ),

    sQTL_slope_standard_error = as.numeric(
      slope_se
    ),

    # Statistical results
    nominal_p_value = as.numeric(
      pval_nominal
    ),

    permutation_p_value = as.numeric(
      pval_perm
    ),

    beta_approximated_p_value = as.numeric(
      pval_beta
    ),

    q_value = as.numeric(
      qval
    ),

    nominal_p_value_threshold = as.numeric(
      pval_nominal_threshold
    )
  ) %>%

  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    !is.na(original_splicing_feature_id),
    original_splicing_feature_id != ""
  ) %>%

  dplyr::arrange(
    q_value,
    nominal_p_value
  )

# Inspect the cleaned splicing table
dim(s1_sGene_clean_all_rows)

head(s1_sGene_clean_all_rows,10)


# Examine repeated splicing-feature identifiers
# Count the number of rows representing each complete splicing feature
s1_sGene_feature_counts <- s1_sGene_clean_all_rows %>%
  dplyr::count(
    original_splicing_feature_id,
    name = "Number_of_Rows"
  ) %>%

  dplyr::arrange(
    dplyr::desc(
      Number_of_Rows
    )
  )


# Display the features occurring more than once
s1_sGene_repeated_feature_counts <- s1_sGene_feature_counts %>%
  dplyr::filter(
    Number_of_Rows > 1
  )

s1_sGene_repeated_feature_counts


# Summarize the repeated identifiers
s1_sGene_repeated_feature_summary <- data.frame(

  Total_sQTL_Rows = nrow(
    s1_sGene_clean_all_rows
  ),

  Unique_Splicing_Features = dplyr::n_distinct(
    s1_sGene_clean_all_rows$original_splicing_feature_id
  ),

  Repeated_Splicing_Features = nrow(
    s1_sGene_repeated_feature_counts
  ),

  Rows_Belonging_to_Repeated_Features = sum(
    s1_sGene_repeated_feature_counts$Number_of_Rows
  ),

  Extra_Rows_From_Repetition = sum(
    s1_sGene_repeated_feature_counts$Number_of_Rows - 1
  ),

  Maximum_Rows_Per_Feature = max(
    s1_sGene_feature_counts$Number_of_Rows
  )
)

s1_sGene_repeated_feature_summary

# Examine the details of repeated splicing features
s1_sGene_repeated_feature_details <- s1_sGene_clean_all_rows %>%

  dplyr::filter(
    original_splicing_feature_id %in%
      s1_sGene_repeated_feature_counts$original_splicing_feature_id
  ) %>%

  dplyr::group_by(
    original_splicing_feature_id,
    ensembl_gene_id,
    s1_gene_symbol
  ) %>%

  dplyr::summarise(

    Number_of_Rows = dplyr::n(),

    Number_of_Primary_Variants = dplyr::n_distinct(
      primary_variant
    ),

    Number_of_Distinct_Slopes = dplyr::n_distinct(
      sQTL_slope
    ),

    Smallest_Q_Value = min(
      q_value,
      na.rm = TRUE
    ),

    Largest_Q_Value = max(
      q_value,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  dplyr::arrange(
    dplyr::desc(
      Number_of_Rows
    ),
    Smallest_Q_Value
  )

s1_sGene_repeated_feature_details


# Confirm the structure of all repeated splicing features
# Count repeated features that have more than one primary variant
sum(s1_sGene_repeated_feature_details$Number_of_Primary_Variants > 1)

# Count repeated features that have more than one sQTL slope
sum(s1_sGene_repeated_feature_details$Number_of_Distinct_Slopes > 1)

# Display any repeated features with different variants or slopes
s1_sGene_repeated_feature_details %>%
  dplyr::filter(
    Number_of_Primary_Variants > 1 |
      Number_of_Distinct_Slopes > 1
  )

# Order rows from the smallest to the largest q-value
s1_sGene_unique_features <- s1_sGene_clean_all_rows %>%
  dplyr::arrange(q_value,nominal_p_value) %>%

  # Retain the first row for each splicing-feature identifier
  # Because the table was ordered by q-value, this keeps the
  # record with the smallest q-value
  dplyr::distinct(
    original_splicing_feature_id,
    .keep_all = TRUE
  )


# Confirm the resulting dimensions
dim(s1_sGene_unique_features)

# Compare the number of rows with the number of unique feature IDs
data.frame(
  Number_of_Rows = nrow(
    s1_sGene_unique_features
  ),

  Number_of_Unique_Splicing_Features = dplyr::n_distinct(
    s1_sGene_unique_features$original_splicing_feature_id
  ),

  Remaining_Duplicate_Feature_IDs = sum(
    duplicated(
      s1_sGene_unique_features$original_splicing_feature_id
    )
  )
)

# Create one row per significant sGene
# For comparison with the TCF4 peak-associated gene sets, we need a gene-level table
# Collapse significant splicing features to their parent genes
s1_sGene_gene_summary <- s1_sGene_unique_features %>%

  dplyr::group_by(
    ensembl_gene_id
  ) %>%

  dplyr::summarise(

    # Retain the gene symbol
    s1_gene_symbol = dplyr::first(
      s1_gene_symbol
    ),

    # Count the number of different significant splicing features
    # assigned to this gene
    Number_of_Significant_Splicing_Features = dplyr::n_distinct(
      original_splicing_feature_id
    ),

    # Record the feature with the smallest q-value
    Strongest_Splicing_Feature_ID = dplyr::first(
      original_splicing_feature_id
    ),

    # Record the primary variant associated with that feature
    Strongest_Primary_Variant = dplyr::first(
      primary_variant
    ),

    # Record the slope for that feature–variant association
    Strongest_sQTL_Slope = dplyr::first(
      sQTL_slope
    ),

    # Record the smallest q-value found for the gene
    Minimum_sQTL_Q_Value = min(
      q_value,
      na.rm = TRUE
    ),

    # Every gene in this table has at least one significant
    # splicing-QTL association
    Significant_sGene = TRUE,

    .groups = "drop"
  )

# Inspect the gene-level sGene table
dim(s1_sGene_gene_summary)

head(s1_sGene_gene_summary,10)

# Confirm that every Ensembl gene appears only once
sum(duplicated(s1_sGene_gene_summary$ensembl_gene_id))


# Compare TCF4-associated genes with fetal-brain sGenes
# Count the overlap for each complete TCF4 gene set
s1_sGene_overlap_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_TCF4_Peak_Associated_Genes = c(
    nrow(npc_all_genes),
    nrow(mcclay_all_genes),
    nrow(forrest_all_genes)
  ),

  Overlap_With_Significant_sGenes = c(

    sum(
      npc_all_genes$ensembl_gene_id %in%
        s1_sGene_gene_summary$ensembl_gene_id
    ),

    sum(
      mcclay_all_genes$ensembl_gene_id %in%
        s1_sGene_gene_summary$ensembl_gene_id
    ),

    sum(
      forrest_all_genes$ensembl_gene_id %in%
        s1_sGene_gene_summary$ensembl_gene_id
    )
  )
)


# Calculate the percentage of each TCF4 gene set that is an sGene
s1_sGene_overlap_summary <- s1_sGene_overlap_summary %>%

  dplyr::mutate(

    Percent_of_TCF4_Set_in_sGene_List = round(
      100 *
        Overlap_With_Significant_sGenes /
        Total_TCF4_Peak_Associated_Genes,
      digits = 2
    )
  )


# Display the result
s1_sGene_overlap_summary

# Join the fetal-brain sGene information to the TCF4 membership table
tcf4_membership_sGene <- tcf4_gene_membership %>%

  dplyr::left_join(
    s1_sGene_gene_summary,
    by = "ensembl_gene_id"
  ) %>%

  dplyr::mutate(

    # Genes absent from the significant sGene table receive FALSE
    Significant_sGene = tidyr::replace_na(
      Significant_sGene,
      FALSE
    ),

    # Genes without a significant splicing feature receive zero
    Number_of_Significant_Splicing_Features =
      tidyr::replace_na(
        Number_of_Significant_Splicing_Features,
        0L
      )
  )

# Calculate sGene representation within each TCF4 membership group
tcf4_membership_sGene_summary <- tcf4_membership_sGene %>%

  dplyr::group_by(
    Membership_Group
  ) %>%

  dplyr::summarise(

    # Total number of TCF4-associated genes in the group
    Total_TCF4_Genes = dplyr::n(),

    # Number that are significant fetal-brain sGenes
    Significant_sGenes = sum(
      Significant_sGene
    ),

    # Percentage of TCF4-associated genes that are sGenes
    Percent_Significant_sGenes = round(
      100 *
        Significant_sGenes /
        Total_TCF4_Genes,
      digits = 2
    ),

    # Total number of significant splicing features represented
    # among all genes in this membership group
    Total_Significant_Splicing_Features = sum(
      Number_of_Significant_Splicing_Features
    ),

    # Average number of significant splicing features per sGene
    Mean_Splicing_Features_Per_sGene = round(
      mean(
        Number_of_Significant_Splicing_Features[
          Significant_sGene
        ]
      ),
      digits = 2
    ),

    # Largest number of significant features found in one gene
    Maximum_Splicing_Features_in_One_Gene = max(
      Number_of_Significant_Splicing_Features
    ),

    .groups = "drop"
  )


# Display the summary
tcf4_membership_sGene_summary


# Prrmutation analysis-TablsS1-sGene #
# Within-TCF4 permutation analysis for significant fetal-brain sGenes
# Set the order of the seven TCF4 membership groups
# This ensures that tables and permutation results use the same order
tcf4_membership_sGene <- tcf4_membership_sGene %>%
  dplyr::mutate(
    Membership_Group = factor(
      Membership_Group,
      levels = c(
        "NPC only",
        "McClay only",
        "Forrest only",
        "NPC + McClay",
        "NPC + Forrest",
        "McClay + Forrest",
        "All three"
      )
    )
  )


# Check the information that will enter the permutation analysis
data.frame(
  Total_TCF4_Genes = nrow(
    tcf4_membership_sGene
  ),

  Total_Significant_sGenes = sum(
    tcf4_membership_sGene$Significant_sGene
  ),

  Total_Not_Significant_sGenes = sum(
    !tcf4_membership_sGene$Significant_sGene
  )
)

# Create one membership column for each of the seven groups
# Each row represents one TCF4-associated gene.
# Each column represents one membership group.
# A value of 1 means that the gene belongs to that group.
# A value of 0 means that it does not belong to that group.
# The "- 1" removes the intercept column so that all seven
# membership groups receive their own separate columns.
sGene_membership_matrix <- model.matrix(
  ~ Membership_Group - 1,
  data = tcf4_membership_sGene)


# Replace the automatically generated column names
# with the simpler membership-group names
colnames(sGene_membership_matrix) <- levels(
  tcf4_membership_sGene$Membership_Group)


# Inspect the dimensions of the matrix
dim(sGene_membership_matrix)

# Display its first ten rows
head(sGene_membership_matrix,10)

# Convert the logical sGene status into numeric values
# TRUE becomes 1
# FALSE becomes 0
observed_sGene_status <- as.integer(
  tcf4_membership_sGene$Significant_sGene)


# Count the actual number of significant sGenes
# observed in every membership group
observed_sGene_counts <- colSums(
  sGene_membership_matrix *
    observed_sGene_status
)


# Display the observed counts
observed_sGene_counts

# Set the seed so that the random analysis can be reproduced
# Running the code again with the same seed will give the same results
set.seed(12345)


# Define the number of random permutations
number_of_permutations <- 10000


# In every permutation:
# 1. Randomly rearrange the sGene labels across the TCF4 genes.
# 2. Preserve the total number of significant sGenes.
# 3. Preserve the seven membership-group sizes.
# 4. Count the shuffled sGenes in each membership group.
# replicate() repeats these steps 10,000 times.
permuted_sGene_counts <- replicate(
  number_of_permutations,
  {

    # Randomly rearrange the observed 1 and 0 labels
    shuffled_sGene_status <- sample(
      observed_sGene_status,
      replace = FALSE
    )


    # Count shuffled sGenes in each membership group
    colSums(
      sGene_membership_matrix *
        shuffled_sGene_status
    )
  }
)

# Confirm the dimensions of the permutation result
dim(permuted_sGene_counts)

# Calculate the average random sGene count for every group
# mean_random_sGene_counts tells us how many sGenes are normally expected in
# each group if sGene status has no relationship with TCF4 membership.
mean_random_sGene_counts <- rowMeans(permuted_sGene_counts)


# Calculate the standard deviation of the random counts
# apply(..., MARGIN = 1, ...) performs the calculation
# separately across every row of the permutation matrix.
random_sGene_standard_deviation <- apply(
  permuted_sGene_counts,
  MARGIN = 1,
  FUN = stats::sd
)

# Calculate empirical P values
# Compare every permuted count with the corresponding observed count
# ">=" asks whether a random permutation produced at least as many
# sGenes as the real data
# This is a one-sided enrichment test
permutation_at_least_as_large_as_observed <- sweep(
  permuted_sGene_counts,
  MARGIN = 1,
  STATS = observed_sGene_counts,
  FUN = ">="
)


# Calculate the empirical P value
# The added 1 prevents an empirical P value of exactly zero.
empirical_sGene_p_values <- (
  rowSums(permutation_at_least_as_large_as_observed) + 1) / (number_of_permutations + 1)

# Create the complete permutation result table
# Combine observed counts, random expectations and empirical P values
sGene_membership_permutation_results <- data.frame(

  Membership_Group = names(
    observed_sGene_counts
  ),

  Total_Genes_in_Group = colSums(
    sGene_membership_matrix
  ),

  Observed_sGenes = as.numeric(
    observed_sGene_counts
  ),

  Mean_Random_sGenes = as.numeric(
    mean_random_sGene_counts
  ),

  Random_Standard_Deviation = as.numeric(
    random_sGene_standard_deviation
  ),

  Empirical_P_Value = as.numeric(
    empirical_sGene_p_values
  )
)

# Add effect-size measurements and FDR
# Calculate the difference between observed and random counts,
# the observed-to-random ratio and multiple-testing-adjusted FDR
sGene_membership_permutation_results <-
  sGene_membership_permutation_results %>%

  dplyr::mutate(

    # Positive values indicate more sGenes than expected.
    # Negative values indicate fewer sGenes than expected.
    Observed_Minus_Random = round(
      Observed_sGenes - Mean_Random_sGenes,
      digits = 2
    ),

    # A ratio above 1 indicates more sGenes than expected.
    # A ratio below 1 indicates fewer sGenes than expected.
    Observed_to_Random_Ratio = round(
      Observed_sGenes / Mean_Random_sGenes,
      digits = 3
    ),

    # Correct the seven empirical P values using
    # the Benjamini-Hochberg method.
    Empirical_FDR = p.adjust(
      Empirical_P_Value,
      method = "BH"
    ),

    # TRUE means the enrichment survives FDR correction.
    Significant_After_FDR = Empirical_FDR < 0.05
  )

# Round display columns and order the results
# Round values for easier reading and arrange from the
# strongest to the weakest enrichment evidence
sGene_membership_permutation_results <-
  sGene_membership_permutation_results %>%

  dplyr::mutate(

    Mean_Random_sGenes = round(
      Mean_Random_sGenes,
      digits = 2
    ),

    Random_Standard_Deviation = round(
      Random_Standard_Deviation,
      digits = 2
    )
  ) %>%

  dplyr::arrange(
    Empirical_FDR,
    Empirical_P_Value
  )


# Display the final results
sGene_membership_permutation_results


# Joint Tri1-Tri2 eGene permutation analysis
# Confirm the five TCF4 membership groups are in the desired order
tcf4_selected_trimester_eGene <-
  tcf4_selected_trimester_eGene %>%
  dplyr::mutate(
    Selected_Membership_Group = factor(
      Selected_Membership_Group,

      levels = c(
        "NPC only",
        "McClay only",
        "Forrest only",
        "McClay + Forrest",
        "All three"
      )
    )
  )

# Create the membership matrix
# Each row represents one TCF4-associated gene
# Each column represents one selected membership group.
# 1 = member in that group/ 0 = non-menber
# The "- 1" removes the intercept so that all five groups
trimester_membership_matrix <- model.matrix(
  ~ Selected_Membership_Group - 1,
  data = tcf4_selected_trimester_eGene
)

# Replace the automatically generated column names
colnames(trimester_membership_matrix) <- levels(
  tcf4_selected_trimester_eGene$
    Selected_Membership_Group
)

# Check the dimensions
dim(trimester_membership_matrix)

# Every row should sum to exactly 1
table(rowSums(trimester_membership_matrix))

# Create the paired Tri1–Tri2 status matrix
# Convert the logical detection labels to numeric values
# TRUE = 1/ FALSE = 0
trimester_eGene_status_matrix <- cbind(
  Tri1 = as.integer(tcf4_selected_trimester_eGene$Detected_Tri1_eGene),

  Tri2 = as.integer(tcf4_selected_trimester_eGene$Detected_Tri2_eGene))


# Display the first ten genes
head(trimester_eGene_status_matrix,10)

# Create the five permutation outcomes
trimester_eGene_outcome_matrix <- cbind(

  # Detected in the Tri1 significant eGene list
  Tri1_eGene =trimester_eGene_status_matrix[, "Tri1"],

  # Detected in the Tri2 significant eGene list
  Tri2_eGene =trimester_eGene_status_matrix[, "Tri2"],

  # Detected in both separate trimester analyses
  Both_Trimesters = as.integer(
    trimester_eGene_status_matrix[, "Tri1"] == 1 &
      trimester_eGene_status_matrix[, "Tri2"] == 1
  ),

  # Detected in Tri1 but not Tri2
  Tri1_Only = as.integer(
    trimester_eGene_status_matrix[, "Tri1"] == 1 &
      trimester_eGene_status_matrix[, "Tri2"] == 0
  ),

  # Detected in Tri2 but not Tri1
  Tri2_Only = as.integer(
    trimester_eGene_status_matrix[, "Tri1"] == 0 &
      trimester_eGene_status_matrix[, "Tri2"] == 1
  )
)

# Verify the outcome totals
colSums(trimester_eGene_outcome_matrix)

# Calculate the observed counts
# Count how many genes in each membership group are positive
# for each of the five trimester outcomes
observed_trimester_eGene_counts <- crossprod(
  trimester_membership_matrix,
  trimester_eGene_outcome_matrix)


observed_trimester_eGene_counts

# Run 10,000 joint permutations
# Set a random seed so the results can be reproduced
set.seed(12345)


# Define the number of permutations
number_of_permutations <- 10000


# In each permutation:
# 1. Randomly rearrange the gene rows.
# 2. Keep each gene's paired Tri1-Tri2 status together.
# 3. Recreate the five outcomes.
# 4. Count randomized positive labels in each TCF4 group.
permuted_trimester_eGene_counts <- replicate(
  number_of_permutations,

  {

    # Generate a random ordering of all retained genes

    shuffled_gene_rows <- sample(
      seq_len(
        nrow(
          trimester_eGene_status_matrix
        )
      ),
      replace = FALSE
    )


    # Move each gene's complete Tri1-Tri2 pattern together

    shuffled_trimester_status <-
      trimester_eGene_status_matrix[
        shuffled_gene_rows,
        ,
        drop = FALSE
      ]


    # Recreate the five outcomes from the shuffled status

    shuffled_trimester_outcomes <- cbind(

      Tri1_eGene =
        shuffled_trimester_status[, "Tri1"],

      Tri2_eGene =
        shuffled_trimester_status[, "Tri2"],

      Both_Trimesters = as.integer(
        shuffled_trimester_status[, "Tri1"] == 1 &
          shuffled_trimester_status[, "Tri2"] == 1
      ),

      Tri1_Only = as.integer(
        shuffled_trimester_status[, "Tri1"] == 1 &
          shuffled_trimester_status[, "Tri2"] == 0
      ),

      Tri2_Only = as.integer(
        shuffled_trimester_status[, "Tri1"] == 0 &
          shuffled_trimester_status[, "Tri2"] == 1
      )
    )


    # Count positive labels in each TCF4 group

    crossprod(
      trimester_membership_matrix,
      shuffled_trimester_outcomes
    )
  }
)

dim(permuted_trimester_eGene_counts) # Check the permutation-array dimensions

# Calculate random expectations
# Calculate the average randomized count for every comparison
mean_random_trimester_eGene_counts <- apply(
  permuted_trimester_eGene_counts,
  MARGIN = c(1, 2),
  FUN = mean
)


# Calculate the standard deviation of the randomized counts
random_trimester_eGene_standard_deviation <- apply(
  permuted_trimester_eGene_counts,
  MARGIN = c(1, 2),
  FUN = stats::sd
)

# Calculate empirical P values
# Repeat the observed-count matrix across all permutations
observed_trimester_eGene_count_array <- array(
  observed_trimester_eGene_counts,
  dim = dim(
    permuted_trimester_eGene_counts
  )
)


# Test whether each randomized count was at least as large
# as the corresponding observed count.
# This is a one-sided enrichment test.
random_at_least_as_large_as_observed <-
  permuted_trimester_eGene_counts >=
  observed_trimester_eGene_count_array


# Count permutations at least as large as observed
number_of_random_results_at_least_as_large <- apply(
  random_at_least_as_large_as_observed,
  MARGIN = c(1, 2),
  FUN = sum
)


# Calculate empirical P values using the added-one correction
empirical_trimester_eGene_p_values <- (
  number_of_random_results_at_least_as_large + 1) / (number_of_permutations + 1)

# Create the result table
# Create all 25 combinations:
# 5 membership groups × 5 trimester outcomes
trimester_eGene_permutation_results <- expand.grid(
  Membership_Group = rownames(observed_trimester_eGene_counts),
  Trimester_Outcome = colnames(observed_trimester_eGene_counts),
  stringsAsFactors = FALSE)


# Add observed and randomized statistics
trimester_eGene_permutation_results <-
  trimester_eGene_permutation_results %>%
  dplyr::mutate(

    Total_Genes_in_Group = rep(
      colSums(trimester_membership_matrix),
      times = ncol(observed_trimester_eGene_counts)),

    Observed_eGenes = as.numeric(observed_trimester_eGene_counts),

    Mean_Random_eGenes = as.numeric(mean_random_trimester_eGene_counts),

    Random_Standard_Deviation = as.numeric(random_trimester_eGene_standard_deviation),

    Empirical_P_Value = as.numeric(empirical_trimester_eGene_p_values))

# Add enrichment effect sizes
trimester_eGene_permutation_results <-
  trimester_eGene_permutation_results %>%

  dplyr::mutate(

    # Positive values indicate more genes than expected
    Observed_Minus_Random = round(
      Observed_eGenes -
        Mean_Random_eGenes,
      digits = 2
    ),

    # Values above 1 indicate enrichment
    Observed_to_Random_Ratio = round(
      Observed_eGenes /
        Mean_Random_eGenes,
      digits = 3
    )
  )

# Calculate FDR within each outcome
# Correct the five membership tests separately
# within each trimester outcome
trimester_eGene_permutation_results <-
  trimester_eGene_permutation_results %>%
  dplyr::group_by(Trimester_Outcome) %>%

  dplyr::mutate(Empirical_FDR_Within_Outcome = p.adjust(
      Empirical_P_Value,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

# Calculate FDR across all 25 tests
# Correct all membership groups and outcomes together.
# Use this global FDR for the primary interpretation.
trimester_eGene_permutation_results <-
  trimester_eGene_permutation_results %>%

  dplyr::mutate(

    Empirical_FDR_Across_All_Tests = p.adjust(
      Empirical_P_Value,
      method = "BH"
    ),

    Significant_After_Global_FDR =
      Empirical_FDR_Across_All_Tests < 0.05,

    Mean_Random_eGenes = round(
      Mean_Random_eGenes,
      digits = 2
    ),

    Random_Standard_Deviation = round(
      Random_Standard_Deviation,
      digits = 2
    )
  ) %>%

  dplyr::arrange(
    Empirical_FDR_Across_All_Tests,
    Empirical_P_Value
  ) %>%

  tibble::as_tibble()

# Display all 25 results
trimester_eGene_permutation_results %>%
  print(n = Inf,width = Inf)
