# TCF4 x developing-brain cross-ancestry atlas integration
# Paper: Cross-ancestry atlas of gene, isoform, and splicing regulation in the
# developing human brain (Science 2024; eadh0829)
# https://pubmed.ncbi.nlm.nih.gov/38781368/
# The paper analyzed RNA-seq on GRCh37 with GENCODE v29lift37.
# This paper has 7 supplementary tables
# THIS CODE IS FOR TABLES2

# Load libraries
library(readxl)  # imports Excel workbooks into R
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(ggplot2)

# Import hg19 data annotated peaks
npc_peak_anno <- read.csv("npc_hg19_peak_annotation.csv")
mcclay_peak_anno <- read.csv("mcclay_hg19_peak_annotation.csv")
forrest_peak_anno <- read.csv("forrest_hg19_peak_annotation.csv")

# Import ancestry-specific fetal-brain eGene results from Table S2
table_s2_file <- file.choose()

# Display all sheet names in the workbook
# This confirms the exact names before importing them
readxl::excel_sheets(table_s2_file)

# Import significant European-ancestry eGenes
s2_EUR_eGene <- readxl::read_excel(
  path = table_s2_file,
  sheet = "ST2-1-EUR-eGenes"
)


# Import significant Admixed-American-ancestry eGenes
s2_AMR_eGene <- readxl::read_excel(
  path = table_s2_file,
  sheet = "ST2-2-AMR-eGenes"
)


# Import significant African-ancestry eGenes
s2_AFR_eGene <- readxl::read_excel(
  path = table_s2_file,
  sheet = "ST2-3-AFR-eGenes"
)

# Check the number of rows and columns
dim(s2_EUR_eGene)
dim(s2_AMR_eGene)
dim(s2_AFR_eGene)


# Display all column names
colnames(s2_EUR_eGene)
colnames(s2_AMR_eGene)
colnames(s2_AFR_eGene)


# Display the first six rows
head(s2_EUR_eGene)
head(s2_AMR_eGene)
head(s2_AFR_eGene)

# Create a clearly labeled EUR eGene table
s2_EUR_eGene_clean <- s2_EUR_eGene %>%
  dplyr::transmute(
    # Remove any Ensembl version suffix if one is present
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      pid
    ),

    ancestry = "EUR",
    primary_variant = sid,
    distance_to_TSS_bp = dist,
    eQTL_slope = slope,
    nominal_p_value = npval,
    permutation_p_value = ppval,
    beta_approximated_p_value = bpval,
    q_value = qval,
    nominal_p_value_threshold = pval_nominal_threshold,
    number_of_variants_tested = nvar) %>%

  # Put the strongest result first if an identifier is repeated
  dplyr::arrange(
    q_value,
    nominal_p_value) %>%

  # Retain one row per stable Ensembl gene
  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

# Create a clearly labeled AMR eGene table
s2_AMR_eGene_clean <- s2_AMR_eGene %>%
  dplyr::transmute(
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      pid
    ),
    ancestry = "AMR",
    primary_variant = sid,
    distance_to_TSS_bp = dist,
    eQTL_slope = slope,
    nominal_p_value = npval,
    permutation_p_value = ppval,
    beta_approximated_p_value = bpval,
    q_value = qval,
    nominal_p_value_threshold = pval_nominal_threshold,
    number_of_variants_tested = nvar ) %>%

  dplyr::arrange(
    q_value,
    nominal_p_value ) %>%

  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )


# Create a clearly labeled AFR eGene table
s2_AFR_eGene_clean <- s2_AFR_eGene %>%
  dplyr::transmute(
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      pid
    ),

    ancestry = "AFR",
    primary_variant = sid,
    distance_to_TSS_bp = dist,
    eQTL_slope = slope,
    nominal_p_value = npval,
    permutation_p_value = ppval,
    beta_approximated_p_value = bpval,
    q_value = qval,
    nominal_p_value_threshold = pval_nominal_threshold,
    number_of_variants_tested = nvar) %>%

  dplyr::arrange(
    q_value,
    nominal_p_value) %>%

  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

# Summarize the number of unique eGenes in each ancestry
s2_eGene_cleaning_summary <- data.frame(
  Ancestry = c(
    "EUR",
    "AMR",
    "AFR"
  ),

  Original_Rows = c(
    nrow(s2_EUR_eGene),
    nrow(s2_AMR_eGene),
    nrow(s2_AFR_eGene)
  ),

  Unique_eGenes = c(
    nrow(s2_EUR_eGene_clean),
    nrow(s2_AMR_eGene_clean),
    nrow(s2_AFR_eGene_clean)
  ),

  Remaining_Duplicate_Gene_IDs = c(

    sum(
      duplicated(
        s2_EUR_eGene_clean$ensembl_gene_id
      )
    ),

    sum(
      duplicated(
        s2_AMR_eGene_clean$ensembl_gene_id
      )
    ),

    sum(
      duplicated(
        s2_AFR_eGene_clean$ensembl_gene_id
      )
    )
  ),

  Missing_Gene_IDs = c(

    sum(
      is.na(s2_EUR_eGene_clean$ensembl_gene_id) |
        s2_EUR_eGene_clean$ensembl_gene_id == ""
    ),

    sum(
      is.na(s2_AMR_eGene_clean$ensembl_gene_id) |
        s2_AMR_eGene_clean$ensembl_gene_id == ""
    ),

    sum(
      is.na(s2_AFR_eGene_clean$ensembl_gene_id) |
        s2_AFR_eGene_clean$ensembl_gene_id == ""
    )
  )
)


s2_eGene_cleaning_summary


#create TCF4 peak-associated gene sets
# Create the NPC gene list
# ENSEMBL contains the Ensembl gene identifier assigned to each peak.
# SYMBOL contains the corresponding gene symbol
npc_all_genes <- npc_peak_anno %>%
  dplyr::transmute(
    # Remove an Ensembl version suffix if one is present
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      trimws(
        as.character(
          ENSEMBL
        )
      )
    ),

    gene_symbol = trimws(
      as.character(
        SYMBOL
      )
    )
  ) %>%

  # Remove rows without a valid Ensembl gene identifier
  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    grepl(
      "^ENSG[0-9]+$",
      ensembl_gene_id
    )
  ) %>%

  # A gene may have multiple TCF4 peaks.
  # Keep only one row per Ensembl gene.
  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )


# Create the McClay gene list
mcclay_all_genes <- mcclay_peak_anno %>%

  dplyr::transmute(

    ensembl_gene_id = sub(
      "\\..*$",
      "",
      trimws(
        as.character(
          ENSEMBL
        )
      )
    ),

    gene_symbol = trimws(
      as.character(
        SYMBOL
      )
    )
  ) %>%

  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    grepl(
      "^ENSG[0-9]+$",
      ensembl_gene_id
    )
  ) %>%

  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

# Create the Forrest gene list
forrest_all_genes <- forrest_peak_anno %>%
  dplyr::transmute(
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      trimws(
        as.character(
          ENSEMBL
        )
      )
    ),

    gene_symbol = trimws(
      as.character(
        SYMBOL
      )
    )
  ) %>%

  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    grepl(
      "^ENSG[0-9]+$",
      ensembl_gene_id
    )
  ) %>%

  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )


# Display the first six genes from each dataset
head(npc_all_genes)
head(mcclay_all_genes)
head(forrest_all_genes)

# Count unique peak-associated genes
tcf4_gene_set_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),
  Unique_Peak_Associated_Genes = c(
    nrow(npc_all_genes),
    nrow(mcclay_all_genes),
    nrow(forrest_all_genes)
  )
)


tcf4_gene_set_summary

# Create an NPC membership table
npc_membership <- npc_all_genes %>%
  dplyr::transmute(ensembl_gene_id,Present_NPC = TRUE)


# Create a McClay membership table
mcclay_membership <- mcclay_all_genes %>%
  dplyr::transmute(ensembl_gene_id,Present_McClay = TRUE)


# Create a Forrest membership table
forrest_membership <- forrest_all_genes %>%
  dplyr::transmute(ensembl_gene_id, Present_Forrest = TRUE)

# Each table now contains one Ensembl gene ID and a TRUE value indicating that it belongs to that TCF4 dataset
# Now, Combine all TCF4 genes
# full_join() retains every gene found in any of the three datasets
tcf4_gene_membership <- npc_membership %>%
  dplyr::full_join(
    mcclay_membership,
    by = "ensembl_gene_id"
  ) %>%

  dplyr::full_join(
    forrest_membership,
    by = "ensembl_gene_id"
  )

# Convert missing membership values to FALSE
# Genes absent from a dataset currently have NA. Convert those missing values to FALSE
tcf4_gene_membership <- tcf4_gene_membership %>%
  dplyr::mutate(Present_NPC = tidyr::replace_na(Present_NPC,FALSE),

    Present_McClay = tidyr::replace_na(Present_McClay,FALSE),

    Present_Forrest = tidyr::replace_na(Present_Forrest,FALSE)
  )

# Assign the seven mutually exclusive membership groups
tcf4_gene_membership <- tcf4_gene_membership %>%

  dplyr::mutate(

    Membership_Group = dplyr::case_when(

      Present_NPC &
        !Present_McClay &
        !Present_Forrest ~
        "NPC only",

      !Present_NPC &
        Present_McClay &
        !Present_Forrest ~
        "McClay only",

      !Present_NPC &
        !Present_McClay &
        Present_Forrest ~
        "Forrest only",

      Present_NPC &
        Present_McClay &
        !Present_Forrest ~
        "NPC + McClay",

      Present_NPC &
        !Present_McClay &
        Present_Forrest ~
        "NPC + Forrest",

      !Present_NPC &
        Present_McClay &
        Present_Forrest ~
        "McClay + Forrest",

      Present_NPC &
        Present_McClay &
        Present_Forrest ~
        "All three"
    )
  )

# Add gene symbols
# Combine the gene labels from the three datasets
tcf4_gene_labels <- dplyr::bind_rows(
  npc_all_genes,
  mcclay_all_genes,
  forrest_all_genes
) %>%

  # Retain one label per Ensembl gene
  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )


# Add the gene symbol to the membership table
tcf4_gene_membership <- tcf4_gene_membership %>%

  dplyr::left_join(
    tcf4_gene_labels,
    by = "ensembl_gene_id"
  )

# Set the membership-group order
tcf4_gene_membership <- tcf4_gene_membership %>%
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

# Verify the reconstructed groups
tcf4_membership_summary <- tcf4_gene_membership %>%
  dplyr::group_by(
    Membership_Group
  ) %>%
  dplyr::summarise(
    Number_of_Genes = dplyr::n(),
    .groups = "drop"
  )

tcf4_membership_summary



# Compare TCF4-associated genes with ancestry-specific fetal-brain eGenes
# Calculate overlap with NPC gene set
s2_NPC_eGene_overlap <- data.frame(
  TCF4_Dataset = "NPC",
  Ancestry = c(
    "EUR",
    "AMR",
    "AFR"
  ),

  Total_TCF4_Genes = nrow(
    npc_all_genes
  ),

  Significant_eGene_Overlap = c(

    sum(npc_all_genes$ensembl_gene_id %in%
        s2_EUR_eGene_clean$ensembl_gene_id
    ),

    sum(npc_all_genes$ensembl_gene_id %in%
        s2_AMR_eGene_clean$ensembl_gene_id
    ),

    sum(npc_all_genes$ensembl_gene_id %in%
        s2_AFR_eGene_clean$ensembl_gene_id
    )
  )
)

# Calculate overlap with McClay gene set
s2_McClay_eGene_overlap <- data.frame(
  TCF4_Dataset = "McClay",
  Ancestry = c(
    "EUR",
    "AMR",
    "AFR"
  ),

  Total_TCF4_Genes = nrow(
    mcclay_all_genes
  ),

  Significant_eGene_Overlap = c(

    sum(mcclay_all_genes$ensembl_gene_id %in%
        s2_EUR_eGene_clean$ensembl_gene_id
    ),

    sum(mcclay_all_genes$ensembl_gene_id %in%
        s2_AMR_eGene_clean$ensembl_gene_id
    ),

    sum(mcclay_all_genes$ensembl_gene_id %in%
        s2_AFR_eGene_clean$ensembl_gene_id
    )
  )
)

# Calculate overlap with Forrest gene set
s2_Forrest_eGene_overlap <- data.frame(
  TCF4_Dataset = "Forrest",
  Ancestry = c(
    "EUR",
    "AMR",
    "AFR"
  ),

  Total_TCF4_Genes = nrow(
    forrest_all_genes
  ),

  Significant_eGene_Overlap = c(

    sum(forrest_all_genes$ensembl_gene_id %in%
        s2_EUR_eGene_clean$ensembl_gene_id
    ),

    sum(forrest_all_genes$ensembl_gene_id %in%
        s2_AMR_eGene_clean$ensembl_gene_id
    ),

    sum(forrest_all_genes$ensembl_gene_id %in%
        s2_AFR_eGene_clean$ensembl_gene_id
    )
  )
)

# Combine and calculate percentages
# Combine the nine comparisons into one table
s2_eGene_overlap_summary <- dplyr::bind_rows(
  s2_NPC_eGene_overlap,
  s2_McClay_eGene_overlap,
  s2_Forrest_eGene_overlap
) %>%

  dplyr::mutate(

    Percent_of_TCF4_Set_in_eGene_List = round(
      100 *
        Significant_eGene_Overlap /
        Total_TCF4_Genes,
      digits = 2
    )
  )


# Display the overall overlap results
s2_eGene_overlap_summary

# Add ancestry-specific eGene labels to all TCF4 genes
# For every gene in the combined TCF4 universe, determine whether
# it is a significant eGene in EUR, AMR and AFR
tcf4_membership_ancestry_eGene <- tcf4_gene_membership %>%
  dplyr::mutate(

    Significant_eGene_EUR =ensembl_gene_id %in%
      s2_EUR_eGene_clean$ensembl_gene_id,

    Significant_eGene_AMR =ensembl_gene_id %in%
      s2_AMR_eGene_clean$ensembl_gene_id,

    Significant_eGene_AFR =ensembl_gene_id %in%
      s2_AFR_eGene_clean$ensembl_gene_id
  )

# Count the number of ancestry panels detecting each gene
tcf4_membership_ancestry_eGene <-
  tcf4_membership_ancestry_eGene %>%
  dplyr::mutate(

    # Add the three TRUE/FALSE values
    # TRUE is treated as 1
    # FALSE is treated as 0
    # The result ranges from 0 to 3
    Number_of_Ancestries_With_eGene = rowSums(
      cbind(Significant_eGene_EUR,
            Significant_eGene_AMR,
            Significant_eGene_AFR
      )
    ),

    # Identify genes detected in at least two ancestry panels
    eGene_in_at_least_2_Ancestries = Number_of_Ancestries_With_eGene >= 2,

    # Identify genes detected in all three ancestry panels
    eGene_in_all_3_Ancestries = Number_of_Ancestries_With_eGene == 3
  )


# Summarize the seven TCF4 membership groups
# Calculate ancestry-specific eGene representation in every
# mutually exclusive TCF4 membership group
# Display all rows and columns without truncation
s2_membership_eGene_summary %>%

  dplyr::select(
    Membership_Group,
    Total_TCF4_Genes,

    EUR_eGenes,
    Percent_EUR_eGenes,

    AMR_eGenes,
    Percent_AMR_eGenes,

    AFR_eGenes,
    Percent_AFR_eGenes,

    eGenes_in_at_least_2_Ancestries,
    Percent_in_at_least_2_Ancestries,

    eGenes_in_all_3_Ancestries,
    Percent_in_all_3_Ancestries
  ) %>%

  print(
    n = Inf,
    width = Inf
  )

#########################
# Joint cross-ancestry eGene permutation analysis

# Set the order of the seven TCF4 membership groups

tcf4_membership_ancestry_eGene <-
  tcf4_membership_ancestry_eGene %>%

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


# Create the membership matrix
# Each row represents one TCF4-associated gene.
# Each column represents one TCF4 membership group.
# The "-1" removes the intercept so that all seven groups
# receive their own separate columns.
ancestry_eGene_membership_matrix <- model.matrix(
  ~ Membership_Group - 1,
  data = tcf4_membership_ancestry_eGene
)


# Replace the automatically generated column names
colnames(ancestry_eGene_membership_matrix) <- levels(
  tcf4_membership_ancestry_eGene$Membership_Group)

# Check its dimensions
dim(ancestry_eGene_membership_matrix)

# Create the three ancestry-status columns
# Convert TRUE and FALSE ancestry labels to 1 and 0
ancestry_eGene_status_matrix <- cbind(
  EUR = as.integer(tcf4_membership_ancestry_eGene$Significant_eGene_EUR),

  AMR = as.integer(tcf4_membership_ancestry_eGene$Significant_eGene_AMR),

  AFR = as.integer(tcf4_membership_ancestry_eGene$Significant_eGene_AFR)
)


# Inspect the first ten genes
# each row is a TCF4 gene
head(ancestry_eGene_status_matrix,10)

# Create the five outcomes
# Create the three individual ancestry outcomes and the
# two cross-ancestry detection outcomes
ancestry_eGene_outcome_matrix <- cbind(

  EUR_eGene =ancestry_eGene_status_matrix[, "EUR"],

  AMR_eGene =ancestry_eGene_status_matrix[, "AMR"],

  AFR_eGene =ancestry_eGene_status_matrix[, "AFR"],

  At_least_2_Ancestries =
    as.integer(
      rowSums(
        ancestry_eGene_status_matrix
      ) >= 2
    ),

  All_3_Ancestries =
    as.integer(
      rowSums(
        ancestry_eGene_status_matrix
      ) == 3
    )
)

# Count the total number of TCF4 genes positive for each outcome
colSums(ancestry_eGene_outcome_matrix)

# Calculate the observed counts
# crossprod() counts positive eGene labels in every
# membership group for every outcome.
# The resulting matrix has:
#   - seven rows representing membership groups
#   - five columns representing ancestry outcomes
observed_ancestry_eGene_counts <- crossprod(
  ancestry_eGene_membership_matrix,
  ancestry_eGene_outcome_matrix
)


# Display the observed counts
observed_ancestry_eGene_counts

# Run 10,000 joint permutations
# Set the random seed for reproducibility
set.seed(12345)


# Define the number of permutations
number_of_permutations <- 10000


# In each permutation:
# 1. Randomly rearrange the gene rows.
# 2. Keep each gene's EUR, AMR and AFR labels together.
# 3. Recalculate the at-least-two and all-three outcomes.
# 4. Count the randomized positive labels in each TCF4 group.
permuted_ancestry_eGene_counts <- replicate(
  number_of_permutations,

  {

    # Randomly rearrange gene row positions
    shuffled_gene_rows <- sample(
      seq_len(
        nrow(
          ancestry_eGene_status_matrix
        )
      ),
      replace = FALSE
    )


    # Move the complete EUR-AMR-AFR pattern together
    shuffled_ancestry_status <-
      ancestry_eGene_status_matrix[
        shuffled_gene_rows,
        ,
        drop = FALSE
      ]


    # Recreate the five outcomes after shuffling
    shuffled_ancestry_outcomes <- cbind(

      EUR_eGene =
        shuffled_ancestry_status[, "EUR"],

      AMR_eGene =
        shuffled_ancestry_status[, "AMR"],

      AFR_eGene =
        shuffled_ancestry_status[, "AFR"],

      At_least_2_Ancestries =
        as.integer(
          rowSums(
            shuffled_ancestry_status
          ) >= 2
        ),

      All_3_Ancestries =
        as.integer(
          rowSums(
            shuffled_ancestry_status
          ) == 3
        )
    )


    # Count randomized positive labels in each membership group
    crossprod(
      ancestry_eGene_membership_matrix,
      shuffled_ancestry_outcomes
    )
  }
)

# Confirm the permutation-array dimensions
dim(permuted_ancestry_eGene_counts)

# Calculate the random expectations
# Calculate the average random count for every comparison

mean_random_ancestry_eGene_counts <- apply(
  permuted_ancestry_eGene_counts,
  MARGIN = c(1, 2),
  FUN = mean
)


# Calculate the standard deviation of the random counts
random_ancestry_eGene_standard_deviation <- apply(
  permuted_ancestry_eGene_counts,
  MARGIN = c(1, 2),
  FUN = stats::sd
)

# Calculate empirical P values
# Repeat the observed-count matrix across all 10,000 permutations
# so each random count can be compared with its corresponding
# observed count
observed_ancestry_eGene_count_array <- array(
  observed_ancestry_eGene_counts,
  dim = dim(permuted_ancestry_eGene_counts))


# Determine whether each randomized count was at least as large
# as the corresponding observed count.
# This is a one-sided enrichment test.
random_at_least_as_large_as_observed <-
  permuted_ancestry_eGene_counts >=
  observed_ancestry_eGene_count_array


# Count the number of permutations at least as large as observed
number_of_random_results_at_least_as_large <- apply(
  random_at_least_as_large_as_observed,
  MARGIN = c(1, 2),
  FUN = sum
)


# Calculate empirical P values using the added-one correction
empirical_ancestry_eGene_p_values <- (
  number_of_random_results_at_least_as_large + 1) / (number_of_permutations + 1)

# Create the complete result table
# Create one row for every combination of:
#   seven TCF4 membership groups
#   five ancestry outcomes
# This produces 35 result rows.
ancestry_eGene_permutation_results <- expand.grid(

  Membership_Group = rownames(
    observed_ancestry_eGene_counts),

  Ancestry_Outcome = colnames(
    observed_ancestry_eGene_counts),

  stringsAsFactors = FALSE)


# Add observed counts, random counts and empirical P values
ancestry_eGene_permutation_results <-
  ancestry_eGene_permutation_results %>%

  dplyr::mutate(

    Total_Genes_in_Group = rep(
      colSums(
        ancestry_eGene_membership_matrix
      ),
      times = ncol(
        observed_ancestry_eGene_counts
      )
    ),

    Observed_eGenes = as.numeric(
      observed_ancestry_eGene_counts
    ),

    Mean_Random_eGenes = as.numeric(
      mean_random_ancestry_eGene_counts
    ),

    Random_Standard_Deviation = as.numeric(
      random_ancestry_eGene_standard_deviation
    ),

    Empirical_P_Value = as.numeric(
      empirical_ancestry_eGene_p_values
    )
  )

# Add effect-size measurements
ancestry_eGene_permutation_results <-
  ancestry_eGene_permutation_results %>%
  dplyr::mutate(

    Observed_Minus_Random = round(
      Observed_eGenes - Mean_Random_eGenes,
      digits = 2),

    Observed_to_Random_Ratio = round(
      Observed_eGenes / Mean_Random_eGenes,
      digits = 3)
  )

# Calculate FDR within each outcome
# This is the more conservative correction.
# It corrects all seven groups and all five outcomes together.
ancestry_eGene_permutation_results <-
  ancestry_eGene_permutation_results %>%

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
  )

# Display the complete results
ancestry_eGene_permutation_results %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

# Display only globally significant enrichments
ancestry_eGene_significant_permutation_results <-
  ancestry_eGene_permutation_results %>%

  dplyr::filter(
    Significant_After_Global_FDR
  ) %>%

  dplyr::arrange(
    Empirical_FDR_Across_All_Tests,
    Empirical_P_Value
  )


ancestry_eGene_significant_permutation_results %>%
  tibble::as_tibble() %>%
print(n = Inf,width = Inf)
