# TCF4 x developing-brain cross-ancestry atlas integration
# Paper: Cross-ancestry atlas of gene, isoform, and splicing regulation in the
# developing human brain (Science 2024; eadh0829)
# https://pubmed.ncbi.nlm.nih.gov/38781368/
# The paper analyzed RNA-seq on GRCh37 with GENCODE v29lift37.
# This paper has 7 supplementary tables
# THIS CODE IS FOR TABLES3
# Compare trimester-specific fetal-brain eGenes with TCF4-associated genes

# Load libraries
library(readxl)  # imports Excel workbooks into R
library(dplyr)
library(tidyr)
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

# Select the Table S3 workbook
# Choose: science.adh0829_table_s3.xlsx
table_s3_file <- file.choose()


# Create unique TCF4 peak-associated gene sets
# Create the NPC peak-associated gene list
npc_all_genes <- npc_peak_anno %>%
  dplyr::transmute(

    # Extract the stable Ensembl gene identifier
    # and remove a version suffix if one exists
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      trimws(
        as.character(
          ENSEMBL
        )
      )
    ),

    # Retain the corresponding gene symbol
    gene_symbol = trimws(
      as.character(
        SYMBOL
      )
    )
  ) %>%

  # Remove rows without valid Ensembl gene identifiers
  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    grepl(
      "^ENSG[0-9]+$",
      ensembl_gene_id
    )
  ) %>%

  # Multiple peaks may be assigned to the same gene.
  # Retain only one row per Ensembl gene.
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

# Check the resulting gene sets
tcf4_gene_set_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_Annotated_Peaks = c(
    nrow(npc_peak_anno),
    nrow(mcclay_peak_anno),
    nrow(forrest_peak_anno)
  ),

  Unique_Peak_Associated_Genes = c(
    nrow(npc_all_genes),
    nrow(mcclay_all_genes),
    nrow(forrest_all_genes)
  )
)


tcf4_gene_set_summary

# Display all sheet names
readxl::excel_sheets(table_s3_file)

# Import first-trimester significant eGenes
s3_Tri1_eGene <- readxl::read_excel(path = table_s3_file,
  sheet = "ST3-1-tri1-eGene")

# Import second-trimester significant eGenes
s3_Tri2_eGene <- readxl::read_excel(path = table_s3_file,
  sheet = "ST3-2-tri2-eGene")

# Inspect the imported data
dim(s3_Tri1_eGene)
dim(s3_Tri2_eGene)

# Show column names
colnames(s3_Tri1_eGene)
colnames(s3_Tri2_eGene)

# Display the first rows
head(s3_Tri1_eGene)
head(s3_Tri2_eGene)

# Clean the Tri1 eGene table
# Create one clearly labeled row per stable Ensembl gene
s3_Tri1_eGene_clean <- s3_Tri1_eGene %>%
  dplyr::transmute(
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      trimws(
        as.character(
          pid
        )
      )
    ),

    trimester = "Tri1",
    primary_variant = sid,
    distance_to_TSS_bp = dist,
    eQTL_slope = slope,
    nominal_p_value = npval,
    permutation_p_value = ppval,
    beta_approximated_p_value = bpval,
    q_value = qval,
    nominal_p_value_threshold = pval_nominal_threshold,
    number_of_variants_tested = nvar) %>%
  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "") %>%

  # Put the strongest result first if a gene is repeated
  dplyr::arrange(q_value,
    nominal_p_value) %>%

  # Retain one row per gene
  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

# Clean the Tri2 eGene table
s3_Tri2_eGene_clean <- s3_Tri2_eGene %>%
  dplyr::transmute(
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      trimws(
        as.character(
          pid
        )
      )
    ),

    trimester = "Tri2",
    primary_variant = sid,
    distance_to_TSS_bp = dist,
    eQTL_slope = slope,
    nominal_p_value = npval,
    permutation_p_value = ppval,
    beta_approximated_p_value = bpval,
    q_value = qval,
    nominal_p_value_threshold = pval_nominal_threshold,
    number_of_variants_tested = nvar
  ) %>%

  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "" ) %>%

  dplyr::arrange(
    q_value,
    nominal_p_value) %>%

  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

# Verify the cleaned gene list
s3_eGene_cleaning_summary <- data.frame(
  Trimester = c(
    "Tri1",
    "Tri2"
  ),
  Original_Rows = c(
    nrow(s3_Tri1_eGene),
    nrow(s3_Tri2_eGene)
  ),
  Unique_eGenes = c(
    nrow(s3_Tri1_eGene_clean),
    nrow(s3_Tri2_eGene_clean)
  ),

  Remaining_Duplicate_Gene_IDs = c(

    sum(duplicated(s3_Tri1_eGene_clean$ensembl_gene_id)),

    sum(duplicated(s3_Tri2_eGene_clean$ensembl_gene_id))
  ),

  Missing_Gene_IDs = c(

    sum(is.na(s3_Tri1_eGene_clean$ensembl_gene_id) |
        s3_Tri1_eGene_clean$ensembl_gene_id == ""
    ),

    sum(is.na(s3_Tri2_eGene_clean$ensembl_gene_id) |
        s3_Tri2_eGene_clean$ensembl_gene_id == ""
    )
  )
)


s3_eGene_cleaning_summary

# Compare the complete NPC, McClay and Forrest sets
# Create a six-row table:
# NPC-Tri1
# NPC-Tri2
# McClay-Tri1
# McClay-Tri2
# Forrest-Tri1
# Forrest-Tri2

s3_TCF4_eGene_overlap_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC",
    "NPC",
    "McClay",
    "McClay",
    "Forrest",
    "Forrest"
  ),

  Trimester = c(
    "Tri1",
    "Tri2",
    "Tri1",
    "Tri2",
    "Tri1",
    "Tri2"
  ),

  Total_TCF4_Genes = c(
    nrow(npc_all_genes),
    nrow(npc_all_genes),
    nrow(mcclay_all_genes),
    nrow(mcclay_all_genes),
    nrow(forrest_all_genes),
    nrow(forrest_all_genes)
  ),

  Significant_eGene_Overlap = c(

    # NPC and Tri1
    sum(npc_all_genes$ensembl_gene_id %in%
        s3_Tri1_eGene_clean$ensembl_gene_id
    ),

    # NPC and Tri2
    sum(npc_all_genes$ensembl_gene_id %in%
        s3_Tri2_eGene_clean$ensembl_gene_id
    ),

    # McClay and Tri1
    sum(mcclay_all_genes$ensembl_gene_id %in%
        s3_Tri1_eGene_clean$ensembl_gene_id
    ),

    # McClay and Tri2
    sum(mcclay_all_genes$ensembl_gene_id %in%
        s3_Tri2_eGene_clean$ensembl_gene_id
    ),

    # Forrest and Tri1
    sum(forrest_all_genes$ensembl_gene_id %in%
        s3_Tri1_eGene_clean$ensembl_gene_id
    ),

    # Forrest and Tri2
    sum(forrest_all_genes$ensembl_gene_id %in%
        s3_Tri2_eGene_clean$ensembl_gene_id
    )
  )
) %>%

  dplyr::mutate(
    Percent_of_TCF4_Set_in_eGene_List = round(
      100 *Significant_eGene_Overlap /Total_TCF4_Genes,digits = 2))


s3_TCF4_eGene_overlap_summary


# Create the five selected TCF4 membership groups
# Create an NPC membership indicator
npc_membership <- npc_all_genes %>%

  dplyr::transmute(
    ensembl_gene_id,
    Present_NPC = TRUE
  )


# Create a McClay membership indicator
mcclay_membership <- mcclay_all_genes %>%

  dplyr::transmute(
    ensembl_gene_id,
    Present_McClay = TRUE
  )


# Create a Forrest membership indicator
forrest_membership <- forrest_all_genes %>%

  dplyr::transmute(
    ensembl_gene_id,
    Present_Forrest = TRUE
  )


# Combine genes from all three datasets
# full_join() retains every gene found in at least one dataset
tcf4_selected_gene_membership <- npc_membership %>%

  dplyr::full_join(
    mcclay_membership,
    by = "ensembl_gene_id"
  ) %>%

  dplyr::full_join(
    forrest_membership,
    by = "ensembl_gene_id"
  ) %>%

  # Genes absent from a particular dataset currently contain NA.
  # Replace those missing values with FALSE.
  dplyr::mutate(

    Present_NPC = tidyr::replace_na(
      Present_NPC,
      FALSE
    ),

    Present_McClay = tidyr::replace_na(
      Present_McClay,
      FALSE
    ),

    Present_Forrest = tidyr::replace_na(
      Present_Forrest,
      FALSE
    )
  )

# Assign the five selected groups
tcf4_selected_gene_membership <-
  tcf4_selected_gene_membership %>%

  dplyr::mutate(

    Selected_Membership_Group = dplyr::case_when(

      # Present only in NPC

      Present_NPC &
        !Present_McClay &
        !Present_Forrest ~
        "NPC only",

      # Present only in McClay

      !Present_NPC &
        Present_McClay &
        !Present_Forrest ~
        "McClay only",

      # Present only in Forrest

      !Present_NPC &
        !Present_McClay &
        Present_Forrest ~
        "Forrest only",

      # Shared by both SH-SY5Y datasets but absent from NPC

      !Present_NPC &
        Present_McClay &
        Present_Forrest ~
        "McClay + Forrest",

      # Present in NPC, McClay and Forrest

      Present_NPC &
        Present_McClay &
        Present_Forrest ~
        "All three",

      # NPC + McClay and NPC + Forrest genes receive NA
      # because they are intentionally excluded

      TRUE ~ NA_character_
    )
  ) %>%

  # Remove the two excluded membership categories

  dplyr::filter(
    !is.na(
      Selected_Membership_Group
    )
  ) %>%

  # Set the order used in tables and plots

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

# Verify the five groups
tcf4_selected_membership_summary <-
  tcf4_selected_gene_membership %>%

  dplyr::group_by(
    Selected_Membership_Group
  ) %>%

  dplyr::summarise(

    Number_of_Genes = dplyr::n(),

    .groups = "drop"
  )


tcf4_selected_membership_summary


# Count genes retained in the selected five-group universe
nrow(tcf4_selected_gene_membership)

# Confirm that 1,103 genes were excluded
10459 -nrow(tcf4_selected_gene_membership)


# Add the Tri1 and Tri2 labels
tcf4_selected_trimester_eGene <-
  tcf4_selected_gene_membership %>%

  dplyr::mutate(

    Detected_Tri1_eGene =
      ensembl_gene_id %in%
      s3_Tri1_eGene_clean$ensembl_gene_id,

    Detected_Tri2_eGene =
      ensembl_gene_id %in%
      s3_Tri2_eGene_clean$ensembl_gene_id,

    Detected_in_Both_Trimesters =
      Detected_Tri1_eGene &
      Detected_Tri2_eGene,

    Detected_Only_in_Tri1_Analysis =
      Detected_Tri1_eGene &
      !Detected_Tri2_eGene,

    Detected_Only_in_Tri2_Analysis =
      !Detected_Tri1_eGene &
      Detected_Tri2_eGene,

    Not_Detected_in_Either_Analysis =
      !Detected_Tri1_eGene &
      !Detected_Tri2_eGene
  )


# Summarize trimester-specific eGenes in the five selected TCF4 groups
s3_selected_membership_eGene_summary <-
  tcf4_selected_trimester_eGene %>%

  dplyr::group_by(
    Selected_Membership_Group
  ) %>%

  dplyr::summarise(

    # Total genes in this TCF4 membership group
    Total_TCF4_Genes = dplyr::n(),

    # Genes detected in the Tri1 eGene analysis
    Tri1_eGenes = sum(
      Detected_Tri1_eGene
    ),

    Percent_Tri1_eGenes = round(100 * Tri1_eGenes /Total_TCF4_Genes, digits = 2),

    # Genes detected in the Tri2 eGene analysis
    Tri2_eGenes = sum( Detected_Tri2_eGene),

    Percent_Tri2_eGenes = round(100 *Tri2_eGenes /Total_TCF4_Genes,digits = 2),

    # Genes detected in both separate trimester analyses
    Detected_in_Both = sum(Detected_in_Both_Trimesters),

    Percent_Detected_in_Both = round(100 *Detected_in_Both /Total_TCF4_Genes,digits = 2),

    # Genes detected in Tri1 but not Tri2
    Detected_Only_in_Tri1 = sum(
      Detected_Only_in_Tri1_Analysis
    ),

    Percent_Detected_Only_in_Tri1 = round(
      100 * Detected_Only_in_Tri1 /Total_TCF4_Genes,digits = 2),

    # Genes detected in Tri2 but not Tri1
    Detected_Only_in_Tri2 = sum(Detected_Only_in_Tri2_Analysis),

    Percent_Detected_Only_in_Tri2 = round(
      100 * Detected_Only_in_Tri2 / Total_TCF4_Genes, digits = 2),

    # Genes not detected in either significant list
    Detected_in_Neither = sum(Not_Detected_in_Either_Analysis),

    Percent_Detected_in_Neither = round(
      100 * Detected_in_Neither /Total_TCF4_Genes,digits = 2),

    # Descriptive difference between the two detection percentages
    Tri1_Minus_Tri2_Percentage = round(
      Percent_Tri1_eGenes -
        Percent_Tri2_eGenes,
      digits = 2
    ),

    .groups = "drop"
  )

# Verify that the four patterns sum to each group total
s3_selected_membership_eGene_summary <-
  s3_selected_membership_eGene_summary %>%
    dplyr::mutate(
    Sum_of_Detection_Patterns =
      Detected_Only_in_Tri1 +
      Detected_Only_in_Tri2 +
      Detected_in_Both +
      Detected_in_Neither,

    Pattern_Counts_Match_Total =
      Sum_of_Detection_Patterns ==
      Total_TCF4_Genes
  )

# Display the complete result
s3_selected_membership_eGene_summary %>%
  tibble::as_tibble() %>%

  print(n = Inf, width = Inf)
