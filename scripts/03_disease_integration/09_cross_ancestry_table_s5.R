# TCF4 x developing-brain cross-ancestry atlas integration
# Paper: Cross-ancestry atlas of gene, isoform, and splicing regulation in the
# developing human brain (Science 2024; eadh0829)
# https://pubmed.ncbi.nlm.nih.gov/38781368/
# The paper analyzed RNA-seq on GRCh37 with GENCODE v29lift37.
# This paper has 7 supplementary tables
# THIS CODE IS FOR TABLES5
# ST5-1a preFOCUS

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

# Select the Table S5 workbook
# Choose: science.adh0829_table_s5.xlsx
table_s5_file <- file.choose()


# Import the schizophrenia isoTWAS pre-FOCUS results
s5_isoTWAS_preFOCUS <- read_excel(
  path = table_s5_file,
  sheet = "ST5-1a-isoTWAS-preFOCUS")

# check the imported dataset
dim(s5_isoTWAS_preFOCUS)
colnames(s5_isoTWAS_preFOCUS)
head(s5_isoTWAS_preFOCUS)

# Clean the isoTWAS transcript results
s5_isoTWAS_preFOCUS_clean <- s5_isoTWAS_preFOCUS |>
  dplyr::transmute(
    #Psychiatric disorder
    Disorder = Trait,

    # Stable Ensembl gene identifier
    # The Gene column currently does not appear to contain version numbers,
    # but this removes them if any are present elsewhere in the dataset
    ensembl_gene_id = sub(
      pattern = "\\..*$",
      replacement = "",
      x = Gene
    ),

    # Gene symbol
    gene_symbol = HGNC,

    # Original transcript identifier, including its version number
    original_transcript_id = Transcript,

    # Stable transcript identifier without the version suffix
    # Example: ENST00000418034.1 becomes ENST00000418034
    ensembl_transcript_id = sub(
      pattern = "\\..*$",
      replacement = "",
      x = Transcript
    ),

    # Genomic and transcript annotation
    chromosome = Chromosome,
    transcript_start = Start,
    transcript_end = End,
    transcript_biotype = Biotype,

    # IsoTWAS association statistics
    isoTWAS_Z_Score = Z,
    isoTWAS_P_Value = P,
    isoTWAS_Permutation_P_Value = Permutation.P,

    # Strongest GWAS variant associated with this transcript result
    Top_GWAS_SNP = Top.GWAS.SNP,

    # P value of the top GWAS variant
    Top_GWAS_SNP_P_Value = Top.GWAS.SNP.1,

    # Study screening and confirmation P values
    Screening_P_Value = Screening.P,
    Confirmation_P_Value = Confirmation.P
  ) |>

  # Remove rows without usable gene or transcript identifiers
  dplyr::filter(
    !is.na(
      ensembl_gene_id
    ),
    ensembl_gene_id != "",
    !is.na(
      ensembl_transcript_id
    ),
    ensembl_transcript_id != ""
  )

# Inspect the cleaned table
dim(s5_isoTWAS_preFOCUS_clean)
head(s5_isoTWAS_preFOCUS_clean)


# Summarize genes, transcripts and significance
# Z > 0: greater genetically predicted transcript abundance is associated with greater schizophrenia risk.
# Z < 0: greater genetically predicted transcript abundance is associated with lower schizophrenia risk.
s5_isoTWAS_preFOCUS_summary <- data.frame(

  Total_isoTWAS_Rows = nrow(s5_isoTWAS_preFOCUS_clean),

  Unique_Ensembl_Genes = dplyr::n_distinct(
    s5_isoTWAS_preFOCUS_clean$ensembl_gene_id
  ),

  Unique_Gene_Symbols = dplyr::n_distinct(
    s5_isoTWAS_preFOCUS_clean$gene_symbol
  ),

  Unique_Transcripts = dplyr::n_distinct(
    s5_isoTWAS_preFOCUS_clean$ensembl_transcript_id
  ),

  Positive_Z_Score_Transcripts = sum(
    s5_isoTWAS_preFOCUS_clean$isoTWAS_Z_Score > 0,
    na.rm = TRUE
  ),

  Negative_Z_Score_Transcripts = sum(
    s5_isoTWAS_preFOCUS_clean$isoTWAS_Z_Score < 0,
    na.rm = TRUE
  ),

  Permutation_P_Below_0_05 = sum(
    s5_isoTWAS_preFOCUS_clean$isoTWAS_Permutation_P_Value < 0.05,
    na.rm = TRUE
  ),

  Smallest_Permutation_P = min(
    s5_isoTWAS_preFOCUS_clean$isoTWAS_Permutation_P_Value,
    na.rm = TRUE
  ),

  Largest_Permutation_P = max(
    s5_isoTWAS_preFOCUS_clean$isoTWAS_Permutation_P_Value,
    na.rm = TRUE
  )
)

s5_isoTWAS_preFOCUS_summary

# Check whether individual transcripts appear more than once:
# Count the number of rows for each stable transcript identifier
s5_isoTWAS_transcript_counts <- s5_isoTWAS_preFOCUS_clean |>
  dplyr::count(
    ensembl_transcript_id,
    name = "Number_of_Rows"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      Number_of_Rows
    )
  )

# Display only transcripts represented by multiple rows
s5_isoTWAS_repeated_transcripts <- s5_isoTWAS_transcript_counts |>
  dplyr::filter(
    Number_of_Rows > 1
  )

s5_isoTWAS_repeated_transcripts


# Check how many associated transcripts each gene has
# Count distinct associated transcripts for every gene
s5_isoTWAS_transcripts_per_gene <- s5_isoTWAS_preFOCUS_clean |>
  dplyr::group_by(
    ensembl_gene_id,
    gene_symbol
  ) |>
  dplyr::summarise(

    Number_of_Associated_Transcripts = dplyr::n_distinct(
      ensembl_transcript_id
    ),

    .groups = "drop"
  ) |>

  dplyr::arrange(
    dplyr::desc(
      Number_of_Associated_Transcripts
    )
  )

head(s5_isoTWAS_transcripts_per_gene,20)

# Why are nine transcript IDs repeated?
# we found nine stable transcript IDs represented twice
# Inspect all rows belonging to repeated stable transcript identifiers
s5_isoTWAS_repeated_transcript_details <- s5_isoTWAS_preFOCUS_clean |>
  dplyr::filter(
    ensembl_transcript_id %in%
      s5_isoTWAS_repeated_transcripts$ensembl_transcript_id
  ) |>

  dplyr::arrange(
    ensembl_transcript_id,
    isoTWAS_Permutation_P_Value
  ) |>

  dplyr::select(
    ensembl_gene_id,
    gene_symbol,
    original_transcript_id,
    ensembl_transcript_id,
    isoTWAS_Z_Score,
    isoTWAS_P_Value,
    isoTWAS_Permutation_P_Value,
    Top_GWAS_SNP,
    Top_GWAS_SNP_P_Value,
    Screening_P_Value,
    Confirmation_P_Value
  )

s5_isoTWAS_repeated_transcript_details

# we have 271 Ensembl genes but only 251 gene symbols. First, check for missing symbols:
# Identify Ensembl genes without a gene symbol
s5_isoTWAS_missing_gene_symbols <- s5_isoTWAS_preFOCUS_clean |>
  dplyr::filter(
    is.na(
      gene_symbol
    ) |
      gene_symbol == ""
  ) |>

  dplyr::distinct(
    ensembl_gene_id,
    gene_symbol
  )

s5_isoTWAS_missing_gene_symbols

# Then check whether one gene symbol maps to multiple Ensembl genes:
# Identify gene symbols connected to multiple Ensembl gene identifiers
s5_isoTWAS_symbols_with_multiple_gene_IDs <- s5_isoTWAS_preFOCUS_clean |>
  dplyr::filter(
    !is.na(
      gene_symbol
    ),
    gene_symbol != ""
  ) |>

  dplyr::group_by(
    gene_symbol
  ) |>

  dplyr::summarise(

    Number_of_Ensembl_Genes = dplyr::n_distinct(
      ensembl_gene_id
    ),

    Ensembl_Gene_IDs = paste(
      sort(
        unique(
          ensembl_gene_id
        )
      ),
      collapse = "; "
    ),

    .groups = "drop"
  ) |>

  dplyr::filter(
    Number_of_Ensembl_Genes > 1
  ) |>

  dplyr::arrange(
    dplyr::desc(
      Number_of_Ensembl_Genes
    )
  )

s5_isoTWAS_symbols_with_multiple_gene_IDs

# The repeated transcript rows are not caused by transcript-version removal.
# Each pair has the exact same original transcript ID, gene, Z score,
# ordinary P value, top GWAS SNP and top-SNP P value. only permutation P-value differs
# I found 21 Ensembl genes without gene symbols. So, We should use ensembl_gene_id for
# all TCF4 overlaps and treat gene_symbol only as a readable label.
# Count missing and empty gene symbols separately
s5_isoTWAS_gene_symbol_check <- data.frame(

  Unique_Ensembl_Genes = dplyr::n_distinct(
    s5_isoTWAS_preFOCUS_clean$ensembl_gene_id
  ),

  Missing_Gene_Symbols = sum(
    is.na(
      s5_isoTWAS_preFOCUS_clean$gene_symbol
    )
  ),

  Empty_Gene_Symbols = sum(
    !is.na(
      s5_isoTWAS_preFOCUS_clean$gene_symbol
    ) &
      s5_isoTWAS_preFOCUS_clean$gene_symbol == ""
  ),

  Unique_Nonmissing_Gene_Symbols = dplyr::n_distinct(
    s5_isoTWAS_preFOCUS_clean$gene_symbol[
      !is.na(
        s5_isoTWAS_preFOCUS_clean$gene_symbol
      ) &
        s5_isoTWAS_preFOCUS_clean$gene_symbol != ""
    ]
  )
)

s5_isoTWAS_gene_symbol_check

# Check whether any Ensembl gene is connected to multiple symbols
s5_isoTWAS_genes_with_multiple_symbols <- s5_isoTWAS_preFOCUS_clean |>
  dplyr::filter(
    !is.na(gene_symbol),gene_symbol != "") |>

  dplyr::group_by(ensembl_gene_id) |>

  dplyr::summarise(
    Number_of_Gene_Symbols = dplyr::n_distinct(
      gene_symbol
    ),

    Gene_Symbols = paste(
      sort(
        unique(
          gene_symbol
        )
      ),
      collapse = "; "
    ),

    .groups = "drop"
  ) |>

  dplyr::filter(Number_of_Gene_Symbols > 1)


s5_isoTWAS_genes_with_multiple_symbols


# Create one row per Ensembl gene for the TCF4 overlap analysis
s5_isoTWAS_preFOCUS_clean  # what does our clean dataset look like?


s5_isoTWAS_gene_summary <- s5_isoTWAS_preFOCUS_clean |>
  dplyr::group_by(
    ensembl_gene_id
  ) |>
  dplyr::summarise(

    # Retain the available gene symbol
    gene_symbol = dplyr::first(
      gene_symbol
    ),

    # Count distinct associated transcripts
    Number_of_Associated_Transcripts = dplyr::n_distinct(
      ensembl_transcript_id
    ),

    # Count distinct transcripts with positive Z scores
    Number_of_Positive_Z_Transcripts = dplyr::n_distinct(
      ensembl_transcript_id[
        isoTWAS_Z_Score > 0
      ]
    ),

    # Count distinct transcripts with negative Z scores
    Number_of_Negative_Z_Transcripts = dplyr::n_distinct(
      ensembl_transcript_id[
        isoTWAS_Z_Score < 0
      ]
    ),

    # Combine stable transcript identifiers
    Associated_Transcript_IDs = paste(
      sort(
        unique(
          ensembl_transcript_id
        )
      ),
      collapse = "; "
    ),

    # Retain the smallest permutation P value for descriptive ranking
    # (optional:if we decided to report it)
    Minimum_Permutation_P_Value = min(
      isoTWAS_Permutation_P_Value,
      na.rm = TRUE
    ),
    # Identify the transcript with the smallest permutation P value
    Most_Significant_Transcript = ensembl_transcript_id[
      which.min(
        isoTWAS_Permutation_P_Value
      )
    ],

    # Retain the corresponding Z score
    Z_Score_of_Most_Significant_Transcript = isoTWAS_Z_Score[
      which.min(
        isoTWAS_Permutation_P_Value
      )
    ],

    .groups = "drop"
  ) |>

  dplyr::mutate(
    Direction_Pattern = dplyr::case_when(

      Number_of_Positive_Z_Transcripts > 0 &
        Number_of_Negative_Z_Transcripts > 0 ~
        "Both positive and negative",

      Number_of_Positive_Z_Transcripts > 0 ~
        "Positive only",

      Number_of_Negative_Z_Transcripts > 0 ~
        "Negative only",

      TRUE ~"No direction"
    ),

    Significant_SCZ_isoTWAS_Gene = TRUE) |>

  dplyr::arrange(Minimum_Permutation_P_Value)

# check the dimesion our TWAS gene dataset
dim(s5_isoTWAS_gene_summary)

# Confirm that every Ensembl gene appears once
sum(duplicated(s5_isoTWAS_gene_summary$ensembl_gene_id))

# Inspect the first 10 genes
head(s5_isoTWAS_gene_summary,10)

# Assign all TCF4-associated genes to mutually exclusive membership groups
# Create mutually exclusive TCF4 membership groups

# Add an NPC membership indicator
npc_gene_membership <- npc_all_genes |>
  dplyr::transmute(
    ensembl_gene_id,
    NPC_gene_symbol = gene_symbol,
    In_NPC = TRUE
  )

# Add a McClay membership indicator
mcclay_gene_membership <- mcclay_all_genes |>
  dplyr::transmute(
    ensembl_gene_id,
    McClay_gene_symbol = gene_symbol,
    In_McClay = TRUE
  )

# Add a Forrest membership indicator
forrest_gene_membership <- forrest_all_genes |>
  dplyr::transmute(
    ensembl_gene_id,
    Forrest_gene_symbol = gene_symbol,
    In_Forrest = TRUE
  )

# Combine all three gene sets
tcf4_gene_membership <- npc_gene_membership |>
  dplyr::full_join(
    mcclay_gene_membership,
    by = "ensembl_gene_id"
  ) |>

  dplyr::full_join(
    forrest_gene_membership,
    by = "ensembl_gene_id"
  ) |>

  dplyr::mutate(

    # Recover the available gene symbol
    gene_symbol = dplyr::coalesce(
      NPC_gene_symbol,
      McClay_gene_symbol,
      Forrest_gene_symbol
    ),

    # Genes absent from a dataset receive FALSE
    In_NPC = tidyr::replace_na(
      In_NPC,
      FALSE
    ),

    In_McClay = tidyr::replace_na(
      In_McClay,
      FALSE
    ),

    In_Forrest = tidyr::replace_na(
      In_Forrest,
      FALSE
    ),

    # Assign every gene to one mutually exclusive category
    Membership_Group = dplyr::case_when(

      In_NPC &
        !In_McClay &
        !In_Forrest ~
        "NPC only",

      !In_NPC &
        In_McClay &
        !In_Forrest ~
        "McClay only",

      !In_NPC &
        !In_McClay &
        In_Forrest ~
        "Forrest only",

      In_NPC &
        In_McClay &
        !In_Forrest ~
        "NPC + McClay",

      In_NPC &
        !In_McClay &
        In_Forrest ~
        "NPC + Forrest",

      !In_NPC &
        In_McClay &
        In_Forrest ~
        "McClay + Forrest",

      In_NPC &
        In_McClay &
        In_Forrest ~
        "All three"
    )
  ) |>

  dplyr::select(
    ensembl_gene_id,
    gene_symbol,
    In_NPC,
    In_McClay,
    In_Forrest,
    Membership_Group
  )

# Confirm the seven complete membership groups
tcf4_membership_summary <- tcf4_gene_membership |>
  dplyr::count(
    Membership_Group,
    name = "Number_of_Genes"
  )

tcf4_membership_summary


# Compare complete TCF4 datasets with SCZ isoTWAS genes
# Extract the 271 unique SCZ isoTWAS Ensembl gene identifiers
SCZ_isoTWAS_gene_IDs <- unique(
  s5_isoTWAS_gene_summary$ensembl_gene_id)

# Calculate the number and percentage of overlapping genes
complete_TCF4_isoTWAS_overlap_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_TCF4_Peak_Associated_Genes = c(
    nrow(
      npc_all_genes
    ),
    nrow(
      mcclay_all_genes
    ),
    nrow(
      forrest_all_genes
    )
  ),

  SCZ_isoTWAS_Overlap_Genes = c(

    sum(
      npc_all_genes$ensembl_gene_id %in%
        SCZ_isoTWAS_gene_IDs
    ),

    sum(
      mcclay_all_genes$ensembl_gene_id %in%
        SCZ_isoTWAS_gene_IDs
    ),

    sum(
      forrest_all_genes$ensembl_gene_id %in%
        SCZ_isoTWAS_gene_IDs
    )
  )
)

# Calculate the percentage of each TCF4 set found in the SCZ isoTWAS list
complete_TCF4_isoTWAS_overlap_summary <-
  complete_TCF4_isoTWAS_overlap_summary |>

  dplyr::mutate(
    Percent_of_TCF4_Set_in_SCZ_isoTWAS = round(
      100 *
        SCZ_isoTWAS_Overlap_Genes /
        Total_TCF4_Peak_Associated_Genes,
      digits = 2
    )
  )


complete_TCF4_isoTWAS_overlap_summary

# Find the exact NPC and SCZ isoTWAS overlap
NPC_SCZ_isoTWAS_overlap <- npc_all_genes |>

  dplyr::rename(
    TCF4_gene_symbol = gene_symbol
  ) |>

  dplyr::inner_join(
    s5_isoTWAS_gene_summary |>
      dplyr::rename(
        isoTWAS_gene_symbol = gene_symbol
      ),
    by = "ensembl_gene_id"
  ) |>

  dplyr::mutate(
    TCF4_Dataset = "NPC",
    .before = 1
  ) |>

  dplyr::arrange(
    Minimum_Permutation_P_Value
  )


NPC_SCZ_isoTWAS_overlap


# Find the exact McClay and SCZ isoTWAS overlap
McClay_SCZ_isoTWAS_overlap <- mcclay_all_genes |>

  dplyr::rename(
    TCF4_gene_symbol = gene_symbol
  ) |>

  dplyr::inner_join(
    s5_isoTWAS_gene_summary |>
      dplyr::rename(
        isoTWAS_gene_symbol = gene_symbol
      ),
    by = "ensembl_gene_id"
  ) |>

  dplyr::mutate(
    TCF4_Dataset = "McClay",
    .before = 1
  ) |>

  dplyr::arrange(
    Minimum_Permutation_P_Value
  )

McClay_SCZ_isoTWAS_overlap


# Find the exact Forrest and SCZ isoTWAS overlap
Forrest_SCZ_isoTWAS_overlap <- forrest_all_genes |>

  dplyr::rename(
    TCF4_gene_symbol = gene_symbol
  ) |>

  dplyr::inner_join(
    s5_isoTWAS_gene_summary |>
      dplyr::rename(
        isoTWAS_gene_symbol = gene_symbol
      ),
    by = "ensembl_gene_id"
  ) |>

  dplyr::mutate(
    TCF4_Dataset = "Forrest",
    .before = 1
  ) |>

  dplyr::arrange(
    Minimum_Permutation_P_Value
  )


Forrest_SCZ_isoTWAS_overlap

# Combine all complete-dataset overlaps
# A gene can appear more than once if it belongs to multiple TCF4 datasets
complete_TCF4_isoTWAS_overlap_genes <- dplyr::bind_rows(
  NPC_SCZ_isoTWAS_overlap,
  McClay_SCZ_isoTWAS_overlap,
  Forrest_SCZ_isoTWAS_overlap
)


complete_TCF4_isoTWAS_overlap_genes


# Confirm overlap counts using the exact joined tables
complete_TCF4_isoTWAS_join_check <- data.frame(

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Overlap_From_Joined_Table = c(
    nrow(
      NPC_SCZ_isoTWAS_overlap
    ),
    nrow(
      McClay_SCZ_isoTWAS_overlap
    ),
    nrow(
      Forrest_SCZ_isoTWAS_overlap
    )
  )
)

complete_TCF4_isoTWAS_join_check

# Count positive-only, negative-only and mixed-direction genes
complete_TCF4_isoTWAS_direction_summary <-
  complete_TCF4_isoTWAS_overlap_genes |>

  dplyr::count(
    TCF4_Dataset,
    Direction_Pattern,
    name = "Number_of_SCZ_isoTWAS_Genes"
  ) |>

  dplyr::arrange(
    TCF4_Dataset,
    Direction_Pattern
  )

complete_TCF4_isoTWAS_direction_summary


# Exploratory Fisher tests within the TCF4 gene universe
# Create the union of genes found in at least one TCF4 dataset
TCF4_union_gene_IDs <- unique(
  c(
    npc_all_genes$ensembl_gene_id,
    mcclay_all_genes$ensembl_gene_id,
    forrest_all_genes$ensembl_gene_id
  )
)

# Confirm the size of the TCF4 union
length(TCF4_union_gene_IDs)

# Take a quick loot to the genes
NPC_SCZ_isoTWAS_overlap$TCF4_gene_symbol[NPC_SCZ_isoTWAS_overlap$Direction_Pattern == "Positive only"]
NPC_SCZ_isoTWAS_overlap$TCF4_gene_symbol[NPC_SCZ_isoTWAS_overlap$Direction_Pattern == "Negative only"]
NPC_SCZ_isoTWAS_overlap$isoTWAS_gene_symbol
McClay_SCZ_isoTWAS_overlap$isoTWAS_gene_symbol
Forrest_SCZ_isoTWAS_overlap$isoTWAS_gene_symbol




################################################################################
# ST5-1b FOCUS fine-mapping

# Import the fine-mapped schizophrenia isoTWAS worksheet
s5_isoTWAS_finemapped <- read_excel(
  path = table_s5_file,
  sheet = "ST5-1b-isoTWAS-finemapped"
)

# Check the number of rows and columns
dim(s5_isoTWAS_finemapped)

# Display the exact column names
colnames(s5_isoTWAS_finemapped)

# Display the first six rows
head(s5_isoTWAS_finemapped)

# Inspect the type of every column
str(s5_isoTWAS_finemapped)

# Clean the fine-mapped data
# Clean the ST5-1b fine-mapped isoTWAS results
# Create a clean table with clearly labelled columns
s5_isoTWAS_finemapped_clean <- s5_isoTWAS_finemapped |>

  dplyr::transmute(

    # Gene information
    gene_symbol = HGNC,

    # Remove a possible Ensembl gene version number
    # For example:
    # ENSG00000123456.5 becomes ENSG00000123456
    ensembl_gene_id = sub(
      pattern = "\\..*$",
      replacement = "",
      x = Ensembl
    ),

    # Keep the original transcript identifier with its version number
    original_transcript_id = Transcript,

    # Create a stable transcript identifier without its version number
    # For example:
    # ENST00000418034.1 becomes ENST00000418034
    ensembl_transcript_id = sub(
      pattern = "\\..*$",
      replacement = "",
      x = Transcript
    ),

    # Study trait
    trait = Trait,

    # Genomic information
    chromosome = Chromosome,
    gene_start = Start,
    gene_end = End,
    gene_biotype = Biotype,

    # Original isoTWAS statistics
    isoTWAS_Z_score = Z,
    isoTWAS_P_value = P,
    permutation_P_value = `Permutation Pvalue`,

    # Most strongly associated GWAS variant
    top_GWAS_SNP = `Top GWAS SNP`,
    top_GWAS_SNP_P_value = `TOP GWAS SNP Pvalue`,

    # Transcript screening and confirmation statistics
    screening_P_value = `Screening Pvalue`,
    confirmation_P_value = `Confirmation Pvalue`,

    # FOCUS fine-mapping results
    posterior_inclusion_probability = PIP,
    included_in_90_percent_credible_set = `In credible set`,

    # Loss-of-function intolerance annotation
    pLI_score = pLI
  ) |>

  # Create a readable direction label from the isoTWAS Z score
  dplyr::mutate(
    isoTWAS_direction = dplyr::case_when(
      isoTWAS_Z_score > 0 ~ "Positive",
      isoTWAS_Z_score < 0 ~ "Negative",
      isoTWAS_Z_score == 0 ~ "No direction",
      TRUE ~ NA_character_
    )
  )

# Check the cleaned table
# Check its dimensions
dim(s5_isoTWAS_finemapped_clean)

# Display the first six rows
head(s5_isoTWAS_finemapped_clean)

# Confirm the type of every cleaned column
str(s5_isoTWAS_finemapped_clean)


# Summarize the cleaned fine-mapped isoTWAS dataset
# Count the main features in the dataset
s5_isoTWAS_finemapped_summary <- data.frame(

  Total_Fine_Mapped_Rows = nrow(
    s5_isoTWAS_finemapped_clean
  ),

  Unique_Stable_Transcripts = dplyr::n_distinct(
    s5_isoTWAS_finemapped_clean$ensembl_transcript_id,
    na.rm = TRUE
  ),

  Unique_Ensembl_Genes = dplyr::n_distinct(
    s5_isoTWAS_finemapped_clean$ensembl_gene_id,
    na.rm = TRUE
  ),

  Unique_Nonmissing_Gene_Symbols = dplyr::n_distinct(
    s5_isoTWAS_finemapped_clean$gene_symbol,
    na.rm = TRUE
  ),

  Missing_or_Empty_Gene_Symbols = sum(
    is.na(
      s5_isoTWAS_finemapped_clean$gene_symbol
    ) |
      s5_isoTWAS_finemapped_clean$gene_symbol == ""
  ),

  Positive_Z_Transcripts = sum(
    s5_isoTWAS_finemapped_clean$isoTWAS_direction == "Positive",
    na.rm = TRUE
  ),

  Negative_Z_Transcripts = sum(
    s5_isoTWAS_finemapped_clean$isoTWAS_direction == "Negative",
    na.rm = TRUE
  ),

  Transcripts_in_90_Percent_Credible_Set = sum(
    s5_isoTWAS_finemapped_clean$
      included_in_90_percent_credible_set,
    na.rm = TRUE
  )
)

s5_isoTWAS_finemapped_summary


# Confirm that every row concerns schizophrenia
table(s5_isoTWAS_finemapped_clean$trait,useNA = "ifany")

# Count positive and negative isoTWAS associations
table(s5_isoTWAS_finemapped_clean$isoTWAS_direction,useNA = "ifany")

# Check whether every transcript is included in a credible set
table(s5_isoTWAS_finemapped_clean$included_in_90_percent_credible_set,useNA = "ifany")

# Examine the different gene biotypes
table(s5_isoTWAS_finemapped_clean$gene_biotype,useNA = "ifany")

# Examine the PIP values
# Display the minimum, quartiles, median, mean and maximum PIP
summary(
  s5_isoTWAS_finemapped_clean$
    posterior_inclusion_probability
)

# Count transcripts with PIP values close to 1
sum(
  s5_isoTWAS_finemapped_clean$
    posterior_inclusion_probability >= 0.90,
  na.rm = TRUE
)

# Count transcripts with PIP values below 0.01
sum(
  s5_isoTWAS_finemapped_clean$
    posterior_inclusion_probability < 0.01,
  na.rm = TRUE
)

# Create readable PIP categories
# Place transcripts into descriptive PIP categories
s5_isoTWAS_finemapped_clean <- s5_isoTWAS_finemapped_clean |>

  dplyr::mutate(
    PIP_category = dplyr::case_when(

      posterior_inclusion_probability >= 0.90 ~
        "Very high: 0.90–1.00",

      posterior_inclusion_probability >= 0.50 ~
        "High: 0.50–0.89",

      posterior_inclusion_probability >= 0.10 ~
        "Moderate: 0.10–0.49",

      posterior_inclusion_probability >= 0.01 ~
        "Low: 0.01–0.09",

      posterior_inclusion_probability < 0.01 ~
        "Very low: below 0.01",

      TRUE ~
        NA_character_
    )
  )

# Count transcripts in each category
s5_isoTWAS_PIP_category_summary <-
  s5_isoTWAS_finemapped_clean |>

  dplyr::count(
    PIP_category,
    name = "Number_of_Transcripts"
  ) |>

  dplyr::arrange(
    dplyr::desc(
      Number_of_Transcripts
    )
  )

s5_isoTWAS_PIP_category_summary

# Check for repeated transcript IDs
# Determine whether any stable transcript ID appears more than once
s5_isoTWAS_finemapped_transcript_counts <-
  s5_isoTWAS_finemapped_clean |>

  dplyr::count(
    ensembl_transcript_id,
    name = "Number_of_Rows"
  ) |>

  dplyr::arrange(
    dplyr::desc(
      Number_of_Rows
    )
  )

# Keep only transcripts represented by multiple rows
s5_isoTWAS_finemapped_repeated_transcripts <-
  s5_isoTWAS_finemapped_transcript_counts |>

  dplyr::filter(
    Number_of_Rows > 1
  )

s5_isoTWAS_finemapped_repeated_transcripts


# Count fine-mapped transcripts per gene
# Determine whether individual genes contain multiple fine-mapped transcripts
s5_isoTWAS_finemapped_gene_counts <-
  s5_isoTWAS_finemapped_clean |>

  dplyr::group_by(
    ensembl_gene_id,
    gene_symbol
  ) |>

  dplyr::summarise(

    Number_of_Fine_Mapped_Transcripts = dplyr::n_distinct(
      ensembl_transcript_id
    ),

    Highest_PIP = max(
      posterior_inclusion_probability,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) |>

  dplyr::arrange(
    dplyr::desc(
      Number_of_Fine_Mapped_Transcripts
    ),
    dplyr::desc(
      Highest_PIP
    )
  )

s5_isoTWAS_finemapped_gene_counts

# Count how many genes have one, two or three fine-mapped transcripts
s5_isoTWAS_transcripts_per_gene_distribution <-
  s5_isoTWAS_finemapped_gene_counts |>

  dplyr::count(
    Number_of_Fine_Mapped_Transcripts,
    name = "Number_of_Genes"
  ) |>

  dplyr::arrange(
    Number_of_Fine_Mapped_Transcripts
  )

s5_isoTWAS_transcripts_per_gene_distribution

# overlap with the three TCF4 datasets
# Identify fine-mapped SCZ transcripts in TCF4-associated genes
# Identify fine-mapped transcripts belonging to NPC TCF4-associated genes
NPC_SCZ_isoTWAS_finemapped_overlap <-
  npc_all_genes |>

  dplyr::select(
    ensembl_gene_id
  ) |>

  dplyr::distinct() |>

  dplyr::inner_join(
    s5_isoTWAS_finemapped_clean,
    by = "ensembl_gene_id"
  ) |>

  dplyr::arrange(
    dplyr::desc(
      posterior_inclusion_probability
    ),
    isoTWAS_P_value
  )


# Identify fine-mapped transcripts belonging to McClay TCF4-associated genes
McClay_SCZ_isoTWAS_finemapped_overlap <-
  mcclay_all_genes |>

  dplyr::select(
    ensembl_gene_id
  ) |>

  dplyr::distinct() |>

  dplyr::inner_join(
    s5_isoTWAS_finemapped_clean,
    by = "ensembl_gene_id"
  ) |>

  dplyr::arrange(
    dplyr::desc(
      posterior_inclusion_probability
    ),
    isoTWAS_P_value
  )


# Identify fine-mapped transcripts belonging to Forrest TCF4-associated genes
Forrest_SCZ_isoTWAS_finemapped_overlap <-
  forrest_all_genes |>

  dplyr::select(
    ensembl_gene_id
  ) |>

  dplyr::distinct() |>

  dplyr::inner_join(
    s5_isoTWAS_finemapped_clean,
    by = "ensembl_gene_id"
  ) |>

  dplyr::arrange(
    dplyr::desc(
      posterior_inclusion_probability
    ),
    isoTWAS_P_value
  )

# Create the first TCF4 comparison
# Summarize fine-mapped transcript and gene overlaps
s5_isoTWAS_finemapped_TCF4_overlap_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  Total_TCF4_Genes = c(
    nrow(npc_all_genes),
    nrow(mcclay_all_genes),
    nrow(forrest_all_genes)
  ),

  Fine_Mapped_Transcript_Overlap = c(
    dplyr::n_distinct(
      NPC_SCZ_isoTWAS_finemapped_overlap$
        ensembl_transcript_id
    ),

    dplyr::n_distinct(
      McClay_SCZ_isoTWAS_finemapped_overlap$
        ensembl_transcript_id
    ),

    dplyr::n_distinct(
      Forrest_SCZ_isoTWAS_finemapped_overlap$
        ensembl_transcript_id
    )
  ),

  Fine_Mapped_Gene_Overlap = c(
    dplyr::n_distinct(
      NPC_SCZ_isoTWAS_finemapped_overlap$
        ensembl_gene_id
    ),

    dplyr::n_distinct(
      McClay_SCZ_isoTWAS_finemapped_overlap$
        ensembl_gene_id
    ),

    dplyr::n_distinct(
      Forrest_SCZ_isoTWAS_finemapped_overlap$
        ensembl_gene_id
    )
  ),

  Very_High_PIP_Transcripts = c(
    sum(
      NPC_SCZ_isoTWAS_finemapped_overlap$
        posterior_inclusion_probability >= 0.90
    ),

    sum(
      McClay_SCZ_isoTWAS_finemapped_overlap$
        posterior_inclusion_probability >= 0.90
    ),

    sum(
      Forrest_SCZ_isoTWAS_finemapped_overlap$
        posterior_inclusion_probability >= 0.90
    )
  )
)

s5_isoTWAS_finemapped_TCF4_overlap_summary


# Compare preFOCUS and fine-mapped isoTWAS overlaps
# Create the unique preFOCUS gene list
s5_preFOCUS_gene_IDs <- unique(
  s5_isoTWAS_preFOCUS_clean$ensembl_gene_id
)

# Create the unique fine-mapped gene list
s5_finemapped_gene_IDs <- unique(
  s5_isoTWAS_finemapped_clean$ensembl_gene_id
)


# Find NPC TCF4 genes represented in the preFOCUS results
NPC_preFOCUS_gene_IDs <- intersect(
  npc_all_genes$ensembl_gene_id,
  s5_preFOCUS_gene_IDs
)

# Determine which NPC preFOCUS genes survive fine-mapping
NPC_retained_finemapped_gene_IDs <- intersect(
  NPC_preFOCUS_gene_IDs,
  s5_finemapped_gene_IDs
)


# Find McClay TCF4 genes represented in the preFOCUS results
McClay_preFOCUS_gene_IDs <- intersect(
  mcclay_all_genes$ensembl_gene_id,
  s5_preFOCUS_gene_IDs
)

# Determine which McClay preFOCUS genes survive fine-mapping
McClay_retained_finemapped_gene_IDs <- intersect(
  McClay_preFOCUS_gene_IDs,
  s5_finemapped_gene_IDs
)


# Find Forrest TCF4 genes represented in the preFOCUS results
Forrest_preFOCUS_gene_IDs <- intersect(
  forrest_all_genes$ensembl_gene_id,
  s5_preFOCUS_gene_IDs
)

# Determine which Forrest preFOCUS genes survive fine-mapping
Forrest_retained_finemapped_gene_IDs <- intersect(
  Forrest_preFOCUS_gene_IDs,
  s5_finemapped_gene_IDs
)

# Create a gene-level retention table
# Summarize how many preFOCUS genes remain after fine-mapping
s5_isoTWAS_gene_finemapping_retention_summary <- data.frame(

  TCF4_Dataset = c(
    "NPC",
    "McClay",
    "Forrest"
  ),

  PreFOCUS_Overlap_Genes = c(
    length(
      NPC_preFOCUS_gene_IDs
    ),

    length(
      McClay_preFOCUS_gene_IDs
    ),

    length(
      Forrest_preFOCUS_gene_IDs
    )
  ),

  Retained_Fine_Mapped_Genes = c(
    length(
      NPC_retained_finemapped_gene_IDs
    ),

    length(
      McClay_retained_finemapped_gene_IDs
    ),

    length(
      Forrest_retained_finemapped_gene_IDs
    )
  )
) |>

  dplyr::mutate(
    Percent_PreFOCUS_Genes_Retained = round(
      100 *
        Retained_Fine_Mapped_Genes /
        PreFOCUS_Overlap_Genes,
      digits = 2
    )
  )

s5_isoTWAS_gene_finemapping_retention_summary

#########################
# TEST
# Create independent TCF4 contexts
#########################
# Create independent NPC and SH-SY5Y TCF4 groups

tcf4_gene_context <- tcf4_gene_membership |>

  dplyr::transmute(
    ensembl_gene_id,
    Membership_Group,

    TCF4_Cell_Context = dplyr::case_when(

      # TCF4-associated only in NPC
      Membership_Group == "NPC only" ~
        "NPC only",

      # TCF4-associated in one or both SH-SY5Y datasets,
      # but absent from NPC
      Membership_Group %in% c(
        "McClay only",
        "Forrest only",
        "McClay + Forrest"
      ) ~
        "SH-SY5Y only",

      # TCF4-associated with NPC and at least one SH-SY5Y dataset
      Membership_Group %in% c(
        "NPC + McClay",
        "NPC + Forrest",
        "All three"
      ) ~
        "Shared NPC and SH-SY5Y",

      TRUE ~
        NA_character_
    )
  )

# Start with the preFOCUS genes
# Create one row for each unique preFOCUS SCZ isoTWAS gene
s5_preFOCUS_unique_genes <-
  s5_isoTWAS_preFOCUS_clean |>

  dplyr::distinct(
    ensembl_gene_id
  )

# Add TCF4 cell-context information
# inner_join keeps only preFOCUS genes associated with at least one TCF4 dataset
s5_preFOCUS_TCF4_context <-
  s5_preFOCUS_unique_genes |>

  dplyr::inner_join(
    tcf4_gene_context,
    by = "ensembl_gene_id"
  ) |>

  # Mark whether each initial gene survived FOCUS fine-mapping
  dplyr::mutate(
    Retained_After_FOCUS =
      ensembl_gene_id %in%
      s5_isoTWAS_finemapped_clean$ensembl_gene_id
  )

# Calculate retention in the independent groups
# Calculate fine-mapping retention for each TCF4 context
s5_TCF4_context_retention_summary <-
  s5_preFOCUS_TCF4_context |>

  dplyr::group_by(
    TCF4_Cell_Context
  ) |>

  dplyr::summarise(

    PreFOCUS_Genes = dplyr::n(),

    Retained_Fine_Mapped_Genes = sum(
      Retained_After_FOCUS
    ),

    Not_Retained_Genes = sum(
      !Retained_After_FOCUS
    ),

    Percent_Retained = round(
      100 *
        Retained_Fine_Mapped_Genes /
        PreFOCUS_Genes,
      digits = 2
    ),

    .groups = "drop"
  )

s5_TCF4_context_retention_summary


# Exploratory Fisher comparison
# Keep the two cell-specific groups
s5_cell_specific_retention_data <-
  s5_preFOCUS_TCF4_context |>

  dplyr::filter(
    TCF4_Cell_Context %in% c(
      "NPC only",
      "SH-SY5Y only"
    )
  )

# Create the two-by-two table
s5_cell_specific_retention_table <- table(
  s5_cell_specific_retention_data$TCF4_Cell_Context,
  s5_cell_specific_retention_data$Retained_After_FOCUS
)

s5_cell_specific_retention_table

# Compare the retention proportions
s5_cell_specific_retention_fisher <-
  fisher.test(
    s5_cell_specific_retention_table,
    alternative = "two.sided"
  )

s5_cell_specific_retention_fisher

# Compare PIP strength
# Connect fine-mapped transcripts to their TCF4 cell context
s5_finemapped_TCF4_context <-
  s5_isoTWAS_finemapped_clean |>

  dplyr::inner_join(
    tcf4_gene_context,
    by = "ensembl_gene_id"
  )

# Collapse to one row per gene
# The highest transcript PIP is retained for each gene
s5_finemapped_TCF4_gene_PIP_summary <-
  s5_finemapped_TCF4_context |>

  dplyr::group_by(
    ensembl_gene_id,
    gene_symbol,
    TCF4_Cell_Context
  ) |>

  dplyr::summarise(

    Number_of_Fine_Mapped_Transcripts =
      dplyr::n_distinct(
        ensembl_transcript_id
      ),

    Highest_Transcript_PIP = max(
      posterior_inclusion_probability,
      na.rm = TRUE
    ),

    Has_Very_High_PIP_Transcript = any(
      posterior_inclusion_probability >= 0.90
    ),

    .groups = "drop"
  )

# Summarize PIP by cell context
# Compare PIP evidence among the independent TCF4 contexts
s5_TCF4_context_PIP_summary <-
  s5_finemapped_TCF4_gene_PIP_summary |>

  dplyr::group_by(
    TCF4_Cell_Context
  ) |>

  dplyr::summarise(

    Fine_Mapped_Genes = dplyr::n(),

    Median_Highest_PIP = median(
      Highest_Transcript_PIP,
      na.rm = TRUE
    ),

    Mean_Highest_PIP = mean(
      Highest_Transcript_PIP,
      na.rm = TRUE
    ),

    Genes_With_Very_High_PIP = sum(
      Has_Very_High_PIP_Transcript
    ),

    Percent_With_Very_High_PIP = round(
      100 *
        Genes_With_Very_High_PIP /
        Fine_Mapped_Genes,
      digits = 2
    ),

    .groups = "drop"
  )

s5_TCF4_context_PIP_summary


# Reconstruct the Fisher table with retained genes in the first column
s5_cell_specific_retention_table_clear <- matrix(
  c(

    # NPC-only genes
    11,  # Retained
    21,  # Not retained

    # SH-SY5Y-only genes
    27,  # Retained
    35   # Not retained
  ),

  nrow = 2,
  byrow = TRUE
)

rownames(
  s5_cell_specific_retention_table_clear
) <- c(
  "NPC only",
  "SH-SY5Y only"
)

colnames(
  s5_cell_specific_retention_table_clear
) <- c(
  "Retained",
  "Not retained"
)

s5_cell_specific_retention_table_clear

# Run the same two-sided Fisher test
fisher.test(
  s5_cell_specific_retention_table_clear,
  alternative = "two.sided"
)

#########################
# Examine the detailed TCF4 membership groups

# Calculate retention separately for every original membership group
s5_detailed_membership_retention_summary <-
  s5_preFOCUS_TCF4_context |>

  dplyr::group_by(
    Membership_Group
  ) |>

  dplyr::summarise(

    PreFOCUS_Genes = dplyr::n(),

    Retained_Fine_Mapped_Genes = sum(
      Retained_After_FOCUS
    ),

    Not_Retained_Genes = sum(
      !Retained_After_FOCUS
    ),

    Percent_Retained = round(
      100 *
        Retained_Fine_Mapped_Genes /
        PreFOCUS_Genes,
      digits = 2
    ),

    .groups = "drop"
  ) |>

  dplyr::arrange(
    dplyr::desc(
      Percent_Retained
    )
  )

s5_detailed_membership_retention_summary

#########################

# Test whether Forrest-associated genes have different fine-mapping retention

s5_preFOCUS_Forrest_retention_data <-
  s5_preFOCUS_TCF4_context |>

  dplyr::mutate(
    Forrest_Association = dplyr::case_when(

      Membership_Group %in% c(
        "Forrest only",
        "NPC + Forrest",
        "McClay + Forrest",
        "All three"
      ) ~
        "Forrest associated",

      Membership_Group %in% c(
        "NPC only",
        "McClay only",
        "NPC + McClay"
      ) ~
        "Forrest absent",

      TRUE ~
        NA_character_
    )
  )

# Construct the table with retained genes in the first column
s5_Forrest_retention_table <- matrix(
  c(

    # Forrest-associated genes
    22,  # Retained
    21,  # Not retained

    # Genes without Forrest association
    24,  # Retained
    45   # Not retained
  ),

  nrow = 2,
  byrow = TRUE
)

rownames(
  s5_Forrest_retention_table
) <- c(
  "Forrest associated",
  "Forrest absent"
)

colnames(
  s5_Forrest_retention_table
) <- c(
  "Retained",
  "Not retained"
)

s5_Forrest_retention_table

# Run the two-sided Fisher test
s5_Forrest_retention_fisher <- fisher.test(
  s5_Forrest_retention_table,
  alternative = "two.sided"
)

s5_Forrest_retention_fisher
