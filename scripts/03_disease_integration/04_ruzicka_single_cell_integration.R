####################################################################
# Ruzicka single-cell dataset integration with hg38 TCF4 peaks
####################################################################

# Load libraries
library(GenomicRanges)  # handling genomic intervals and peak overlap analysis
library(rtracklayer)    # importing and exporting BED, GFF, GTF, BigWig, and other genomic files
library(ChIPseeker)     # peak annotation
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # Human hg38 transcript annotation database used by ChIPseeker
library(org.Hs.eg.db)  # Human gene annotation database, used to map Entrez IDs to gene symbols
library(clusterProfiler)  # For GO, KEGG, and pathway enrichment analysis
library(txdbmaker)     # chromosome naming style, such as converting 1 to chr1
library(GenomeInfoDb)  # enrichment plots, including dotplot and barplot
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationDbi)
library(GenomeInfoDb)
library(tibble)
library(pheatmap)

txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Import files
npc_peaks     <- import("NPC_hg38_summit_500bp.bed")
mcclay_peaks  <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
forrest_peaks <- import("Forrest_hg38.bed")

# Check peak counts
peak_count_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Peak_Count = c(
    length(npc_peaks),
    length(mcclay_peaks),
    length(forrest_peaks)
  )
)

peak_count_summary

# Import Ruzickas files
s3_raw <- read_xlsx("Single_cell_multi_cohort_paper_datset/science.adg5136_data_s3.xlsx")
s4_raw <- excel_sheets("Single_cell_multi_cohort_paper_datset/science.adg5136_data_s4.xlsx")
s5_raw <- read_xlsx("Single_cell_multi_cohort_paper_datset/science.adg5136_data_s5.xlsx")
s10_raw <-  excel_sheets("Single_cell_multi_cohort_paper_datset/science.adg5136_data_s10.xlsx")
s12_raw <- read_xlsx("Single_cell_multi_cohort_paper_datset/science.adg5136_data_s12.xlsx")

head(s3_raw)
head(s4_raw)
head(s5_raw)
head(s10_raw)
head(s12_raw)


# Annotate hg38 peaks to nearest genes
npc_anno <- annotatePeak(
  npc_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

mcclay_anno <- annotatePeak(
  mcclay_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

forrest_anno <- annotatePeak(
  forrest_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

## Annotate hg38 TCF4 peaks to nearest genes

npc_anno <- annotatePeak(
  npc_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

mcclay_anno <- annotatePeak(
  mcclay_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

forrest_anno <- annotatePeak(
  forrest_peaks,
  TxDb = txdb_hg38,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

# Convert annotation results to data frames
npc_anno_df <- as.data.frame(npc_anno)
mcclay_anno_df <- as.data.frame(mcclay_anno)
forrest_anno_df <- as.data.frame(forrest_anno)

# Extract unique TCF4-associated genes
npc_genes_hg38 <- npc_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

mcclay_genes_hg38 <- mcclay_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

forrest_genes_hg38 <- forrest_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

# Gene count summary
tcf4_gene_count_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Peak_Count = c(
    length(npc_peaks),
    length(mcclay_peaks),
    length(forrest_peaks)
  ),
  Gene_Count = c(
    length(npc_genes_hg38),
    length(mcclay_genes_hg38),
    length(forrest_genes_hg38)
  )
)

tcf4_gene_count_summary

#############################################################################
# Analysis 1: Overlap of hg38 TCF4 genes with Ruzicka S3 cell-type expression
#############################################################################
# Clean S3 table
s3_clean <- s3_raw %>%
  rename(
    gene = ...1
  ) %>%
  filter(
    !is.na(gene),
    gene != ""
  )

# Check S3 dimensions and cell-type columns
dim(s3_clean)
colnames(s3_clean)


## Convert S3 from wide format to long format
s3_long <- s3_clean %>%
  pivot_longer(
    cols = -gene,
    names_to = "Cell_Type",
    values_to = "Expression"
  )

head(s3_long)
View(s3_long)

# Define expressed genes in each cell type
# Here, Expression > 0 means the gene is detected in that Ruzicka cell type.
# Later, we can also repeat this using Expression > 1 for a stricter threshold.
s3_expressed_by_celltype_gt0 <- s3_long %>%
  filter(
    !is.na(Expression),
    Expression > 0
  ) %>%
  distinct(
    Cell_Type,
    gene
  )

# Count expressed genes per Ruzicka cell type
s3_celltype_gene_counts_gt0 <- s3_expressed_by_celltype_gt0 %>%
  group_by(
    Cell_Type
  ) %>%
  summarise(
    Ruzicka_Expressed_Genes_gt0 = n_distinct(gene),
    .groups = "drop"
  )

s3_celltype_gene_counts_gt0
View(s3_celltype_gene_counts_gt0)

# NPC TCF4 overlap with each Ruzicka cell type
npc_s3_celltype_overlap_gt0 <- s3_expressed_by_celltype_gt0 %>%
  group_by(
    Cell_Type
  ) %>%
  summarise(
    NPC_TCF4_Expressed_Genes = length(
      intersect(
        gene,
        npc_genes_hg38
      )
    ),
    .groups = "drop"
  )


# McClay TCF4 overlap with each Ruzicka cell type
mcclay_s3_celltype_overlap_gt0 <- s3_expressed_by_celltype_gt0 %>%
  group_by(
    Cell_Type
  ) %>%
  summarise(
    McClay_TCF4_Expressed_Genes = length(
      intersect(
        gene,
        mcclay_genes_hg38
      )
    ),
    .groups = "drop"
  )




## Forrest TCF4 overlap with each Ruzicka cell type
forrest_s3_celltype_overlap_gt0 <- s3_expressed_by_celltype_gt0 %>%
  group_by(
    Cell_Type
  ) %>%
  summarise(
    Forrest_TCF4_Expressed_Genes = length(
      intersect(
        gene,
        forrest_genes_hg38
      )
    ),
    .groups = "drop"
  )


# Combine all overlap results into one table
s3_tcf4_celltype_overlap_gt0 <- s3_celltype_gene_counts_gt0 %>%
  left_join(
    npc_s3_celltype_overlap_gt0,
    by = "Cell_Type"
  ) %>%
  left_join(
    mcclay_s3_celltype_overlap_gt0,
    by = "Cell_Type"
  ) %>%
  left_join(
    forrest_s3_celltype_overlap_gt0,
    by = "Cell_Type"
  )

# Add percentages
## These percentages answer:
## Among all TCF4-associated genes in each dataset, what percent are expressed
## in each Ruzicka cell type?
s3_tcf4_celltype_overlap_gt0 <- s3_tcf4_celltype_overlap_gt0 %>%
  mutate(
    NPC_Percent_of_TCF4_Genes = round(
      NPC_TCF4_Expressed_Genes / length(npc_genes_hg38) * 100,
      3
    ),
    McClay_Percent_of_TCF4_Genes = round(
      McClay_TCF4_Expressed_Genes / length(mcclay_genes_hg38) * 100,
      3
    ),
    Forrest_Percent_of_TCF4_Genes = round(
      Forrest_TCF4_Expressed_Genes / length(forrest_genes_hg38) * 100,
      3
    )
  )

# View final table ranked by Forrest overlap
s3_tcf4_celltype_overlap_gt0 <- s3_tcf4_celltype_overlap_gt0 %>%
  arrange(
    desc(Forrest_TCF4_Expressed_Genes)
  )

s3_tcf4_celltype_overlap_gt0

# Clean display table
s3_tcf4_celltype_overlap_display <- s3_tcf4_celltype_overlap_gt0 %>%
  select(
    Cell_Type,
    Ruzicka_Expressed_Genes_gt0,
    NPC_TCF4_Expressed_Genes,
    NPC_Percent_of_TCF4_Genes,
    McClay_TCF4_Expressed_Genes,
    McClay_Percent_of_TCF4_Genes,
    Forrest_TCF4_Expressed_Genes,
    Forrest_Percent_of_TCF4_Genes
  )

s3_tcf4_celltype_overlap_display
View(s3_tcf4_celltype_overlap_display)

# Save S3 cell-type overlap table
write.csv(
  s3_tcf4_celltype_overlap_display,
  "Ruzicka_S3_TCF4_celltype_expression_overlap_gt0.csv",
  row.names = FALSE
)

# Make heatmap input table
s3_heatmap_df <- s3_tcf4_celltype_overlap_display %>%
  select(
    Cell_Type,
    NPC_Percent_of_TCF4_Genes,
    McClay_Percent_of_TCF4_Genes,
    Forrest_Percent_of_TCF4_Genes
  ) %>%
  column_to_rownames(
    var = "Cell_Type"
  )

# Convert to matrix
s3_heatmap_matrix <- as.matrix(
  s3_heatmap_df
)

# Rename columns for cleaner plot
colnames(s3_heatmap_matrix) <- c(
  "NPC",
  "McClay",
  "Forrest"
)

# Draw heatmap
pheatmap(
  s3_heatmap_matrix,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.1f",
  main = "TCF4 Target Gene Expression Across Ruzicka Cell Types",
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 45
)

######################################################################
# Analysis 2: Fisher enrichment of TCF4 genes in Ruzicka S3 cell types
######################################################################

# Make background from all unique genes in S3
s3_background_genes <- s3_clean %>%
  filter(
    !is.na(gene),
    gene != ""
  ) %>%
  pull(gene) %>%
  unique() %>%
  sort()

length(s3_background_genes)


# Keep only TCF4 genes that are present in the S3 background
npc_genes_bg <- intersect(
  npc_genes_hg38,
  s3_background_genes
)

mcclay_genes_bg <- intersect(
  mcclay_genes_hg38,
  s3_background_genes
)

forrest_genes_bg <- intersect(
  forrest_genes_hg38,
  s3_background_genes
)

length(npc_genes_bg)
length(mcclay_genes_bg)
length(forrest_genes_bg)


# Create an empty table to store Fisher results
s3_fisher_enrichment_results <- data.frame()


# Loop through each Ruzicka cell type
for (celltype in unique(s3_expressed_by_celltype_gt0$Cell_Type)) {

  # Genes expressed in this cell type

  celltype_genes <- s3_expressed_by_celltype_gt0 %>%
    filter(
      Cell_Type == celltype
    ) %>%
    pull(gene) %>%
    unique()

  ## Restrict cell-type genes to S3 background

  celltype_genes <- intersect(
    celltype_genes,
    s3_background_genes
  )


  # NPC Fisher test

  npc_a <- length(intersect(npc_genes_bg, celltype_genes))
  npc_b <- length(setdiff(npc_genes_bg, celltype_genes))
  npc_c <- length(setdiff(celltype_genes, npc_genes_bg))
  npc_d <- length(setdiff(s3_background_genes, union(npc_genes_bg, celltype_genes)))

  npc_table <- matrix(
    c(npc_a, npc_b, npc_c, npc_d),
    nrow = 2,
    byrow = TRUE
  )

  npc_fisher <- fisher.test(
    npc_table,
    alternative = "greater"
  )


  # McClay Fisher test

  mcclay_a <- length(intersect(mcclay_genes_bg, celltype_genes))
  mcclay_b <- length(setdiff(mcclay_genes_bg, celltype_genes))
  mcclay_c <- length(setdiff(celltype_genes, mcclay_genes_bg))
  mcclay_d <- length(setdiff(s3_background_genes, union(mcclay_genes_bg, celltype_genes)))

  mcclay_table <- matrix(
    c(mcclay_a, mcclay_b, mcclay_c, mcclay_d),
    nrow = 2,
    byrow = TRUE
  )

  mcclay_fisher <- fisher.test(
    mcclay_table,
    alternative = "greater"
  )


  # Forrest Fisher test

  forrest_a <- length(intersect(forrest_genes_bg, celltype_genes))
  forrest_b <- length(setdiff(forrest_genes_bg, celltype_genes))
  forrest_c <- length(setdiff(celltype_genes, forrest_genes_bg))
  forrest_d <- length(setdiff(s3_background_genes, union(forrest_genes_bg, celltype_genes)))

  forrest_table <- matrix(
    c(forrest_a, forrest_b, forrest_c, forrest_d),
    nrow = 2,
    byrow = TRUE
  )

  forrest_fisher <- fisher.test(
    forrest_table,
    alternative = "greater"
  )


  # Add results for this cell type

  s3_fisher_enrichment_results <- rbind(
    s3_fisher_enrichment_results,
    data.frame(
      Cell_Type = celltype,
      Background_Genes = length(s3_background_genes),
      Celltype_Expressed_Genes = length(celltype_genes),

      NPC_TCF4_Background_Genes = length(npc_genes_bg),
      NPC_Overlap = npc_a,
      NPC_Odds_Ratio = as.numeric(npc_fisher$estimate),
      NPC_P_Value = npc_fisher$p.value,

      McClay_TCF4_Background_Genes = length(mcclay_genes_bg),
      McClay_Overlap = mcclay_a,
      McClay_Odds_Ratio = as.numeric(mcclay_fisher$estimate),
      McClay_P_Value = mcclay_fisher$p.value,

      Forrest_TCF4_Background_Genes = length(forrest_genes_bg),
      Forrest_Overlap = forrest_a,
      Forrest_Odds_Ratio = as.numeric(forrest_fisher$estimate),
      Forrest_P_Value = forrest_fisher$p.value
    )
  )
}


# Adjust p-values for multiple testing
s3_fisher_enrichment_results$NPC_FDR <- p.adjust(
  s3_fisher_enrichment_results$NPC_P_Value,
  method = "BH"
)

s3_fisher_enrichment_results$McClay_FDR <- p.adjust(
  s3_fisher_enrichment_results$McClay_P_Value,
  method = "BH"
)

s3_fisher_enrichment_results$Forrest_FDR <- p.adjust(
  s3_fisher_enrichment_results$Forrest_P_Value,
  method = "BH"
)


# Sort by strongest Forrest enrichment
s3_fisher_enrichment_results <- s3_fisher_enrichment_results %>%
  arrange(
    desc(Forrest_Odds_Ratio)
  )

s3_fisher_enrichment_results

##################################################################
# Analysis 3 corrected: TCF4 overlap with Ruzicka S4 DE genes
##################################################################
s4_sheets <- s4_raw

s4_enrichment_results_corrected <- data.frame()

for(sheet_name in s4_sheets){

  cat("Processing:", sheet_name, "\n")

  ## Read one cell-type sheet
  s4_tmp <- read_xlsx(
    "Single_cell_multi_cohort_paper_datset/science.adg5136_data_s4.xlsx",
    sheet = sheet_name
  )

  ## Background = all genes tested in this specific cell type
  background_genes <- s4_tmp %>%
    filter(
      !is.na(gene),
      gene != ""
    ) %>%
    pull(gene) %>%
    unique()

  ## Significant SCZ DE genes in this cell type
  de_genes <- s4_tmp %>%
    filter(
      !is.na(gene),
      gene != "",
      Meta_adj.P.Val < 0.05
    ) %>%
    pull(gene) %>%
    unique()

  ## Restrict TCF4 genes to this cell-type background
  npc_tcf4_bg <- intersect(
    npc_genes_hg38,
    background_genes
  )

  mcclay_tcf4_bg <- intersect(
    mcclay_genes_hg38,
    background_genes
  )

  forrest_tcf4_bg <- intersect(
    forrest_genes_hg38,
    background_genes
  )

  ## NPC Fisher test
  npc_overlap <- intersect(
    npc_tcf4_bg,
    de_genes
  )

  npc_a <- length(npc_overlap)
  npc_b <- length(npc_tcf4_bg) - npc_a
  npc_c <- length(de_genes) - npc_a
  npc_d <- length(background_genes) - npc_a - npc_b - npc_c

  npc_matrix <- matrix(
    c(
      npc_a, npc_b,
      npc_c, npc_d
    ),
    nrow = 2,
    byrow = TRUE
  )

  npc_fisher <- fisher.test(
    npc_matrix,
    alternative = "greater"
  )

  ## McClay Fisher test
  mcclay_overlap <- intersect(
    mcclay_tcf4_bg,
    de_genes
  )

  mcclay_a <- length(mcclay_overlap)
  mcclay_b <- length(mcclay_tcf4_bg) - mcclay_a
  mcclay_c <- length(de_genes) - mcclay_a
  mcclay_d <- length(background_genes) - mcclay_a - mcclay_b - mcclay_c

  mcclay_matrix <- matrix(
    c(
      mcclay_a, mcclay_b,
      mcclay_c, mcclay_d
    ),
    nrow = 2,
    byrow = TRUE
  )

  mcclay_fisher <- fisher.test(
    mcclay_matrix,
    alternative = "greater"
  )

  ## Forrest Fisher test
  forrest_overlap <- intersect(
    forrest_tcf4_bg,
    de_genes
  )

  forrest_a <- length(forrest_overlap)
  forrest_b <- length(forrest_tcf4_bg) - forrest_a
  forrest_c <- length(de_genes) - forrest_a
  forrest_d <- length(background_genes) - forrest_a - forrest_b - forrest_c

  forrest_matrix <- matrix(
    c(
      forrest_a, forrest_b,
      forrest_c, forrest_d
    ),
    nrow = 2,
    byrow = TRUE
  )

  forrest_fisher <- fisher.test(
    forrest_matrix,
    alternative = "greater"
  )

  ## Save results
  s4_enrichment_results_corrected <- rbind(
    s4_enrichment_results_corrected,
    data.frame(
      Cell_Type = sheet_name,

      Background_Genes = length(background_genes),
      DE_Genes = length(de_genes),

      NPC_TCF4_Background_Genes = length(npc_tcf4_bg),
      NPC_Overlap = npc_a,
      NPC_OR = as.numeric(npc_fisher$estimate),
      NPC_P = npc_fisher$p.value,

      McClay_TCF4_Background_Genes = length(mcclay_tcf4_bg),
      McClay_Overlap = mcclay_a,
      McClay_OR = as.numeric(mcclay_fisher$estimate),
      McClay_P = mcclay_fisher$p.value,

      Forrest_TCF4_Background_Genes = length(forrest_tcf4_bg),
      Forrest_Overlap = forrest_a,
      Forrest_OR = as.numeric(forrest_fisher$estimate),
      Forrest_P = forrest_fisher$p.value
    )
  )
}

## FDR correction

s4_enrichment_results_corrected$NPC_FDR <- p.adjust(
  s4_enrichment_results_corrected$NPC_P,
  method = "BH"
)

s4_enrichment_results_corrected$McClay_FDR <- p.adjust(
  s4_enrichment_results_corrected$McClay_P,
  method = "BH"
)

s4_enrichment_results_corrected$Forrest_FDR <- p.adjust(
  s4_enrichment_results_corrected$Forrest_P,
  method = "BH"
)

## Sort and view

s4_enrichment_results_corrected <- s4_enrichment_results_corrected %>%
  arrange(
    desc(NPC_OR)
  )

s4_enrichment_results_corrected

View(s4_enrichment_results_corrected)

##################################################################
# Analysis 4: TCF4 enrichment in S4 SCZ-up and SCZ-down DE genes
##################################################################
s4_sheets <- s4_raw

s4_directional_enrichment_results <- data.frame()

for(sheet_name in s4_sheets){

  cat("Processing:", sheet_name, "\n")

  # Read one S4 cell-type sheet
  s4_tmp <- read_xlsx(
    "Single_cell_multi_cohort_paper_datset/science.adg5136_data_s4.xlsx",
    sheet = sheet_name
  )

  # Background = all genes tested in this specific cell type
  background_genes <- s4_tmp %>%
    filter(
      !is.na(gene),
      gene != ""
    ) %>%
    pull(gene) %>%
    unique()

  # Restrict TCF4 genes to this cell-type background
  npc_tcf4_bg <- intersect(
    npc_genes_hg38,
    background_genes
  )

  mcclay_tcf4_bg <- intersect(
    mcclay_genes_hg38,
    background_genes
  )

  forrest_tcf4_bg <- intersect(
    forrest_genes_hg38,
    background_genes
  )

  # Define SCZ-upregulated genes
  up_genes <- s4_tmp %>%
    filter(
      !is.na(gene),
      gene != "",
      Meta_adj.P.Val < 0.05,
      Meta_logFC > 0
    ) %>%
    pull(gene) %>%
    unique()

  # Define SCZ-downregulated genes
  down_genes <- s4_tmp %>%
    filter(
      !is.na(gene),
      gene != "",
      Meta_adj.P.Val < 0.05,
      Meta_logFC < 0
    ) %>%
    pull(gene) %>%
    unique()

  # Run analysis separately for UP and DOWN gene sets
  for(direction in c("Upregulated", "Downregulated")){

    if(direction == "Upregulated"){
      de_genes <- up_genes
    }

    if(direction == "Downregulated"){
      de_genes <- down_genes
    }

    ##################################################
    ## NPC Fisher test
    ##################################################
    npc_overlap <- intersect(
      npc_tcf4_bg,
      de_genes
    )

    npc_a <- length(npc_overlap)
    npc_b <- length(npc_tcf4_bg) - npc_a
    npc_c <- length(de_genes) - npc_a
    npc_d <- length(background_genes) - npc_a - npc_b - npc_c

    npc_matrix <- matrix(
      c(
        npc_a, npc_b,
        npc_c, npc_d
      ),
      nrow = 2,
      byrow = TRUE
    )

    npc_fisher <- fisher.test(
      npc_matrix,
      alternative = "greater"
    )

    ##################################################
    ## McClay Fisher test
    ##################################################
    mcclay_overlap <- intersect(
      mcclay_tcf4_bg,
      de_genes
    )

    mcclay_a <- length(mcclay_overlap)
    mcclay_b <- length(mcclay_tcf4_bg) - mcclay_a
    mcclay_c <- length(de_genes) - mcclay_a
    mcclay_d <- length(background_genes) - mcclay_a - mcclay_b - mcclay_c

    mcclay_matrix <- matrix(
      c(
        mcclay_a, mcclay_b,
        mcclay_c, mcclay_d
      ),
      nrow = 2,
      byrow = TRUE
    )

    mcclay_fisher <- fisher.test(
      mcclay_matrix,
      alternative = "greater"
    )

    ##################################################
    ## Forrest Fisher test
    ##################################################
    forrest_overlap <- intersect(
      forrest_tcf4_bg,
      de_genes
    )

    forrest_a <- length(forrest_overlap)
    forrest_b <- length(forrest_tcf4_bg) - forrest_a
    forrest_c <- length(de_genes) - forrest_a
    forrest_d <- length(background_genes) - forrest_a - forrest_b - forrest_c

    forrest_matrix <- matrix(
      c(
        forrest_a, forrest_b,
        forrest_c, forrest_d
      ),
      nrow = 2,
      byrow = TRUE
    )

    forrest_fisher <- fisher.test(
      forrest_matrix,
      alternative = "greater"
    )

    ##################################################
    ## Save results
    ##################################################
    s4_directional_enrichment_results <- rbind(
      s4_directional_enrichment_results,
      data.frame(
        Cell_Type = sheet_name,
        Direction = direction,

        Background_Genes = length(background_genes),
        DE_Genes = length(de_genes),

        NPC_TCF4_Background_Genes = length(npc_tcf4_bg),
        NPC_Overlap = npc_a,
        NPC_OR = as.numeric(npc_fisher$estimate),
        NPC_P = npc_fisher$p.value,

        McClay_TCF4_Background_Genes = length(mcclay_tcf4_bg),
        McClay_Overlap = mcclay_a,
        McClay_OR = as.numeric(mcclay_fisher$estimate),
        McClay_P = mcclay_fisher$p.value,

        Forrest_TCF4_Background_Genes = length(forrest_tcf4_bg),
        Forrest_Overlap = forrest_a,
        Forrest_OR = as.numeric(forrest_fisher$estimate),
        Forrest_P = forrest_fisher$p.value
      )
    )
  }
}

# FDR correction separately for each TCF4 dataset

s4_directional_enrichment_results$NPC_FDR <- p.adjust(
  s4_directional_enrichment_results$NPC_P,
  method = "BH"
)

s4_directional_enrichment_results$McClay_FDR <- p.adjust(
  s4_directional_enrichment_results$McClay_P,
  method = "BH"
)

s4_directional_enrichment_results$Forrest_FDR <- p.adjust(
  s4_directional_enrichment_results$Forrest_P,
  method = "BH"
)

# View final directional enrichment table
s4_directional_enrichment_results <- s4_directional_enrichment_results %>%
  arrange(
    Direction,
    Forrest_FDR
  )

s4_directional_enrichment_results

# Significant directional enrichments only
s4_directional_sig_results <- s4_directional_enrichment_results %>%
  filter(
    NPC_FDR < 0.05 |
      McClay_FDR < 0.05 |
      Forrest_FDR < 0.05
  )

s4_directional_sig_results

##################################################################
# Analysis 5: TCF4 overlap with Ruzicka S5 reproducible DE genes
##################################################################
## S5 already contains reproducible DE genes with a celltype column
## We will compare S5 reproducible genes with NPC, McClay, and Forrest TCF4 genes

## Check S5 structure
head(s5_raw)
colnames(s5_raw)
table(s5_raw$celltype)


## Clean S5 reproducible gene table
s5_clean <- s5_raw %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(celltype),
    celltype != ""
  ) %>%
  mutate(
    gene = trimws(gene),
    celltype = trimws(celltype)
  )

length(unique(s5_clean$gene))
length(unique(s5_clean$celltype))


## Create empty result table
s5_tcf4_enrichment_results <- data.frame()


## Loop through each S5 cell type
for(cell_type_name in unique(s5_clean$celltype)){

  cat("Processing:", cell_type_name, "\n")

  ## S5 reproducible DE genes for this cell type

  s5_cell_genes <- s5_clean %>%
    filter(
      celltype == cell_type_name
    ) %>%
    pull(gene) %>%
    unique()


  ## Use the matching S4 sheet as background
  ## Background = all genes tested for DE in that cell type

  s4_background_tmp <- read_xlsx(
    "Single_cell_multi_cohort_paper_datset/science.adg5136_data_s4.xlsx",
    sheet = cell_type_name
  )

  background_genes <- s4_background_tmp %>%
    filter(
      !is.na(gene),
      gene != ""
    ) %>%
    pull(gene) %>%
    unique()


  ## Restrict S5 genes to background

  s5_cell_genes <- intersect(
    s5_cell_genes,
    background_genes
  )


  ## Restrict TCF4 genes to this cell-type background

  npc_tcf4_bg <- intersect(
    npc_genes_hg38,
    background_genes
  )

  mcclay_tcf4_bg <- intersect(
    mcclay_genes_hg38,
    background_genes
  )

  forrest_tcf4_bg <- intersect(
    forrest_genes_hg38,
    background_genes
  )


  ## NPC Fisher test

  npc_overlap <- intersect(
    npc_tcf4_bg,
    s5_cell_genes
  )

  npc_a <- length(npc_overlap)
  npc_b <- length(npc_tcf4_bg) - npc_a
  npc_c <- length(s5_cell_genes) - npc_a
  npc_d <- length(background_genes) - npc_a - npc_b - npc_c

  npc_matrix <- matrix(
    c(
      npc_a, npc_b,
      npc_c, npc_d
    ),
    nrow = 2,
    byrow = TRUE
  )

  npc_fisher <- fisher.test(
    npc_matrix,
    alternative = "greater"
  )


  ## McClay Fisher test

  mcclay_overlap <- intersect(
    mcclay_tcf4_bg,
    s5_cell_genes
  )

  mcclay_a <- length(mcclay_overlap)
  mcclay_b <- length(mcclay_tcf4_bg) - mcclay_a
  mcclay_c <- length(s5_cell_genes) - mcclay_a
  mcclay_d <- length(background_genes) - mcclay_a - mcclay_b - mcclay_c

  mcclay_matrix <- matrix(
    c(
      mcclay_a, mcclay_b,
      mcclay_c, mcclay_d
    ),
    nrow = 2,
    byrow = TRUE
  )

  mcclay_fisher <- fisher.test(
    mcclay_matrix,
    alternative = "greater"
  )


  ## Forrest Fisher test

  forrest_overlap <- intersect(
    forrest_tcf4_bg,
    s5_cell_genes
  )

  forrest_a <- length(forrest_overlap)
  forrest_b <- length(forrest_tcf4_bg) - forrest_a
  forrest_c <- length(s5_cell_genes) - forrest_a
  forrest_d <- length(background_genes) - forrest_a - forrest_b - forrest_c

  forrest_matrix <- matrix(
    c(
      forrest_a, forrest_b,
      forrest_c, forrest_d
    ),
    nrow = 2,
    byrow = TRUE
  )

  forrest_fisher <- fisher.test(
    forrest_matrix,
    alternative = "greater"
  )


  ## Save results

  s5_tcf4_enrichment_results <- rbind(
    s5_tcf4_enrichment_results,
    data.frame(
      Cell_Type = cell_type_name,

      Background_Genes = length(background_genes),
      S5_Reproducible_Genes = length(s5_cell_genes),

      NPC_TCF4_Background_Genes = length(npc_tcf4_bg),
      NPC_Overlap = npc_a,
      NPC_OR = as.numeric(npc_fisher$estimate),
      NPC_P = npc_fisher$p.value,

      McClay_TCF4_Background_Genes = length(mcclay_tcf4_bg),
      McClay_Overlap = mcclay_a,
      McClay_OR = as.numeric(mcclay_fisher$estimate),
      McClay_P = mcclay_fisher$p.value,

      Forrest_TCF4_Background_Genes = length(forrest_tcf4_bg),
      Forrest_Overlap = forrest_a,
      Forrest_OR = as.numeric(forrest_fisher$estimate),
      Forrest_P = forrest_fisher$p.value
    )
  )
}


## FDR correction
s5_tcf4_enrichment_results$NPC_FDR <- p.adjust(
  s5_tcf4_enrichment_results$NPC_P,
  method = "BH"
)

s5_tcf4_enrichment_results$McClay_FDR <- p.adjust(
  s5_tcf4_enrichment_results$McClay_P,
  method = "BH"
)

s5_tcf4_enrichment_results$Forrest_FDR <- p.adjust(
  s5_tcf4_enrichment_results$Forrest_P,
  method = "BH"
)


## Sort by Forrest enrichment
s5_tcf4_enrichment_results <- s5_tcf4_enrichment_results %>%
  arrange(
    desc(Forrest_OR)
  )

s5_tcf4_enrichment_results


## Significant S5 enrichment results
s5_tcf4_sig_results <- s5_tcf4_enrichment_results %>%
  filter(
    NPC_FDR < 0.05 |
      McClay_FDR < 0.05 |
      Forrest_FDR < 0.05
  )

s5_tcf4_sig_results

##################################################################
# Analysis 6: TCF4 overlap with Ruzicka S10 CUT&Tag binding peaks
# from human PFC neurons
##################################################################
# Gene -level comparison
## Read the TCF4 sheet from S10
s10_tcf4_raw <- read_xlsx(
  "Single_cell_multi_cohort_paper_datset/science.adg5136_data_s10.xlsx",
  sheet = "TCF4.reproducible.bind"
)

head(s10_tcf4_raw)
colnames(s10_tcf4_raw)
dim(s10_tcf4_raw)

## Extract Ruzicka TCF4-associated genes
s10_tcf4_genes <- s10_tcf4_raw %>%
  filter(
    !is.na(geneSymbol),
    geneSymbol != ""
  ) %>%
  pull(geneSymbol) %>%
  unique() %>%
  sort()

length(s10_tcf4_genes)
head(s10_tcf4_genes)

## Compare S10 TCF4 genes with your TCF4 gene sets
npc_s10_tcf4_genes <- intersect(
  npc_genes_hg38,
  s10_tcf4_genes
)

mcclay_s10_tcf4_genes <- intersect(
  mcclay_genes_hg38,
  s10_tcf4_genes
)

forrest_s10_tcf4_genes <- intersect(
  forrest_genes_hg38,
  s10_tcf4_genes
)

s10_gene_overlap_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Your_TCF4_Genes = c(
    length(npc_genes_hg38),
    length(mcclay_genes_hg38),
    length(forrest_genes_hg38)
  ),
  Ruzicka_S10_TCF4_Genes = length(s10_tcf4_genes),
  Overlap_Genes = c(
    length(npc_s10_tcf4_genes),
    length(mcclay_s10_tcf4_genes),
    length(forrest_s10_tcf4_genes)
  )
)

s10_gene_overlap_summary$Percent_of_Ruzicka_S10_Genes <- round(
  s10_gene_overlap_summary$Overlap_Genes /
    s10_gene_overlap_summary$Ruzicka_S10_TCF4_Genes * 100,
  3
)

s10_gene_overlap_summary

# Peak level comparison
# Check if we have peak coordinate or gene coordinates
head(s10_tcf4_raw[, c("seqnames","start","end","geneSymbol","distanceToTSS")])
summary(s10_tcf4_raw$width)

# Create S10 peaks
s10_tcf4_peaks <- GRanges(
  seqnames = s10_tcf4_raw$seqnames,
  ranges = IRanges(
    start = s10_tcf4_raw$start,
    end = s10_tcf4_raw$end
  )
)

length(s10_tcf4_peaks)

# Comapre with our TCF4 peaks
# McClay
mcclay_hits <- findOverlaps(
  mcclay_peaks,
  s10_tcf4_peaks,
  ignore.strand = TRUE
)

length(unique(queryHits(mcclay_hits)))

# NPC
npc_hits <- findOverlaps(
  npc_peaks,
  s10_tcf4_peaks,
  ignore.strand = TRUE
)

length(unique(queryHits(npc_hits)))

# Forrest
forrest_hits <- findOverlaps(
  forrest_peaks,
  s10_tcf4_peaks,
  ignore.strand = TRUE
)

length(unique(queryHits(forrest_hits)))

# Summarize it
s10_peak_overlap_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Total_Peaks = c(
    length(npc_peaks),
    length(mcclay_peaks),
    length(forrest_peaks)
  ),
  Overlapping_Peaks = c(
    length(unique(queryHits(npc_hits))),
    length(unique(queryHits(mcclay_hits))),
    length(unique(queryHits(forrest_hits)))
  )
)

s10_peak_overlap_summary$Percent_Overlap <-
  round(
    s10_peak_overlap_summary$Overlapping_Peaks /
      s10_peak_overlap_summary$Total_Peaks * 100,
    2
  )

s10_peak_overlap_summary