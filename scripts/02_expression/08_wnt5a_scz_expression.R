########################################################################
# Gene expression in patient-derived neural progenitors implicates WNT5A
# signaling in the etiology of schizophrenia
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10947993/
# This table comes from:
# NIHMS1551527-supplement-3.xlsx
# Sheet: Table S3. SCZ vs CTL
# Goal:
# 1. Import the SCZ vs Control differential expression table
# 2. Clean the gene symbols
# 3. Define background, significant, upregulated, and downregulated genes
# 4. Prepare for overlap analysis with TCF4 peak-associated gene sets
########################################################################
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
library(DESeq2)
library(limma)
library(ggplot2)


# Import data
wnt5a_raw <- read_xlsx("NIHMS1551527-supplement-3.xlsx",
                       sheet = "Table S3. SCZ vs CTL",
                       skip = 1)

npc_genes_hg38_df <- read.csv(
  "npc_genes_hg38.csv")

mcclay_genes_hg38_df <- read.csv(
  "mcclay_genes_hg38.csv")

forrest_genes_hg38_df <- read.csv(
  "forrest_genes_hg38.csv")

# Check datasets (Sanity check)
length(rownames(wnt5a_raw))
length(rownames(npc_genes_hg38_df))
length(rownames(mcclay_genes_hg38_df))
length(rownames(forrest_genes_hg38_df))
colnames(wnt5a_raw)
head(wnt5a_raw)
dim(wnt5a_raw)

# Clean datasets
## Remove rows without gene symbols
wnt5a_clean <-  wnt5a_raw %>%
  filter(
    !is.na(Symbol)
  )

## Remove duplicated gene symbols if any
wnt5a_clean <-  wnt5a_clean %>%
  arrange(padj) %>%
  distinct(
    Symbol,
    .keep_all = TRUE
  )

## Check clean table
dim(wnt5a_clean)
head(wnt5a_clean)

# Create gene lists
## Create background gene list
wnt5a_bckground <- unique(wnt5a_clean$Symbol)
length(wnt5a_bckground)

## Create significant gene list
## paper used FDR < 10% but we use FDR < 5%
wnt5a_sig_genes_FDR5 <- unique(
  wnt5a_clean$Symbol[
    wnt5a_clean$padj<0.05
  ]
)

length(wnt5a_sig_genes_FDR5)

## Since number of sig genes with FDR < 5% is low, we choose FDR < 10%
wnt5a_sig_genes_FDR10 <- unique(
  wnt5a_clean$Symbol[
    wnt5a_clean$padj< 0.10
  ]
)

length(wnt5a_sig_genes_FDR10)

##  Find Upregulated in SCZ
wnt5a_up_genes <- unique(
  wnt5a_clean$Symbol[
    wnt5a_clean$padj < 0.10 &
      wnt5a_clean$log2FoldChange > 0
  ]
)

## Downregulated in schizophrenia
wnt5a_down_genes <- unique(
  wnt5a_clean$Symbol[
        wnt5a_clean$padj < 0.10 &
            wnt5a_clean$log2FoldChange < 0
  ]
)

## Check WNT5A gene list size
wnt5a_gene_summary <- data.frame(
  Category = c (
    "Background",
    "Significant_FDR_5%",
    "Significant_FDR_10%",
    "Up_in_SCZ",
    "Down_in_SCZ"
  ),
  Gene_Count = c(
    length(wnt5a_bckground),
    length(wnt5a_sig_genes_FDR5),
    length(wnt5a_sig_genes_FDR10),
    length(wnt5a_up_genes),
    length(wnt5a_down_genes)
  )
)

wnt5a_gene_summary

# Extract gene symbols from TCF4 datasets
tcf4_npc_genes <- unique(
  na.omit(npc_genes_hg38_df$x)
  )

tcf4_mcclay_genes <- unique(
  na.omit(mcclay_genes_hg38_df$x)
)

tcf4_forrest_genes <- unique(
  na.omit(forrest_genes_hg38_df$x)
)

## Check gene list sizes
data.frame(
  TCF4_Dataset = c(
    "NPC_peaks",
    "McClay_peaks",
    "Forrest_peaks"
  ),
  Gene_Count = c(
    length(tcf4_npc_genes),
    length(tcf4_mcclay_genes),
    length(tcf4_forrest_genes)
  )
)

# Restrict TCF4 gene sets to WNT5A paper background
# Fisher test should use the same universe/background
# Here the universe is all genes tested in the WNT5A SCZ vs Control analysis
tcf4_npc_genes_wnt5a_bg <- intersect(
  tcf4_npc_genes,
  wnt5a_bckground
)

tcf4_mcclay_genes_wnt5a_bg <- intersect(
  tcf4_mcclay_genes,
  wnt5a_bckground
)

tcf4_forrest_genes_wnt5a_bg <- intersect(
  tcf4_forrest_genes,
  wnt5a_bckground
)

## Check how many TCF4 genes are present in the WNT5A background
tcf4_wnt5a_background_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC_peaks",
    "McClay_peaks",
    "Forrest_peaks"
  ),
  TCF4_Genes_in_WNT5A_Background = c(
    length(tcf4_npc_genes_wnt5a_bg),
    length(tcf4_mcclay_genes_wnt5a_bg),
    length(tcf4_forrest_genes_wnt5a_bg)
  )
)

tcf4_wnt5a_background_summary

# Raw overlaps with WNT5A significant genes
npc_wnt5a_sig_overlap <- intersect(
  tcf4_npc_genes_wnt5a_bg,
  wnt5a_sig_genes_FDR10
)

mcclay_wnt5a_sig_overlap <- intersect(
  tcf4_mcclay_genes_wnt5a_bg,
  wnt5a_sig_genes_FDR10
)

forrest_wnt5a_sig_overlap <- intersect(
  tcf4_forrest_genes_wnt5a_bg,
  wnt5a_sig_genes_FDR10
)

# Check overlap sizes
wnt5a_raw_overlap_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC_peaks",
    "McClay_peaks",
    "Forrest_peaks"
  ),

  Overlap_with_WNT5A_DEGs = c(
    length(npc_wnt5a_sig_overlap),
    length(mcclay_wnt5a_sig_overlap),
    length(forrest_wnt5a_sig_overlap)
  )
)

wnt5a_raw_overlap_summary
