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


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Import hg38 TCF4 peak files
npc_peaks <- import(
  "NPC_hg38_summit_500bp.bed"
)

mcclay_peaks <- import(
  "McClay_TCF4_11322_consensus_hg38_sorted.bed"
)

forrest_peaks <- import(
  "Forrest_hg38.bed"
)

# Annotate NPC peaks
npc_anno <- annotatePeak(
  npc_peaks,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

npc_df <- as.data.frame(npc_anno)

# Extract unique NPC genes
npc_genes_hg38 <- npc_df %>%
  filter(
    !is.na(SYMBOL)
  ) %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

length(npc_genes_hg38)
head(npc_genes_hg38)
write.csv(npc_genes_hg38, "npc_genes_hg38.csv")

## Annotate McClay peaks
mcclay_anno <- annotatePeak(
  mcclay_peaks,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

mcclay_df <- as.data.frame(mcclay_anno)

mcclay_genes_hg38 <- mcclay_df %>%
  filter(
    !is.na(SYMBOL)
  ) %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

length(mcclay_genes_hg38)

write.csv(mcclay_genes_hg38, "mcclay_genes_hg38.csv")

## Annotate Forrest peaks
forrest_anno <- annotatePeak(
  forrest_peaks,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

forrest_df <- as.data.frame(forrest_anno)

forrest_genes_hg38 <- forrest_df %>%
  filter(
    !is.na(SYMBOL)
  ) %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

length(forrest_genes_hg38)

write.csv(forrest_genes_hg38, "forrest_genes_hg38.csv")

####################################################################
# Generate hg19 TCF4 gene lists for SynGO
####################################################################
# Import hg19 TCF4 peak files
npc_hg19 <- import(
  "NPC_ab21_idr_500bp_summit.bed"
)

forrest_hg19 <- import(
  "Forrest_TCF4_IDR_sorted_hg19.bed"
)

mcclay_hg19_df <- read.csv(
  "McClay_TCF4_11322_consensus_hg19_annotated.csv"
)

## Convert McClay CSV peak table to GRanges

mcclay_hg19 <- makeGRangesFromDataFrame(
  mcclay_hg19_df,
  seqnames.field = "chrom",
  start.field = "start_hg19",
  end.field = "end_hg19",
  keep.extra.columns = TRUE
)

## Annotate NPC peaks

npc_anno <- annotatePeak(
  npc_hg19,
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

npc_anno_df <- as.data.frame(npc_anno)

npc_genes_hg19 <- npc_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

length(npc_genes_hg19)

## Annotate McClay peaks

mcclay_anno <- annotatePeak(
  mcclay_hg19,
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

mcclay_anno_df <- as.data.frame(mcclay_anno)

mcclay_genes_hg19 <- mcclay_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

length(mcclay_genes_hg19)

## Annotate Forrest peaks
forrest_anno <- annotatePeak(
  forrest_hg19,
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

forrest_anno_df <- as.data.frame(forrest_anno)

forrest_genes_hg19 <- forrest_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

length(forrest_genes_hg19)

## Save gene lists for SynGO
write.table(
  npc_genes_hg19,
  "NPC_TCF4_hg19_genes_for_SynGO.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  mcclay_genes_hg19,
  "McClay_TCF4_hg19_genes_for_SynGO.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  forrest_genes_hg19,
  "Forrest_TCF4_hg19_genes_for_SynGO.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)