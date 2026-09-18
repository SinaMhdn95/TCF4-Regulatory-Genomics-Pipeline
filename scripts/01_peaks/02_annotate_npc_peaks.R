# I repeated MACS2 peak calling and IDR analysis using the same criteria as the Forrest study
# and obtained reproducible peak files. MEME-ChIP analysis was then performed, and
# the next step is peak annotation using ChIPseeker.

# Load packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(ggplot2)

# Set TxDb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene  # Loads the hg19 human gene annotation database

# Read IDR peak file
peak_file <- "NPC_ab21_TCF4_IDR_sorted_ver2_ForrestCriterea.bed"

# Converts the BED file into a format ChIPseeker can use
peaks <- readPeakFile(peak_file)

# Annotate peaks (Assigns peaks to promoter, exon, intron, intergenic, etc.)
peak_anno <- annotatePeak(
  peaks,
  tssRegion = c(-3000, 3000),   # Defines promoter as ±3 kb around TSS
  TxDb = txdb,                  # Uses hg19 gene coordinates
  annoDb = "org.Hs.eg.db"
)

# Convert to dataframe (make a table)
peak_df <- as.data.frame(peak_anno)

if (FALSE) {
# Save annotated peaks
write_csv(
  peak_df,
  "TCF4_NPC_ChIPseeker_annotation_new_ForrestCrit.csv"
)
}

# Peak annotation percentages
annotation_summary <- peak_df %>%
  count(annotation, name = "Count") %>%     # Counts how many peaks fall into each region
  mutate(Percent = round(100 * Count / sum(Count), 2)) %>%  # Converts counts into percentages
  arrange(desc(Count))

if (FALSE) {
# Save percentages
write_csv(
  annotation_summary,
  "TCF4_NPC_annotation_percentages.csv"
)
}

# Print summary
annotation_summary

# Annotation plots
plotAnnoPie(peak_anno)   # Makes a pie chart of peak annotations

plotAnnoBar(peak_anno)   # Makes a bar plot of peak annotations

# Distance to TSS
plotDistToTSS(peak_anno)  # Shows distance of peaks to nearest TSS

# TSS enrichment profile
tagMatrix <- getTagMatrix(    # Calculates peak density around TSS
  peaks,
  windows = getPromoters(     # Creates promoter windows around genes
    TxDb = txdb,
    upstream = 3000,
    downstream = 3000
  )
)

tagHeatmap(tagMatrix)        # Makes a heatmap of peaks around TSS

plotAvgProf(                 # Makes average TSS profile plot (Shows average peak enrichment around TSS)
  tagMatrix,
  xlim = c(-3000, 3000),
  xlab = "Genomic Region Around TSS",
  ylab = "Peak Frequency"
)

# Split ChIPseeker annotated peaks into promoter and non-promoter groups
peak_df <- as.data.frame(peak_anno)

promoter_peaks <- peak_df %>%
  filter(grepl("Promoter", annotation))

non_promoter_peaks <- peak_df %>%
  filter(!grepl("Promoter", annotation))


# Check numbers
nrow(promoter_peaks)
nrow(non_promoter_peaks)

# Check Ebox overlap analysis
ebox_dir <- "Ebox_predicted"

# E-box overlap analysis using  EBOX_all_reference.bed
# Load E-box reference BED file
ebox_df <- read.table(
  "EBOX_all_reference.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(ebox_df)[1:3] <- c("chr", "start", "end")

# Read peak BED file as dataframe
peaks_df <- read.table(
  peak_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(peaks_df)[1:3] <- c("chr", "start", "end")

# Create consistent peak IDs
peaks_df$peak_id <- paste0(
  peaks_df$chr,
  ":",
  peaks_df$start,
  "-",
  peaks_df$end
)

# Convert peaks and E-box reference to GRanges
peaks_gr <- GRanges(
  seqnames = peaks_df$chr,
  ranges = IRanges(
    start = peaks_df$start,
    end = peaks_df$end
  )
)

ebox_gr <- GRanges(
  seqnames = ebox_df$chr,
  ranges = IRanges(
    start = ebox_df$start,
    end = ebox_df$end
  )
)

# Find peaks overlapping at least one E-box motif
hits <- findOverlaps(peaks_gr, ebox_gr)

peak_hit_index <- unique(queryHits(hits))

peaks_with_ebox <- peaks_df[peak_hit_index, ]

# Save E-box-containing peaks
write.table(
  peaks_with_ebox,
  file = "TCF4_NPC_peaks_with_EBOX_ForrestCriteria.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Check number of E-box peaks
nrow(peaks_with_ebox)

# Annotate E-box-containing peaks
ebox_peak_file <- "TCF4_NPC_peaks_with_EBOX_ForrestCriteria.bed"

ebox_peaks <- readPeakFile(ebox_peak_file)

ebox_peak_anno <- annotatePeak(
  ebox_peaks,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

ebox_anno_df <- as.data.frame(ebox_peak_anno)

write_csv(
  ebox_anno_df,
  "TCF4_NPC_EBOX_peaks_ChIPseeker_annotation_ForrestCriteria.csv"
)

# ============================================================
# Promoter vs non-promoter E-box comparison
# ============================================================

# Convert full annotated peaks to GRanges
anno_gr <- GRanges(
  seqnames = peak_df$seqnames,
  ranges = IRanges(
    start = peak_df$start,
    end = peak_df$end
  )
)

# Convert E-box-containing peaks to GRanges
ebox_peak_gr <- GRanges(
  seqnames = peaks_with_ebox$chr,
  ranges = IRanges(
    start = peaks_with_ebox$start,
    end = peaks_with_ebox$end
  )
)

# Find overlaps
ebox_hits <- findOverlaps(anno_gr, ebox_peak_gr)

# Create E-box status column
peak_df$ebox_status <- "No_Ebox"

peak_df$ebox_status[
  unique(queryHits(ebox_hits))
] <- "Ebox"

# Mark promoter vs non-promoter
peak_df$region_type <- ifelse(
  grepl("Promoter", peak_df$annotation),
  "Promoter",
  "Non_promoter"
)

# Build 2x2 table
ebox_table <- table(
  peak_df$region_type,
  peak_df$ebox_status
)

ebox_table

# Percent within promoter/non-promoter groups
ebox_percent <- prop.table(
  ebox_table,
  margin = 1
) * 100

ebox_percent

# Fisher exact test
ebox_fisher <- fisher.test(ebox_table)

ebox_fisher

# ============================================================
# Extract gene sets
# ============================================================

# All TCF4 peak genes
gene_all <- unique(peak_df$SYMBOL)
gene_all <- gene_all[!is.na(gene_all)]

write.csv(
  data.frame(Gene = gene_all),
  "All_TCF4_peak_genes_ForrestCriteria.csv",
  row.names = FALSE
)

# E-box peak genes
gene_ebox <- unique(ebox_anno_df$SYMBOL)
gene_ebox <- gene_ebox[!is.na(gene_ebox)]

write.csv(
  data.frame(Gene = gene_ebox),
  "EBOX_TCF4_peak_genes_ForrestCriteria.csv",
  row.names = FALSE
)

# Promoter + E-box genes
promoter_ebox_genes <- unique(
  peak_df$SYMBOL[
    peak_df$ebox_status == "Ebox" &
      peak_df$region_type == "Promoter"
  ]
)

promoter_ebox_genes <- promoter_ebox_genes[!is.na(promoter_ebox_genes)]

write.csv(
  data.frame(Gene = promoter_ebox_genes),
  "Promoter_EBOX_TCF4_peak_genes_ForrestCriteria.csv",
  row.names = FALSE
)

# Non-promoter + E-box genes
nonpromoter_ebox_genes <- unique(
  peak_df$SYMBOL[
    peak_df$ebox_status == "Ebox" &
      peak_df$region_type == "Non_promoter"
  ]
)

nonpromoter_ebox_genes <- nonpromoter_ebox_genes[!is.na(nonpromoter_ebox_genes)]

write.csv(
  data.frame(Gene = nonpromoter_ebox_genes),
  "NonPromoter_EBOX_TCF4_peak_genes_ForrestCriteria.csv",
  row.names = FALSE
)

# ============================================================
# Final summary table
# ============================================================

summary_df <- data.frame(
  Metric = c(
    "Total TCF4 peaks",
    "Peaks with E-box",
    "Percent peaks with E-box",
    "All peak genes",
    "E-box peak genes",
    "Promoter + E-box genes",
    "Non-promoter + E-box genes",
    "Fisher odds ratio",
    "Fisher p-value"
  ),
  Value = c(
    nrow(peaks_df),
    nrow(peaks_with_ebox),
    round(100 * nrow(peaks_with_ebox) / nrow(peaks_df), 2),
    length(gene_all),
    length(gene_ebox),
    length(promoter_ebox_genes),
    length(nonpromoter_ebox_genes),
    round(ebox_fisher$estimate, 3),
    signif(ebox_fisher$p.value, 3)
  )
)

write_csv(
  summary_df,
  "TCF4_NPC_EBOX_analysis_summary_ForrestCriteria.csv"
)

summary_df
