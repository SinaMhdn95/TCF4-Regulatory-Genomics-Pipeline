# Peak-level and gene-level overlap analysis
#########################################################
# 1. McClay TCF4 consensus vs NPC TCF4 vs ENCODE NPC CTCF
# Genome build: hg38 / GRCh38
#########################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(  #Bioconductor packages
  "GenomeInfoDb",
  "GenomicRanges",
  "rtracklayer",
  "ChIPseeker",
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene"
))

BiocManager::install("GenomeInfoDb")
BiocManager::install("LOLA")

# We want to find the mechanism of TCF4 regulation in NPC and SH-SY5Y cells
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

# Import BED files
mcclay_peaks <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
npc_peaks <- import("NPC_hg38_summit_500bp.bed")
forrest_peak <- import("Forrest_hg38.bed")

# Read CTCF narrowPeak/BED file safely
# We can't import CTCF file easily because our CTCF file is not a simple 3-column BED file
# When we use import() rtracklayer tried to interpret the file as a standard BED/narrowPeak
# format and expected certain columns to be integers but Columns 7 to 9 contain decimal numbers
# So, we read it as a table and we only need three columns for overlap analysis chr, start and end.
ctcf_df <- read.table(
  "GSE123202_ENCFF102XIH_conservative_idr_thresholded_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

# Convert CTCF BED coordinates to GRanges
# Standardize the ranges into standard genomic interval objects
# BED start is 0-based, but GRanges is 1-based, so add +1 to start.
ctcf_peaks <- GRanges(
  seqnames = ctcf_df[, 1],
  ranges = IRanges(
    start = ctcf_df[, 2] + 1,
    end = ctcf_df[, 3]
  )
)

length(ctcf_peaks)

# Check chromosome naming
# All should be UCSC style: chr1, chr2, chrX, etc.
seqlevelsStyle(mcclay_peaks)
seqlevelsStyle(npc_peaks)
seqlevelsStyle(ctcf_peaks)
seqlevelsStyle(forrest_peak)

# Keep only standard chromosomes
# This removes random, unplaced, and alternative contigs
standard_chr <- paste0("chr", c(1:22, "X", "Y"))

mcclay_standard <- intersect(seqlevels(mcclay_peaks), standard_chr)
npc_standard <- intersect(seqlevels(npc_peaks), standard_chr)
ctcf_standard <- intersect(seqlevels(ctcf_peaks), standard_chr)
forrest_standard <- intersect(seqlevels(forrest_peak), standard_chr)

mcclay_peaks <- keepSeqlevels(
  mcclay_peaks,
  mcclay_standard,
  pruning.mode = "coarse"
)

npc_peaks <- keepSeqlevels(
  npc_peaks,
  npc_standard,
  pruning.mode = "coarse"
)

ctcf_peaks <- keepSeqlevels(
  ctcf_peaks,
  ctcf_standard,
  pruning.mode = "coarse"
)

forrest_peak <- keepSeqlevels(
  forrest_peak,
  forrest_standard,
  pruning.mode = "coarse"
)

length(forrest_peak)

# Peak-level overlaps
mcclay_ctcf <- subsetByOverlaps(mcclay_peaks, ctcf_peaks)   # McClay TCF4 peaks overlapping CTCF
npc_ctcf <- subsetByOverlaps(npc_peaks, ctcf_peaks)         # NPC TCF4 peaks overlapping CTCF
mcclay_npc <- subsetByOverlaps(mcclay_peaks, npc_peaks)     # McClay TCF4 peaks overlapping NPC TCF4 peaks
forrest_ctcf <- subsetByOverlaps(forrest_peak, ctcf_peaks)  # Forrest TCF4 peaks overlapping CTCF
forrest_npc <- subsetByOverlaps(forrest_peak, npc_peaks)  # Forrest TCF4 peaks overlapping NPC TCF4
forrest_mcclay <- subsetByOverlaps(forrest_peak, mcclay_peaks) #McClay TCF4 peaks overlapping Forrest

# Triple overlap:Peaks shared by McClay TCF4, NPC TCF4, and CTCF
triple_overlap <- subsetByOverlaps(mcclay_npc, ctcf_peaks)
four_overlap <- subsetByOverlaps(forrest_peak, triple_overlap)

# Count total peaks after cleaning
length(mcclay_peaks)
length(npc_peaks)
length(ctcf_peaks)
length(mcclay_ctcf)
length(npc_ctcf)
length(mcclay_npc)
length(triple_overlap)
length(forrest_ctcf)
length(forrest_npc)
length(forrest_mcclay)

# Summary table of peak overlaps
peak_overlap_summary <- data.frame(

  Comparison = c(
    "McClay_TCF4_vs_CTCF",
    "NPC_TCF4_vs_CTCF",
    "McClay_TCF4_vs_NPC_TCF4",
    "McClay_TCF4_vs_NPC_TCF4_vs_CTCF",
    "Forrest_vs_CTCF",
    "Forrest_vs_McClay_TCF4",
    "Forrest_vs_NPC_TCF4"
  ),

  Peak_Count = c(
    length(mcclay_ctcf),
    length(npc_ctcf),
    length(mcclay_npc),
    length(triple_overlap),
    length(forrest_ctcf),
    length(forrest_mcclay),
    length(forrest_npc)
  ),

  Percent_of_First_Dataset = c(
    length(mcclay_ctcf) / length(mcclay_peaks) * 100,
    length(npc_ctcf) / length(npc_peaks) * 100,
    length(mcclay_npc) / length(mcclay_peaks) * 100,
    length(triple_overlap) / length(mcclay_peaks) * 100,
    length(forrest_ctcf) / length(forrest_peak) * 100,
    length(forrest_mcclay) / length(mcclay_peaks) * 100,
    length(forrest_npc) / length(npc_peaks) * 100
  )
)

peak_overlap_summary

# Annotate original peak sets
# This assigns peaks to promoters, exons, introns, intergenic regions, and nearest genes
mcclay_anno <- annotatePeak(
  mcclay_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_anno <- annotatePeak(
  npc_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

ctcf_anno <- annotatePeak(
  ctcf_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_anno <- annotatePeak(
  forrest_peak,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)
# Annotate overlap peak sets
# These are more biologically specific than annotating all peaks
mcclay_ctcf_anno <- annotatePeak(
  mcclay_ctcf,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_ctcf_anno <- annotatePeak(
  npc_ctcf,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

triple_anno <- annotatePeak(
  triple_overlap,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

# Convert annotations to data frames
mcclay_anno_df <- as.data.frame(mcclay_anno)
npc_anno_df <- as.data.frame(npc_anno)
ctcf_anno_df <- as.data.frame(ctcf_anno)
mcclay_ctcf_anno_df <- as.data.frame(mcclay_ctcf_anno)
npc_ctcf_anno_df <- as.data.frame(npc_ctcf_anno)
triple_anno_df <- as.data.frame(triple_anno)
forrest_anno_df <- as.data.frame(forrest_anno)

# Extract gene symbols
mcclay_genes <- unique(na.omit(mcclay_anno_df$SYMBOL))
npc_genes <- unique(na.omit(npc_anno_df$SYMBOL))
ctcf_genes <- unique(na.omit(ctcf_anno_df$SYMBOL))
mcclay_ctcf_genes <- unique(na.omit(mcclay_ctcf_anno_df$SYMBOL))
npc_ctcf_genes <- unique(na.omit(npc_ctcf_anno_df$SYMBOL))
triple_genes <- unique(na.omit(triple_anno_df$SYMBOL))
forrest_genes <- unique(na.omit(forrest_anno_df$SYMBOL))

# Gene-level overlap summary
gene_overlap_summary <- data.frame(

  Comparison = c(
    "McClay_genes_vs_CTCF_genes",
    "NPC_genes_vs_CTCF_genes",
    "McClay_genes_vs_NPC_genes",
    "McClay_genes_vs_NPC_genes_vs_CTCF_genes",
    "Forrest_genes_vs_CTCF_genes",
    "Forrest_genes_vs_NPC_genes",
    "Forrest_genes_vs_McClay_genes"
  ),

  Overlap_Genes = c(
    length(intersect(mcclay_genes, ctcf_genes)),
    length(intersect(npc_genes, ctcf_genes)),
    length(intersect(mcclay_genes, npc_genes)),
    length(Reduce(intersect,
                  list(mcclay_genes,
                       npc_genes,
                       ctcf_genes))),
    length(intersect(forrest_genes, ctcf_genes)),
    length(intersect(forrest_genes, npc_genes)),
    length(intersect(forrest_genes, mcclay_genes))
  ),

  First_Gene_Set_Size = c(
    length(mcclay_genes),
    length(npc_genes),
    length(mcclay_genes),
    length(mcclay_genes),
    length(forrest_genes),
    length(forrest_genes),
    length(forrest_genes)
  ),

  Percent_of_First_Gene_Set = round(c(

    length(intersect(mcclay_genes, ctcf_genes)) /
      length(mcclay_genes) * 100,

    length(intersect(npc_genes, ctcf_genes)) /
      length(npc_genes) * 100,

    length(intersect(mcclay_genes, npc_genes)) /
      length(mcclay_genes) * 100,

    length(Reduce(intersect,
                  list(mcclay_genes,
                       npc_genes,
                       ctcf_genes))) /
      length(mcclay_genes) * 100,

    length(intersect(forrest_genes, ctcf_genes)) /
      length(forrest_genes) * 100,

    length(intersect(forrest_genes, npc_genes)) /
      length(forrest_genes) * 100,

    length(intersect(forrest_genes, mcclay_genes)) /
      length(forrest_genes) * 100

  ), 2),

  Second_Gene_Set_Size = c(
    length(ctcf_genes),
    length(ctcf_genes),
    length(npc_genes),
    length(ctcf_genes),
    length(ctcf_genes),
    length(npc_genes),
    length(mcclay_genes)
  ),

  Percent_of_Second_Gene_Set = round(c(

    length(intersect(mcclay_genes, ctcf_genes)) /
      length(ctcf_genes) * 100,

    length(intersect(npc_genes, ctcf_genes)) /
      length(ctcf_genes) * 100,

    length(intersect(mcclay_genes, npc_genes)) /
      length(npc_genes) * 100,

    length(Reduce(intersect,
                  list(mcclay_genes,
                       npc_genes,
                       ctcf_genes))) /
      length(ctcf_genes) * 100,

    length(intersect(forrest_genes, ctcf_genes)) /
      length(ctcf_genes) * 100,

    length(intersect(forrest_genes, npc_genes)) /
      length(npc_genes) * 100,

    length(intersect(forrest_genes, mcclay_genes)) /
      length(mcclay_genes) * 100

  ), 2)

)

gene_overlap_summary

write.csv(
  gene_overlap_summary,
  "gene_overlap_summary_with_percentages.csv",
  row.names = FALSE
)


length(intersect(forrest_genes, mcclay_genes))

############################################################
# 2. FGFR1 and NOTCH peaks in SCZ vs NPC TCF4 vs Dr.MCclay's
# Genome build: hg19 / GRCh19
############################################################

# Import data
control_fgfr1 <- import("FGFR1_NOTCH/GSM2439179_C4-R1.bed.gz")
patient_fgfr1 <- import("FGFR1_NOTCH/GSM2439181_P3-R1.bed.gz")
patient_notch <- import("FGFR1_NOTCH/GSM2439182_P3-N.bed.gz")

length(control_fgfr1)
length(patient_fgfr1)
length(patient_notch)

# Check chromosome style
seqlevels(control_fgfr1)[1:10]
seqlevels(patient_fgfr1)[1:10]
seqlevels(patient_notch)[1:10]

# Download chain file to do liftOver for fgfr1 peaks
download.file(
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
  destfile = "hg19ToHg38.over.chain.gz",
  mode = "wb"
)

R.utils::gunzip("hg19ToHg38.over.chain.gz", remove = FALSE)  #unzip the chain file (a translation map between two genome assemblies)

# Import chain
chain <- import.chain("hg19ToHg38.over.chain")

# LiftOver
control_fgfr1_hg38_list <- liftOver(control_fgfr1, chain)
patient_fgfr1_hg38_list <- liftOver(patient_fgfr1, chain)
patient_notch_hg38_list <- liftOver(patient_notch, chain)

length(control_fgfr1_hg38_list)
length(patient_fgfr1_hg38_list)
length(patient_notch_hg38_list)

# Make one-to-one clean version
keep_one_to_one <- function(grl) {
unlist(grl[elementNROWS(grl) == 1])
}

control_fgfr1_hg38_clean <- keep_one_to_one(control_fgfr1_hg38_list)
patient_fgfr1_hg38_clean <- keep_one_to_one(patient_fgfr1_hg38_list)
patient_notch_hg38_clean <- keep_one_to_one(patient_notch_hg38_list)

length(control_fgfr1_hg38_clean)
length(patient_fgfr1_hg38_clean)
length(patient_notch_hg38_clean)

# Export hg38 FGFR1 and NOTCH peaks
export(control_fgfr1_hg38_clean,
       "FGFR1_NOTCH/Control_FGFR1_hg38_clean.bed")

export(patient_fgfr1_hg38_clean,
       "FGFR1_NOTCH/Patient_FGFR1_hg38_clean.bed")

export(patient_notch_hg38_clean,
       "FGFR1_NOTCH/Patient_NOTCH_hg38_clean.bed")


length(mcclay_peaks)
length(npc_peaks)
length(forrest_peak)

length(control_fgfr1_hg38_clean)
length(patient_fgfr1_hg38_clean)
length(patient_notch_hg38_clean)

# Check the peak overlaps with NPC, Forrst and McClay's peaks
overlap_summary <- function(query, subject, query_name, subject_name) {
  ov <- findOverlaps(query, subject)

  data.frame(
    Comparison = paste(query_name, "vs", subject_name),
    Query_Peaks = length(query),
    Subject_Peaks = length(subject),
    Overlap_Peaks = length(unique(queryHits(ov))),
    Percent_Query_Overlap = length(unique(queryHits(ov))) / length(query) * 100,
    Subject_Overlapped = length(unique(subjectHits(ov))),
    Percent_Subject_Overlap = length(unique(subjectHits(ov))) / length(subject) * 100
  )
}

results <- rbind(
  overlap_summary(mcclay_peaks, control_fgfr1_hg38_clean, "McClay_TCF4", "Control_FGFR1"),
  overlap_summary(mcclay_peaks, patient_fgfr1_hg38_clean, "McClay_TCF4", "Patient_FGFR1"),
  overlap_summary(mcclay_peaks, patient_notch_hg38_clean, "McClay_TCF4", "Patient_NOTCH"),

  overlap_summary(npc_peaks, control_fgfr1_hg38_clean, "NPC_TCF4", "Control_FGFR1"),
  overlap_summary(npc_peaks, patient_fgfr1_hg38_clean, "NPC_TCF4", "Patient_FGFR1"),
  overlap_summary(npc_peaks, patient_notch_hg38_clean, "NPC_TCF4", "Patient_NOTCH"),

  overlap_summary(forrest_peak, control_fgfr1_hg38_clean, "Forrest_TCF4", "Control_FGFR1"),
  overlap_summary(forrest_peak, patient_fgfr1_hg38_clean, "Forrest_TCF4", "Patient_FGFR1"),
  overlap_summary(forrest_peak, patient_notch_hg38_clean, "Forrest_TCF4", "Patient_NOTCH")
)

results


# Double check manually (outside function made by AI too see if the code AI gave us worked right?!)
# Find overlaps
hits <- findOverlaps(
  npc_peaks,
  control_fgfr1_hg38_clean
)

# Number of unique NPC peaks that overlap FGFR1
overlap_peaks <- length(
  unique(queryHits(hits))
)

# Number of unique FGFR1 peaks overlapped by NPC peaks
subject_overlapped <- length(
  unique(subjectHits(hits))
)

# Percentage of NPC peaks overlapping FGFR1
percent_query_overlap <-
  overlap_peaks / length(npc_peaks) * 100

# Percentage of FGFR1 peaks overlapped by NPC peaks
percent_subject_overlap <-
  subject_overlapped / length(control_fgfr1_hg38_clean) * 100

# Print results
cat("NPC peaks:", length(npc_peaks), "\n")
cat("FGFR1 peaks:", length(control_fgfr1_hg38_clean), "\n")
cat("Overlapping NPC peaks:", overlap_peaks, "\n")
cat("Percent NPC overlap:", percent_query_overlap, "\n")
cat("Overlapping FGFR1 peaks:", subject_overlapped, "\n")
cat("Percent FGFR1 overlap:", percent_subject_overlap, "\n")

results[results$Comparison == "NPC_TCF4 vs Control_FGFR1", ]


# Number of overlap pairs
length(hits)

# Number of unique query peaks
length(unique(queryHits(hits)))

# Number of unique subject peaks
length(unique(subjectHits(hits)))


# Gene-level overlap: TCF4 vs FGFR1 / NOTCH
# Annotate FGFR1 and NOTCH hg38 peak sets
control_fgfr1_anno <- annotatePeak(
  control_fgfr1_hg38_clean,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

patient_fgfr1_anno <- annotatePeak(
  patient_fgfr1_hg38_clean,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

patient_notch_anno <- annotatePeak(
  patient_notch_hg38_clean,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

# Convert to data frames
control_fgfr1_anno_df <- as.data.frame(control_fgfr1_anno)
patient_fgfr1_anno_df <- as.data.frame(patient_fgfr1_anno)
patient_notch_anno_df <- as.data.frame(patient_notch_anno)

# Extract unique gene symbols
control_fgfr1_genes <- unique(na.omit(control_fgfr1_anno_df$SYMBOL))
patient_fgfr1_genes <- unique(na.omit(patient_fgfr1_anno_df$SYMBOL))
patient_notch_genes <- unique(na.omit(patient_notch_anno_df$SYMBOL))

# Check gene set sizes
length(mcclay_genes)
length(npc_genes)
length(forrest_genes)

length(control_fgfr1_genes)
length(patient_fgfr1_genes)
length(patient_notch_genes)

# Gene overlap summary table
gene_overlap_fgfr1_notch <- data.frame(

  Comparison = c(
    "McClay_genes_vs_Control_FGFR1_genes",
    "McClay_genes_vs_Patient_FGFR1_genes",
    "McClay_genes_vs_Patient_NOTCH_genes",

    "NPC_genes_vs_Control_FGFR1_genes",
    "NPC_genes_vs_Patient_FGFR1_genes",
    "NPC_genes_vs_Patient_NOTCH_genes",

    "Forrest_genes_vs_Control_FGFR1_genes",
    "Forrest_genes_vs_Patient_FGFR1_genes",
    "Forrest_genes_vs_Patient_NOTCH_genes"
  ),

  Overlap_Genes = c(
    length(intersect(mcclay_genes, control_fgfr1_genes)),
    length(intersect(mcclay_genes, patient_fgfr1_genes)),
    length(intersect(mcclay_genes, patient_notch_genes)),

    length(intersect(npc_genes, control_fgfr1_genes)),
    length(intersect(npc_genes, patient_fgfr1_genes)),
    length(intersect(npc_genes, patient_notch_genes)),

    length(intersect(forrest_genes, control_fgfr1_genes)),
    length(intersect(forrest_genes, patient_fgfr1_genes)),
    length(intersect(forrest_genes, patient_notch_genes))
  ),

  First_Gene_Set_Size = c(
    length(mcclay_genes),
    length(mcclay_genes),
    length(mcclay_genes),

    length(npc_genes),
    length(npc_genes),
    length(npc_genes),

    length(forrest_genes),
    length(forrest_genes),
    length(forrest_genes)
  ),

  Second_Gene_Set_Size = c(
    length(control_fgfr1_genes),
    length(patient_fgfr1_genes),
    length(patient_notch_genes),

    length(control_fgfr1_genes),
    length(patient_fgfr1_genes),
    length(patient_notch_genes),

    length(control_fgfr1_genes),
    length(patient_fgfr1_genes),
    length(patient_notch_genes)
  )
)

# Add percentages
gene_overlap_fgfr1_notch$Percent_of_First_Gene_Set <-
  round(
    gene_overlap_fgfr1_notch$Overlap_Genes /
      gene_overlap_fgfr1_notch$First_Gene_Set_Size * 100,
    2
  )

gene_overlap_fgfr1_notch$Percent_of_Second_Gene_Set <-
  round(
    gene_overlap_fgfr1_notch$Overlap_Genes /
      gene_overlap_fgfr1_notch$Second_Gene_Set_Size * 100,
    2
  )

gene_overlap_fgfr1_notch


# Fisher exact test for gene-level overlaps
# Background universe: final_ref_gene_2026_filtered.csv

# Import background universe
background_gene <- read.csv(
  "final_ref_gene_2026_filtered.csv",
  stringsAsFactors = FALSE)

# Check column names first
colnames(background_gene)

# Use gene column as background universe
allgenes <- unique(na.omit(background_gene$gene))

# Clean background gene symbols
allgenes <- unique(toupper(allgenes))

length(allgenes)
head(allgenes)

# Clean all gene sets
mcclay_genes_clean <- intersect(unique(toupper(mcclay_genes)), allgenes)
npc_genes_clean <- intersect(unique(toupper(npc_genes)), allgenes)
forrest_genes_clean <- intersect(unique(toupper(forrest_genes)), allgenes)

control_fgfr1_genes_clean <- intersect(unique(toupper(control_fgfr1_genes)), allgenes)
patient_fgfr1_genes_clean <- intersect(unique(toupper(patient_fgfr1_genes)), allgenes)
patient_notch_genes_clean <- intersect(unique(toupper(patient_notch_genes)), allgenes)

# Fisher-test function
fisher_gene_overlap <- function(set1, set2, universe, set1_name, set2_name) {

  set1 <- intersect(unique(set1), universe)
  set2 <- intersect(unique(set2), universe)

  a <- length(intersect(set1, set2))
  b <- length(setdiff(set1, set2))
  c <- length(setdiff(set2, set1))
  d <- length(setdiff(universe, union(set1, set2)))

  fisher_matrix <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE
  )

  rownames(fisher_matrix) <- c(set1_name, paste0("Not_", set1_name))
  colnames(fisher_matrix) <- c(set2_name, paste0("Not_", set2_name))

  fisher_result <- fisher.test(
    fisher_matrix,
    alternative = "greater"
  )

  data.frame(
    Comparison = paste(set1_name, "vs", set2_name),
    Set1_Size = length(set1),
    Set2_Size = length(set2),
    Overlap = a,
    Set1_Only = b,
    Set2_Only = c,
    Neither = d,
    Odds_Ratio = unname(fisher_result$estimate),
    P_value = fisher_result$p.value
  )
}

fisher_results_fgfr1_notch <- rbind(

  fisher_gene_overlap(
    mcclay_genes_clean,
    control_fgfr1_genes_clean,
    allgenes,
    "McClay_TCF4",
    "Control_FGFR1"
  ),

  fisher_gene_overlap(
    mcclay_genes_clean,
    patient_fgfr1_genes_clean,
    allgenes,
    "McClay_TCF4",
    "Patient_FGFR1"
  ),

  fisher_gene_overlap(
    mcclay_genes_clean,
    patient_notch_genes_clean,
    allgenes,
    "McClay_TCF4",
    "Patient_NOTCH"
  ),

  fisher_gene_overlap(
    npc_genes_clean,
    control_fgfr1_genes_clean,
    allgenes,
    "NPC_TCF4",
    "Control_FGFR1"
  ),

  fisher_gene_overlap(
    npc_genes_clean,
    patient_fgfr1_genes_clean,
    allgenes,
    "NPC_TCF4",
    "Patient_FGFR1"
  ),

  fisher_gene_overlap(
    npc_genes_clean,
    patient_notch_genes_clean,
    allgenes,
    "NPC_TCF4",
    "Patient_NOTCH"
  ),

  fisher_gene_overlap(
    forrest_genes_clean,
    control_fgfr1_genes_clean,
    allgenes,
    "Forrest_TCF4",
    "Control_FGFR1"
  ),

  fisher_gene_overlap(
    forrest_genes_clean,
    patient_fgfr1_genes_clean,
    allgenes,
    "Forrest_TCF4",
    "Patient_FGFR1"
  ),

  fisher_gene_overlap(
    forrest_genes_clean,
    patient_notch_genes_clean,
    allgenes,
    "Forrest_TCF4",
    "Patient_NOTCH"
  )
)

fisher_results_fgfr1_notch

############################################################
# NPC TCF4 vs CTCF vs FGFR1 vs NOTCH
# Peak-level and gene-level overlap
############################################################

# 1. Peak-level overlap
# Pairwise overlaps with NPC TCF4
npc_ctcf_peaks <- subsetByOverlaps(
  npc_peaks,
  ctcf_peaks
)

npc_control_fgfr1_peaks <- subsetByOverlaps(
  npc_peaks,
  control_fgfr1_hg38_clean
)

npc_patient_fgfr1_peaks <- subsetByOverlaps(
  npc_peaks,
  patient_fgfr1_hg38_clean
)

npc_patient_notch_peaks <- subsetByOverlaps(
  npc_peaks,
  patient_notch_hg38_clean
)

# Triple overlaps: NPC TCF4 + CTCF + each factor
npc_ctcf_control_fgfr1_peaks <- subsetByOverlaps(
  npc_ctcf_peaks,
  control_fgfr1_hg38_clean
)

npc_ctcf_patient_fgfr1_peaks <- subsetByOverlaps(
  npc_ctcf_peaks,
  patient_fgfr1_hg38_clean
)

npc_ctcf_patient_notch_peaks <- subsetByOverlaps(
  npc_ctcf_peaks,
  patient_notch_hg38_clean
)

# Four-way overlap: NPC TCF4 + CTCF + Patient FGFR1 + Patient NOTCH
npc_ctcf_patient_fgfr1_notch_peaks <- subsetByOverlaps(
  npc_ctcf_patient_fgfr1_peaks,
  patient_notch_hg38_clean
)

# Peak summary table
npc_ctcf_fgfr1_notch_peak_summary <- data.frame(
  Comparison = c(
    "NPC_TCF4_vs_CTCF",
    "NPC_TCF4_vs_Control_FGFR1",
    "NPC_TCF4_vs_Patient_FGFR1",
    "NPC_TCF4_vs_Patient_NOTCH",
    "NPC_TCF4_vs_CTCF_vs_Control_FGFR1",
    "NPC_TCF4_vs_CTCF_vs_Patient_FGFR1",
    "NPC_TCF4_vs_CTCF_vs_Patient_NOTCH",
    "NPC_TCF4_vs_CTCF_vs_Patient_FGFR1_vs_Patient_NOTCH"
  ),

  Peak_Count = c(
    length(npc_ctcf_peaks),
    length(npc_control_fgfr1_peaks),
    length(npc_patient_fgfr1_peaks),
    length(npc_patient_notch_peaks),
    length(npc_ctcf_control_fgfr1_peaks),
    length(npc_ctcf_patient_fgfr1_peaks),
    length(npc_ctcf_patient_notch_peaks),
    length(npc_ctcf_patient_fgfr1_notch_peaks)
  ),

  Percent_of_NPC_TCF4 = round(c(
    length(npc_ctcf_peaks) / length(npc_peaks) * 100,
    length(npc_control_fgfr1_peaks) / length(npc_peaks) * 100,
    length(npc_patient_fgfr1_peaks) / length(npc_peaks) * 100,
    length(npc_patient_notch_peaks) / length(npc_peaks) * 100,
    length(npc_ctcf_control_fgfr1_peaks) / length(npc_peaks) * 100,
    length(npc_ctcf_patient_fgfr1_peaks) / length(npc_peaks) * 100,
    length(npc_ctcf_patient_notch_peaks) / length(npc_peaks) * 100,
    length(npc_ctcf_patient_fgfr1_notch_peaks) / length(npc_peaks) * 100
  ), 2),

  Percent_of_NPC_CTCF = round(c(
    NA,
    NA,
    NA,
    NA,
    length(npc_ctcf_control_fgfr1_peaks) / length(npc_ctcf_peaks) * 100,
    length(npc_ctcf_patient_fgfr1_peaks) / length(npc_ctcf_peaks) * 100,
    length(npc_ctcf_patient_notch_peaks) / length(npc_ctcf_peaks) * 100,
    length(npc_ctcf_patient_fgfr1_notch_peaks) / length(npc_ctcf_peaks) * 100
  ), 2)
)

npc_ctcf_fgfr1_notch_peak_summary


# 2. Gene-level overlap
# Pairwise shared genes
npc_ctcf_shared_genes <- intersect(
  npc_genes,
  ctcf_genes
)

npc_control_fgfr1_shared_genes <- intersect(
  npc_genes,
  control_fgfr1_genes
)

npc_patient_fgfr1_shared_genes <- intersect(
  npc_genes,
  patient_fgfr1_genes
)

npc_patient_notch_shared_genes <- intersect(
  npc_genes,
  patient_notch_genes
)

# Triple shared genes
npc_ctcf_control_fgfr1_shared_genes <- Reduce(
  intersect,
  list(npc_genes, ctcf_genes, control_fgfr1_genes)
)

npc_ctcf_patient_fgfr1_shared_genes <- Reduce(
  intersect,
  list(npc_genes, ctcf_genes, patient_fgfr1_genes)
)

npc_ctcf_patient_notch_shared_genes <- Reduce(
  intersect,
  list(npc_genes, ctcf_genes, patient_notch_genes)
)

# Four-way shared genes
npc_ctcf_patient_fgfr1_notch_shared_genes <- Reduce(
  intersect,
  list(npc_genes, ctcf_genes, patient_fgfr1_genes, patient_notch_genes)
)

# Gene summary table
npc_ctcf_fgfr1_notch_gene_summary <- data.frame(
  Comparison = c(
    "NPC_genes_vs_CTCF_genes",
    "NPC_genes_vs_Control_FGFR1_genes",
    "NPC_genes_vs_Patient_FGFR1_genes",
    "NPC_genes_vs_Patient_NOTCH_genes",
    "NPC_genes_vs_CTCF_vs_Control_FGFR1_genes",
    "NPC_genes_vs_CTCF_vs_Patient_FGFR1_genes",
    "NPC_genes_vs_CTCF_vs_Patient_NOTCH_genes",
    "NPC_genes_vs_CTCF_vs_Patient_FGFR1_vs_Patient_NOTCH_genes"
  ),

  Gene_Count = c(
    length(npc_ctcf_shared_genes),
    length(npc_control_fgfr1_shared_genes),
    length(npc_patient_fgfr1_shared_genes),
    length(npc_patient_notch_shared_genes),
    length(npc_ctcf_control_fgfr1_shared_genes),
    length(npc_ctcf_patient_fgfr1_shared_genes),
    length(npc_ctcf_patient_notch_shared_genes),
    length(npc_ctcf_patient_fgfr1_notch_shared_genes)
  ),

  Percent_of_NPC_Genes = round(c(
    length(npc_ctcf_shared_genes) / length(npc_genes) * 100,
    length(npc_control_fgfr1_shared_genes) / length(npc_genes) * 100,
    length(npc_patient_fgfr1_shared_genes) / length(npc_genes) * 100,
    length(npc_patient_notch_shared_genes) / length(npc_genes) * 100,
    length(npc_ctcf_control_fgfr1_shared_genes) / length(npc_genes) * 100,
    length(npc_ctcf_patient_fgfr1_shared_genes) / length(npc_genes) * 100,
    length(npc_ctcf_patient_notch_shared_genes) / length(npc_genes) * 100,
    length(npc_ctcf_patient_fgfr1_notch_shared_genes) / length(npc_genes) * 100
  ), 2),

  Percent_of_NPC_CTCF_Genes = round(c(
    NA,
    NA,
    NA,
    NA,
    length(npc_ctcf_control_fgfr1_shared_genes) / length(npc_ctcf_shared_genes) * 100,
    length(npc_ctcf_patient_fgfr1_shared_genes) / length(npc_ctcf_shared_genes) * 100,
    length(npc_ctcf_patient_notch_shared_genes) / length(npc_ctcf_shared_genes) * 100,
    length(npc_ctcf_patient_fgfr1_notch_shared_genes) / length(npc_ctcf_shared_genes) * 100
  ), 2)
)

npc_ctcf_fgfr1_notch_gene_summary


############################################################
# 2. H3K27ac peaks vs NPC TCF4 vs Dr.MCclay's
# Genome build: hg38 / GRCh38
############################################################
# Import BED files
mcclay_peaks <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
npc_peaks <- import("NPC_hg38_summit_500bp.bed")
forrest_peak <- import("Forrest_hg38.bed")


# Since  file has a decimal value in the BED score column, but
# rtracklayer::import() expects that column to be an integer
# We import it in another way
H3K27ac_df <- read.table(
  "GSE96178_ENCFF407DWP_replicated_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

# Keep only chr, start, end
H3K27ac_df <- H3K27ac_df[, 1:3]
colnames(H3K27ac_df) <- c("chr", "start", "end")

H3K27ac_df$start <- as.numeric(H3K27ac_df$start)
H3K27ac_df$end <- as.numeric(H3K27ac_df$end)


# BED is 0-based; GRanges is 1-based
H3K27ac_peak <- GRanges(
  seqnames = H3K27ac_df$chr,
  ranges = IRanges(
    start = H3K27ac_df$start + 1,
    end = H3K27ac_df$end
  )
)

seqlevelsStyle(H3K27ac_peak) <- "UCSC"
H3K27ac_peak <- keepStandardChromosomes(H3K27ac_peak, pruning.mode = "coarse")
H3K27ac_peak <- sort(H3K27ac_peak)

length(H3K27ac_peak)

# Set universe background gene
background_df <- read.csv("final_ref_gene_2026_filtered.csv")

background_genes <- unique(
  background_df$gene
)

background_genes <- background_genes[
  !is.na(background_genes) &
    background_genes != ""
]

length(background_genes)

# Clean all peak files
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

seqlevelsStyle(mcclay_peaks) <- "UCSC"
mcclay_peaks <- keepStandardChromosomes(mcclay_peaks, pruning.mode = "coarse")
mcclay_peaks <- sort(mcclay_peaks)

seqlevelsStyle(npc_peaks) <- "UCSC"
npc_peaks <- keepStandardChromosomes(npc_peaks, pruning.mode = "coarse")
npc_peaks <- sort(npc_peaks)

seqlevelsStyle(forrest_peak) <- "UCSC"
forrest_peak <- keepStandardChromosomes(forrest_peak, pruning.mode = "coarse")
forrest_peak <- sort(forrest_peak)

seqlevelsStyle(H3K27ac_peak) <- "UCSC"
H3K27ac_peak <- keepStandardChromosomes(H3K27ac_peak, pruning.mode = "coarse")
H3K27ac_peak <- sort(H3K27ac_peak)

length(mcclay_peaks)
length(npc_peaks)
length(forrest_peak)
length(H3K27ac_peak)


# Peak-level overlap: TCF4 datasets vs H3K27ac
# NPC TCF4 vs H3K27ac
npc_h3k27ac_hits <- findOverlaps(npc_peaks, H3K27ac_peak)

npc_h3k27ac_query_peaks <- unique(queryHits(npc_h3k27ac_hits))
npc_h3k27ac_subject_peaks <- unique(subjectHits(npc_h3k27ac_hits))

length(npc_h3k27ac_query_peaks)
length(npc_h3k27ac_subject_peaks)

# McClay TCF4 vs H3K27ac
mcclay_h3k27ac_hits <- findOverlaps(mcclay_peaks, H3K27ac_peak)

mcclay_h3k27ac_query_peaks <- unique(queryHits(mcclay_h3k27ac_hits))
mcclay_h3k27ac_subject_peaks <- unique(subjectHits(mcclay_h3k27ac_hits))

length(mcclay_h3k27ac_query_peaks)
length(mcclay_h3k27ac_subject_peaks)

# Forrest TCF4 vs H3K27ac
forrest_h3k27ac_hits <- findOverlaps(forrest_peak, H3K27ac_peak)

forrest_h3k27ac_query_peaks <- unique(queryHits(forrest_h3k27ac_hits))
forrest_h3k27ac_subject_peaks <- unique(subjectHits(forrest_h3k27ac_hits))

length(forrest_h3k27ac_query_peaks)
length(forrest_h3k27ac_subject_peaks)

# Peak-level overlap summary table
h3k27ac_peak_overlap_summary <- data.frame(
  Comparison = c(
    "NPC_TCF4_vs_H3K27ac",
    "McClay_TCF4_vs_H3K27ac",
    "Forrest_TCF4_vs_H3K27ac"
  ),
  TCF4_Total_Peaks = c(
    length(npc_peaks),
    length(mcclay_peaks),
    length(forrest_peak)
  ),
  H3K27ac_Total_Peaks = c(
    length(H3K27ac_peak),
    length(H3K27ac_peak),
    length(H3K27ac_peak)
  ),
  TCF4_Overlapping_Peaks = c(
    length(npc_h3k27ac_query_peaks),
    length(mcclay_h3k27ac_query_peaks),
    length(forrest_h3k27ac_query_peaks)
  ),
  H3K27ac_Overlapping_Peaks = c(
    length(npc_h3k27ac_subject_peaks),
    length(mcclay_h3k27ac_subject_peaks),
    length(forrest_h3k27ac_subject_peaks)
  )
)

h3k27ac_peak_overlap_summary$Percent_TCF4_Peaks_Overlapping_H3K27ac <-
  round(
    h3k27ac_peak_overlap_summary$TCF4_Overlapping_Peaks /
      h3k27ac_peak_overlap_summary$TCF4_Total_Peaks * 100,
    3
  )

h3k27ac_peak_overlap_summary$Percent_H3K27ac_Peaks_Overlapping_TCF4 <-
  round(
    h3k27ac_peak_overlap_summary$H3K27ac_Overlapping_Peaks /
      h3k27ac_peak_overlap_summary$H3K27ac_Total_Peaks * 100,
    3
  )

h3k27ac_peak_overlap_summary


# Create H3K27ac-overlapping TCF4 peak objects
npc_TCF4_H3K27ac_peaks <- npc_peaks[npc_h3k27ac_query_peaks]

mcclay_TCF4_H3K27ac_peaks <- mcclay_peaks[mcclay_h3k27ac_query_peaks]

forrest_TCF4_H3K27ac_peaks <- forrest_peak[forrest_h3k27ac_query_peaks]

length(npc_TCF4_H3K27ac_peaks)
length(mcclay_TCF4_H3K27ac_peaks)
length(forrest_TCF4_H3K27ac_peaks)

############################################################
# Gene-level analysis: H3K27ac-overlapping TCF4 genes
# Background: final_ref_gene_2026_filtered.csv
############################################################
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# Import background gene universe
background_df <- read.csv("final_ref_gene_2026_filtered.csv")

background_genes <- unique(background_df$gene)

background_genes <- background_genes[
  !is.na(background_genes) &
    background_genes != ""
]

length(background_genes)


# Annotate all TCF4 peak sets to genes
npc_all_anno <- annotatePeak(
  npc_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_all_anno_df <- as.data.frame(npc_all_anno)

npc_all_genes <- unique(npc_all_anno_df$SYMBOL)
npc_all_genes <- npc_all_genes[!is.na(npc_all_genes) & npc_all_genes != ""]

mcclay_all_anno <- annotatePeak(
  mcclay_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_all_anno_df <- as.data.frame(mcclay_all_anno)

mcclay_all_genes <- unique(mcclay_all_anno_df$SYMBOL)
mcclay_all_genes <- mcclay_all_genes[!is.na(mcclay_all_genes) & mcclay_all_genes != ""]

forrest_all_anno <- annotatePeak(
  forrest_peak,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_all_anno_df <- as.data.frame(forrest_all_anno)

forrest_all_genes <- unique(forrest_all_anno_df$SYMBOL)
forrest_all_genes <- forrest_all_genes[!is.na(forrest_all_genes) & forrest_all_genes != ""]

# Annotate H3K27ac-overlapping TCF4 peaks to genes
npc_h3k27ac_anno <- annotatePeak(
  npc_TCF4_H3K27ac_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_h3k27ac_anno_df <- as.data.frame(npc_h3k27ac_anno)

npc_h3k27ac_genes <- unique(npc_h3k27ac_anno_df$SYMBOL)
npc_h3k27ac_genes <- npc_h3k27ac_genes[
  !is.na(npc_h3k27ac_genes) &
    npc_h3k27ac_genes != ""
]

mcclay_h3k27ac_anno <- annotatePeak(
  mcclay_TCF4_H3K27ac_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_h3k27ac_anno_df <- as.data.frame(mcclay_h3k27ac_anno)

mcclay_h3k27ac_genes <- unique(mcclay_h3k27ac_anno_df$SYMBOL)
mcclay_h3k27ac_genes <- mcclay_h3k27ac_genes[
  !is.na(mcclay_h3k27ac_genes) &
    mcclay_h3k27ac_genes != ""
]

forrest_h3k27ac_anno <- annotatePeak(
  forrest_TCF4_H3K27ac_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_h3k27ac_anno_df <- as.data.frame(forrest_h3k27ac_anno)

forrest_h3k27ac_genes <- unique(forrest_h3k27ac_anno_df$SYMBOL)
forrest_h3k27ac_genes <- forrest_h3k27ac_genes[
  !is.na(forrest_h3k27ac_genes) &
    forrest_h3k27ac_genes != ""
]


# Restrict gene sets to background universe
npc_all_genes_bg <- intersect(npc_all_genes, background_genes)
mcclay_all_genes_bg <- intersect(mcclay_all_genes, background_genes)
forrest_all_genes_bg <- intersect(forrest_all_genes, background_genes)

npc_h3k27ac_genes_bg <- intersect(npc_h3k27ac_genes, background_genes)
mcclay_h3k27ac_genes_bg <- intersect(mcclay_h3k27ac_genes, background_genes)
forrest_h3k27ac_genes_bg <- intersect(forrest_h3k27ac_genes, background_genes)

length(npc_h3k27ac_genes)
length(npc_h3k27ac_genes_bg)

length(mcclay_h3k27ac_genes)
length(mcclay_h3k27ac_genes_bg)

length(forrest_h3k27ac_genes)
length(forrest_h3k27ac_genes_bg)


# Gene-level summary table
h3k27ac_gene_overlap_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Total_TCF4_Genes = c(
    length(npc_all_genes),
    length(mcclay_all_genes),
    length(forrest_all_genes)
  ),
  Total_TCF4_Genes_In_Background = c(
    length(npc_all_genes_bg),
    length(mcclay_all_genes_bg),
    length(forrest_all_genes_bg)
  ),
  TCF4_H3K27ac_Overlap_Genes = c(
    length(npc_h3k27ac_genes),
    length(mcclay_h3k27ac_genes),
    length(forrest_h3k27ac_genes)
  ),
  TCF4_H3K27ac_Overlap_Genes_In_Background = c(
    length(npc_h3k27ac_genes_bg),
    length(mcclay_h3k27ac_genes_bg),
    length(forrest_h3k27ac_genes_bg)
  ),
  Background_Genes = length(background_genes)
)

h3k27ac_gene_overlap_summary$Percent_TCF4_Genes_With_H3K27ac <- round(
  h3k27ac_gene_overlap_summary$TCF4_H3K27ac_Overlap_Genes /
    h3k27ac_gene_overlap_summary$Total_TCF4_Genes * 100,
  3
)

h3k27ac_gene_overlap_summary

# Gene overlap between H3K27ac-overlapping TCF4 datasets
npc_mcclay_h3k27ac_genes <- intersect(
  npc_h3k27ac_genes_bg,
  mcclay_h3k27ac_genes_bg
)

npc_forrest_h3k27ac_genes <- intersect(
  npc_h3k27ac_genes_bg,
  forrest_h3k27ac_genes_bg
)

mcclay_forrest_h3k27ac_genes <- intersect(
  mcclay_h3k27ac_genes_bg,
  forrest_h3k27ac_genes_bg
)

all_three_h3k27ac_genes <- Reduce(
  intersect,
  list(
    npc_h3k27ac_genes_bg,
    mcclay_h3k27ac_genes_bg,
    forrest_h3k27ac_genes_bg
  )
)

length(npc_mcclay_h3k27ac_genes)
length(npc_forrest_h3k27ac_genes)
length(mcclay_forrest_h3k27ac_genes)
length(all_three_h3k27ac_genes)

# Fisher exact tests
# NPC vs McClay
a <- length(intersect(npc_h3k27ac_genes_bg, mcclay_h3k27ac_genes_bg))
b <- length(setdiff(npc_h3k27ac_genes_bg, mcclay_h3k27ac_genes_bg))
c <- length(setdiff(mcclay_h3k27ac_genes_bg, npc_h3k27ac_genes_bg))
d <- length(setdiff(background_genes, union(npc_h3k27ac_genes_bg, mcclay_h3k27ac_genes_bg)))

npc_mcclay_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
npc_mcclay_fisher <- fisher.test(npc_mcclay_table)

# NPC vs Forrest
a <- length(intersect(npc_h3k27ac_genes_bg, forrest_h3k27ac_genes_bg))
b <- length(setdiff(npc_h3k27ac_genes_bg, forrest_h3k27ac_genes_bg))
c <- length(setdiff(forrest_h3k27ac_genes_bg, npc_h3k27ac_genes_bg))
d <- length(setdiff(background_genes, union(npc_h3k27ac_genes_bg, forrest_h3k27ac_genes_bg)))

npc_forrest_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
npc_forrest_fisher <- fisher.test(npc_forrest_table)

# McClay vs Forrest
a <- length(intersect(mcclay_h3k27ac_genes_bg, forrest_h3k27ac_genes_bg))
b <- length(setdiff(mcclay_h3k27ac_genes_bg, forrest_h3k27ac_genes_bg))
c <- length(setdiff(forrest_h3k27ac_genes_bg, mcclay_h3k27ac_genes_bg))
d <- length(setdiff(background_genes, union(mcclay_h3k27ac_genes_bg, forrest_h3k27ac_genes_bg)))

mcclay_forrest_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
mcclay_forrest_fisher <- fisher.test(mcclay_forrest_table)

# Fisher summary table
h3k27ac_gene_fisher_summary <- data.frame(
  Comparison = c(
    "NPC_H3K27ac_genes_vs_McClay_H3K27ac_genes",
    "NPC_H3K27ac_genes_vs_Forrest_H3K27ac_genes",
    "McClay_H3K27ac_genes_vs_Forrest_H3K27ac_genes"
  ),
  Dataset_1_Genes = c(
    length(npc_h3k27ac_genes_bg),
    length(npc_h3k27ac_genes_bg),
    length(mcclay_h3k27ac_genes_bg)
  ),
  Dataset_2_Genes = c(
    length(mcclay_h3k27ac_genes_bg),
    length(forrest_h3k27ac_genes_bg),
    length(forrest_h3k27ac_genes_bg)
  ),
  Overlap_Genes = c(
    length(npc_mcclay_h3k27ac_genes),
    length(npc_forrest_h3k27ac_genes),
    length(mcclay_forrest_h3k27ac_genes)
  ),
  Percent_of_Dataset_1 = c(
    round(length(npc_mcclay_h3k27ac_genes) / length(npc_h3k27ac_genes_bg) * 100, 3),
    round(length(npc_forrest_h3k27ac_genes) / length(npc_h3k27ac_genes_bg) * 100, 3),
    round(length(mcclay_forrest_h3k27ac_genes) / length(mcclay_h3k27ac_genes_bg) * 100, 3)
  ),
  Percent_of_Dataset_2 = c(
    round(length(npc_mcclay_h3k27ac_genes) / length(mcclay_h3k27ac_genes_bg) * 100, 3),
    round(length(npc_forrest_h3k27ac_genes) / length(forrest_h3k27ac_genes_bg) * 100, 3),
    round(length(mcclay_forrest_h3k27ac_genes) / length(forrest_h3k27ac_genes_bg) * 100, 3)
  ),
  Background_Genes = length(background_genes),
  Fisher_P_Value = c(
    npc_mcclay_fisher$p.value,
    npc_forrest_fisher$p.value,
    mcclay_forrest_fisher$p.value
  ),
  Odds_Ratio = c(
    unname(npc_mcclay_fisher$estimate),
    unname(npc_forrest_fisher$estimate),
    unname(mcclay_forrest_fisher$estimate)
  )
)

h3k27ac_gene_fisher_summary

############################################################
# H3K4me3 peaks vs NPC TCF4 vs Dr. McClay's vs Forrest
# Genome build: hg38 / GRCh38
############################################################
# Import BED files
mcclay_peaks <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
npc_peaks <- import("NPC_hg38_summit_500bp.bed")
forrest_peak <- import("Forrest_hg38.bed")

# Import H3K4me3 BED file safely
H3K4me3_df <- read.table(
  "GSE176914_ENCFF907ZRF_replicated_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

H3K4me3_df <- H3K4me3_df[, 1:3]
colnames(H3K4me3_df) <- c("chr", "start", "end")

H3K4me3_df$start <- as.numeric(H3K4me3_df$start)
H3K4me3_df$end <- as.numeric(H3K4me3_df$end)

H3K4me3_df <- H3K4me3_df[
  !is.na(H3K4me3_df$chr) &
    !is.na(H3K4me3_df$start) &
    !is.na(H3K4me3_df$end) &
    H3K4me3_df$end > H3K4me3_df$start,
]

H3K4me3_peak <- GRanges(
  seqnames = H3K4me3_df$chr,
  ranges = IRanges(
    start = H3K4me3_df$start + 1,
    end = H3K4me3_df$end
  )
)

seqlevelsStyle(H3K4me3_peak) <- "UCSC"
H3K4me3_peak <- keepStandardChromosomes(H3K4me3_peak, pruning.mode = "coarse")
H3K4me3_peak <- sort(H3K4me3_peak)

length(H3K4me3_peak)

# Import background gene universe
background_df <- read.csv("final_ref_gene_2026_filtered.csv")

background_genes <- unique(background_df$gene)

background_genes <- background_genes[
  !is.na(background_genes) &
    background_genes != ""
]

length(background_genes)


# Clean all peak files
seqlevelsStyle(mcclay_peaks) <- "UCSC"
mcclay_peaks <- keepStandardChromosomes(mcclay_peaks, pruning.mode = "coarse")
mcclay_peaks <- sort(mcclay_peaks)

seqlevelsStyle(npc_peaks) <- "UCSC"
npc_peaks <- keepStandardChromosomes(npc_peaks, pruning.mode = "coarse")
npc_peaks <- sort(npc_peaks)

seqlevelsStyle(forrest_peak) <- "UCSC"
forrest_peak <- keepStandardChromosomes(forrest_peak, pruning.mode = "coarse")
forrest_peak <- sort(forrest_peak)

seqlevelsStyle(H3K4me3_peak) <- "UCSC"
H3K4me3_peak <- keepStandardChromosomes(H3K4me3_peak, pruning.mode = "coarse")
H3K4me3_peak <- sort(H3K4me3_peak)

length(mcclay_peaks)
length(npc_peaks)
length(forrest_peak)
length(H3K4me3_peak)


# Peak-level overlap: TCF4 datasets vs H3K4me3
# NPC TCF4 vs H3K4me3
npc_h3k4me3_hits <- findOverlaps(npc_peaks, H3K4me3_peak)

npc_h3k4me3_query_peaks <- unique(queryHits(npc_h3k4me3_hits))
npc_h3k4me3_subject_peaks <- unique(subjectHits(npc_h3k4me3_hits))

length(npc_h3k4me3_query_peaks)
length(npc_h3k4me3_subject_peaks)

# McClay TCF4 vs H3K4me3
mcclay_h3k4me3_hits <- findOverlaps(mcclay_peaks, H3K4me3_peak)

mcclay_h3k4me3_query_peaks <- unique(queryHits(mcclay_h3k4me3_hits))
mcclay_h3k4me3_subject_peaks <- unique(subjectHits(mcclay_h3k4me3_hits))

length(mcclay_h3k4me3_query_peaks)
length(mcclay_h3k4me3_subject_peaks)

# Forrest TCF4 vs H3K4me3
forrest_h3k4me3_hits <- findOverlaps(forrest_peak, H3K4me3_peak)

forrest_h3k4me3_query_peaks <- unique(queryHits(forrest_h3k4me3_hits))
forrest_h3k4me3_subject_peaks <- unique(subjectHits(forrest_h3k4me3_hits))

length(forrest_h3k4me3_query_peaks)
length(forrest_h3k4me3_subject_peaks)

# Peak-level overlap summary table
h3k4me3_peak_overlap_summary <- data.frame(
  Comparison = c(
    "NPC_TCF4_vs_H3K4me3",
    "McClay_TCF4_vs_H3K4me3",
    "Forrest_TCF4_vs_H3K4me3"
  ),
  TCF4_Total_Peaks = c(
    length(npc_peaks),
    length(mcclay_peaks),
    length(forrest_peak)
  ),
  H3K4me3_Total_Peaks = c(
    length(H3K4me3_peak),
    length(H3K4me3_peak),
    length(H3K4me3_peak)
  ),
  TCF4_Overlapping_Peaks = c(
    length(npc_h3k4me3_query_peaks),
    length(mcclay_h3k4me3_query_peaks),
    length(forrest_h3k4me3_query_peaks)
  ),
  H3K4me3_Overlapping_Peaks = c(
    length(npc_h3k4me3_subject_peaks),
    length(mcclay_h3k4me3_subject_peaks),
    length(forrest_h3k4me3_subject_peaks)
  )
)

h3k4me3_peak_overlap_summary$Percent_TCF4_Peaks_Overlapping_H3K4me3 <-
  round(
    h3k4me3_peak_overlap_summary$TCF4_Overlapping_Peaks /
      h3k4me3_peak_overlap_summary$TCF4_Total_Peaks * 100,
    3
  )

h3k4me3_peak_overlap_summary$Percent_H3K4me3_Peaks_Overlapping_TCF4 <-
  round(
    h3k4me3_peak_overlap_summary$H3K4me3_Overlapping_Peaks /
      h3k4me3_peak_overlap_summary$H3K4me3_Total_Peaks * 100,
    3
  )

h3k4me3_peak_overlap_summary


# Create H3K4me3-overlapping TCF4 peak objects
npc_TCF4_H3K4me3_peaks <- npc_peaks[npc_h3k4me3_query_peaks]

mcclay_TCF4_H3K4me3_peaks <- mcclay_peaks[mcclay_h3k4me3_query_peaks]

forrest_TCF4_H3K4me3_peaks <- forrest_peak[forrest_h3k4me3_query_peaks]

length(npc_TCF4_H3K4me3_peaks)
length(mcclay_TCF4_H3K4me3_peaks)
length(forrest_TCF4_H3K4me3_peaks)

# Gene-level analysis: H3K4me3-overlapping TCF4 genes
# Annotate all TCF4 peak sets to genes
npc_all_anno <- annotatePeak(
  npc_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_all_anno_df <- as.data.frame(npc_all_anno)

npc_all_genes <- unique(npc_all_anno_df$SYMBOL)
npc_all_genes <- npc_all_genes[!is.na(npc_all_genes) & npc_all_genes != ""]

mcclay_all_anno <- annotatePeak(
  mcclay_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_all_anno_df <- as.data.frame(mcclay_all_anno)

mcclay_all_genes <- unique(mcclay_all_anno_df$SYMBOL)
mcclay_all_genes <- mcclay_all_genes[!is.na(mcclay_all_genes) & mcclay_all_genes != ""]

forrest_all_anno <- annotatePeak(
  forrest_peak,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_all_anno_df <- as.data.frame(forrest_all_anno)

forrest_all_genes <- unique(forrest_all_anno_df$SYMBOL)
forrest_all_genes <- forrest_all_genes[!is.na(forrest_all_genes) & forrest_all_genes != ""]

# Annotate H3K4me3-overlapping TCF4 peaks to genes
npc_h3k4me3_anno <- annotatePeak(
  npc_TCF4_H3K4me3_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_h3k4me3_anno_df <- as.data.frame(npc_h3k4me3_anno)

npc_h3k4me3_genes <- unique(npc_h3k4me3_anno_df$SYMBOL)
npc_h3k4me3_genes <- npc_h3k4me3_genes[
  !is.na(npc_h3k4me3_genes) &
    npc_h3k4me3_genes != ""
]

mcclay_h3k4me3_anno <- annotatePeak(
  mcclay_TCF4_H3K4me3_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_h3k4me3_anno_df <- as.data.frame(mcclay_h3k4me3_anno)

mcclay_h3k4me3_genes <- unique(mcclay_h3k4me3_anno_df$SYMBOL)
mcclay_h3k4me3_genes <- mcclay_h3k4me3_genes[
  !is.na(mcclay_h3k4me3_genes) &
    mcclay_h3k4me3_genes != ""
]

forrest_h3k4me3_anno <- annotatePeak(
  forrest_TCF4_H3K4me3_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_h3k4me3_anno_df <- as.data.frame(forrest_h3k4me3_anno)

forrest_h3k4me3_genes <- unique(forrest_h3k4me3_anno_df$SYMBOL)
forrest_h3k4me3_genes <- forrest_h3k4me3_genes[
  !is.na(forrest_h3k4me3_genes) &
    forrest_h3k4me3_genes != ""]

# Restrict gene sets to background universe
npc_all_genes_bg <- intersect(npc_all_genes, background_genes)
mcclay_all_genes_bg <- intersect(mcclay_all_genes, background_genes)
forrest_all_genes_bg <- intersect(forrest_all_genes, background_genes)

npc_h3k4me3_genes_bg <- intersect(npc_h3k4me3_genes, background_genes)
mcclay_h3k4me3_genes_bg <- intersect(mcclay_h3k4me3_genes, background_genes)
forrest_h3k4me3_genes_bg <- intersect(forrest_h3k4me3_genes, background_genes)

length(npc_h3k4me3_genes)
length(npc_h3k4me3_genes_bg)

length(mcclay_h3k4me3_genes)
length(mcclay_h3k4me3_genes_bg)

length(forrest_h3k4me3_genes)
length(forrest_h3k4me3_genes_bg)


# Gene-level summary table
h3k4me3_gene_overlap_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Total_TCF4_Genes = c(
    length(npc_all_genes),
    length(mcclay_all_genes),
    length(forrest_all_genes)
  ),
  Total_TCF4_Genes_In_Background = c(
    length(npc_all_genes_bg),
    length(mcclay_all_genes_bg),
    length(forrest_all_genes_bg)
  ),
  TCF4_H3K4me3_Overlap_Genes = c(
    length(npc_h3k4me3_genes),
    length(mcclay_h3k4me3_genes),
    length(forrest_h3k4me3_genes)
  ),
  TCF4_H3K4me3_Overlap_Genes_In_Background = c(
    length(npc_h3k4me3_genes_bg),
    length(mcclay_h3k4me3_genes_bg),
    length(forrest_h3k4me3_genes_bg)
  ),
  Background_Genes = length(background_genes)
)

h3k4me3_gene_overlap_summary$Percent_TCF4_Genes_With_H3K4me3 <- round(
  h3k4me3_gene_overlap_summary$TCF4_H3K4me3_Overlap_Genes /
    h3k4me3_gene_overlap_summary$Total_TCF4_Genes * 100,
  3
)

h3k4me3_gene_overlap_summary


# Gene overlap between H3K4me3-overlapping TCF4 datasets
npc_mcclay_h3k4me3_genes <- intersect(
  npc_h3k4me3_genes_bg,
  mcclay_h3k4me3_genes_bg
)

npc_forrest_h3k4me3_genes <- intersect(
  npc_h3k4me3_genes_bg,
  forrest_h3k4me3_genes_bg
)

mcclay_forrest_h3k4me3_genes <- intersect(
  mcclay_h3k4me3_genes_bg,
  forrest_h3k4me3_genes_bg
)

all_three_h3k4me3_genes <- Reduce(
  intersect,
  list(
    npc_h3k4me3_genes_bg,
    mcclay_h3k4me3_genes_bg,
    forrest_h3k4me3_genes_bg
  )
)

length(npc_mcclay_h3k4me3_genes)
length(npc_forrest_h3k4me3_genes)
length(mcclay_forrest_h3k4me3_genes)
length(all_three_h3k4me3_genes)

############################################################
# H3K27me3 peaks vs NPC TCF4 vs Dr. McClay's vs Forrest
# Genome build: hg38 / GRCh38
############################################################
# Import BED files
mcclay_peaks <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
npc_peaks <- import("NPC_hg38_summit_500bp.bed")
forrest_peak <- import("Forrest_hg38.bed")

# Import H3K27me3 BED file safely
H3K27me3_df <- read.table(
  gzfile("ENCFF056AFA.bed.gz"),
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

H3K27me3_df <- H3K27me3_df[, 1:3]
colnames(H3K27me3_df) <- c("chr", "start", "end")

H3K27me3_df$start <- as.numeric(H3K27me3_df$start)
H3K27me3_df$end <- as.numeric(H3K27me3_df$end)

H3K27me3_df <- H3K27me3_df[
  !is.na(H3K27me3_df$chr) &
    !is.na(H3K27me3_df$start) &
    !is.na(H3K27me3_df$end) &
    H3K27me3_df$end > H3K27me3_df$start,
]

H3K27me3_peak <- GRanges(
  seqnames = H3K27me3_df$chr,
  ranges = IRanges(
    start = H3K27me3_df$start + 1,
    end = H3K27me3_df$end
  )
)

seqlevelsStyle(H3K27me3_peak) <- "UCSC"
H3K27me3_peak <- keepStandardChromosomes(H3K27me3_peak, pruning.mode = "coarse")
H3K27me3_peak <- sort(H3K27me3_peak)

length(H3K27me3_peak)

# Import background gene universe
background_df <- read.csv("final_ref_gene_2026_filtered.csv")

background_genes <- unique(background_df$gene)

background_genes <- background_genes[
  !is.na(background_genes) &
    background_genes != ""
]

length(background_genes)

# Clean all peak files
seqlevelsStyle(mcclay_peaks) <- "UCSC"
mcclay_peaks <- keepStandardChromosomes(mcclay_peaks, pruning.mode = "coarse")
mcclay_peaks <- sort(mcclay_peaks)

seqlevelsStyle(npc_peaks) <- "UCSC"
npc_peaks <- keepStandardChromosomes(npc_peaks, pruning.mode = "coarse")
npc_peaks <- sort(npc_peaks)

seqlevelsStyle(forrest_peak) <- "UCSC"
forrest_peak <- keepStandardChromosomes(forrest_peak, pruning.mode = "coarse")
forrest_peak <- sort(forrest_peak)

seqlevelsStyle(H3K27me3_peak) <- "UCSC"
H3K27me3_peak <- keepStandardChromosomes(H3K27me3_peak, pruning.mode = "coarse")
H3K27me3_peak <- sort(H3K27me3_peak)

length(mcclay_peaks)
length(npc_peaks)
length(forrest_peak)
length(H3K27me3_peak)

# Peak-level overlap: TCF4 datasets vs H3K27me3

# NPC TCF4 vs H3K27me3
npc_h3k27me3_hits <- findOverlaps(npc_peaks, H3K27me3_peak)

npc_h3k27me3_query_peaks <- unique(queryHits(npc_h3k27me3_hits))
npc_h3k27me3_subject_peaks <- unique(subjectHits(npc_h3k27me3_hits))

length(npc_h3k27me3_query_peaks)
length(npc_h3k27me3_subject_peaks)

# McClay TCF4 vs H3K27me3
mcclay_h3k27me3_hits <- findOverlaps(mcclay_peaks, H3K27me3_peak)

mcclay_h3k27me3_query_peaks <- unique(queryHits(mcclay_h3k27me3_hits))
mcclay_h3k27me3_subject_peaks <- unique(subjectHits(mcclay_h3k27me3_hits))

length(mcclay_h3k27me3_query_peaks)
length(mcclay_h3k27me3_subject_peaks)

# Forrest TCF4 vs H3K27me3
forrest_h3k27me3_hits <- findOverlaps(forrest_peak, H3K27me3_peak)

forrest_h3k27me3_query_peaks <- unique(queryHits(forrest_h3k27me3_hits))
forrest_h3k27me3_subject_peaks <- unique(subjectHits(forrest_h3k27me3_hits))

length(forrest_h3k27me3_query_peaks)
length(forrest_h3k27me3_subject_peaks)

# Peak-level overlap summary table
h3k27me3_peak_overlap_summary <- data.frame(
  Comparison = c(
    "NPC_TCF4_vs_H3K27me3",
    "McClay_TCF4_vs_H3K27me3",
    "Forrest_TCF4_vs_H3K27me3"
  ),
  TCF4_Total_Peaks = c(
    length(npc_peaks),
    length(mcclay_peaks),
    length(forrest_peak)
  ),
  H3K27me3_Total_Peaks = c(
    length(H3K27me3_peak),
    length(H3K27me3_peak),
    length(H3K27me3_peak)
  ),
  TCF4_Overlapping_Peaks = c(
    length(npc_h3k27me3_query_peaks),
    length(mcclay_h3k27me3_query_peaks),
    length(forrest_h3k27me3_query_peaks)
  ),
  H3K27me3_Overlapping_Peaks = c(
    length(npc_h3k27me3_subject_peaks),
    length(mcclay_h3k27me3_subject_peaks),
    length(forrest_h3k27me3_subject_peaks)
  )
)

h3k27me3_peak_overlap_summary$Percent_TCF4_Peaks_Overlapping_H3K27me3 <-
  round(
    h3k27me3_peak_overlap_summary$TCF4_Overlapping_Peaks /
      h3k27me3_peak_overlap_summary$TCF4_Total_Peaks * 100,
    3
  )

h3k27me3_peak_overlap_summary$Percent_H3K27me3_Peaks_Overlapping_TCF4 <-
  round(
    h3k27me3_peak_overlap_summary$H3K27me3_Overlapping_Peaks /
      h3k27me3_peak_overlap_summary$H3K27me3_Total_Peaks * 100,
    3
  )

h3k27me3_peak_overlap_summary

# Create H3K27me3-overlapping TCF4 peak objects
npc_TCF4_H3K27me3_peaks <- npc_peaks[npc_h3k27me3_query_peaks]

mcclay_TCF4_H3K27me3_peaks <- mcclay_peaks[mcclay_h3k27me3_query_peaks]

forrest_TCF4_H3K27me3_peaks <- forrest_peak[forrest_h3k27me3_query_peaks]

length(npc_TCF4_H3K27me3_peaks)
length(mcclay_TCF4_H3K27me3_peaks)
length(forrest_TCF4_H3K27me3_peaks)

# Gene-level analysis: H3K27me3-overlapping TCF4 genes
# Annotate H3K27me3-overlapping TCF4 peaks to genes

npc_h3k27me3_anno <- annotatePeak(
  npc_TCF4_H3K27me3_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_h3k27me3_anno_df <- as.data.frame(npc_h3k27me3_anno)

npc_h3k27me3_genes <- unique(npc_h3k27me3_anno_df$SYMBOL)
npc_h3k27me3_genes <- npc_h3k27me3_genes[
  !is.na(npc_h3k27me3_genes) &
    npc_h3k27me3_genes != ""
]


mcclay_h3k27me3_anno <- annotatePeak(
  mcclay_TCF4_H3K27me3_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_h3k27me3_anno_df <- as.data.frame(mcclay_h3k27me3_anno)

mcclay_h3k27me3_genes <- unique(mcclay_h3k27me3_anno_df$SYMBOL)
mcclay_h3k27me3_genes <- mcclay_h3k27me3_genes[
  !is.na(mcclay_h3k27me3_genes) &
    mcclay_h3k27me3_genes != ""
]


forrest_h3k27me3_anno <- annotatePeak(
  forrest_TCF4_H3K27me3_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_h3k27me3_anno_df <- as.data.frame(forrest_h3k27me3_anno)

forrest_h3k27me3_genes <- unique(forrest_h3k27me3_anno_df$SYMBOL)
forrest_h3k27me3_genes <- forrest_h3k27me3_genes[
  !is.na(forrest_h3k27me3_genes) &
    forrest_h3k27me3_genes != ""
]


# Simple gene count summary

h3k27me3_gene_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4_H3K27me3",
    "McClay_TCF4_H3K27me3",
    "Forrest_TCF4_H3K27me3"
  ),
  Gene_Count = c(
    length(npc_h3k27me3_genes),
    length(mcclay_h3k27me3_genes),
    length(forrest_h3k27me3_genes)
  )
)

h3k27me3_gene_summary


# Pairwise gene overlaps

npc_mcclay_h3k27me3_genes <- intersect(
  npc_h3k27me3_genes,
  mcclay_h3k27me3_genes
)

npc_forrest_h3k27me3_genes <- intersect(
  npc_h3k27me3_genes,
  forrest_h3k27me3_genes
)

mcclay_forrest_h3k27me3_genes <- intersect(
  mcclay_h3k27me3_genes,
  forrest_h3k27me3_genes
)

all_three_h3k27me3_genes <- Reduce(
  intersect,
  list(
    npc_h3k27me3_genes,
    mcclay_h3k27me3_genes,
    forrest_h3k27me3_genes
  )
)


# Gene overlap summary table

h3k27me3_gene_overlap_simple <- data.frame(
  Comparison = c(
    "NPC vs McClay",
    "NPC vs Forrest",
    "McClay vs Forrest",
    "NPC vs McClay vs Forrest"
  ),
  Overlap_Genes = c(
    length(npc_mcclay_h3k27me3_genes),
    length(npc_forrest_h3k27me3_genes),
    length(mcclay_forrest_h3k27me3_genes),
    length(all_three_h3k27me3_genes)
  )
)

h3k27me3_gene_overlap_simple

############################################################
# EP300 peaks vs NPC TCF4 vs Dr. McClay's vs Forrest
# Genome build: hg38 / GRCh38
############################################################
# Import BED files
mcclay_peaks <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
npc_peaks <- import("NPC_hg38_summit_500bp.bed")
forrest_peak <- import("Forrest_hg38.bed")

# Import EP300 optimal BED file safely
EP300_optimal_df <- read.table(
  "GSE127584_ENCFF459ARL_optimal_idr_thresholded_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

EP300_optimal_df <- EP300_optimal_df[, 1:3]
colnames(EP300_optimal_df) <- c("chr", "start", "end")

EP300_optimal_df$start <- as.numeric(EP300_optimal_df$start)
EP300_optimal_df$end <- as.numeric(EP300_optimal_df$end)

EP300_optimal_df <- EP300_optimal_df[
  !is.na(EP300_optimal_df$chr) &
    !is.na(EP300_optimal_df$start) &
    !is.na(EP300_optimal_df$end) &
    EP300_optimal_df$end > EP300_optimal_df$start,
]

EP300_optimal_peak <- GRanges(
  seqnames = EP300_optimal_df$chr,
  ranges = IRanges(
    start = EP300_optimal_df$start + 1,
    end = EP300_optimal_df$end
  )
)

seqlevelsStyle(EP300_optimal_peak) <- "UCSC"
EP300_optimal_peak <- keepStandardChromosomes(EP300_optimal_peak, pruning.mode = "coarse")
EP300_optimal_peak <- sort(EP300_optimal_peak)

length(EP300_optimal_peak)

# Import EP300 conservative BED file safely
EP300_conservative_df <- read.table(
  "GSE127584_ENCFF747GDL_conservative_idr_thresholded_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

EP300_conservative_df <- EP300_conservative_df[, 1:3]
colnames(EP300_conservative_df) <- c("chr", "start", "end")

EP300_conservative_df$start <- as.numeric(EP300_conservative_df$start)
EP300_conservative_df$end <- as.numeric(EP300_conservative_df$end)

EP300_conservative_df <- EP300_conservative_df[
  !is.na(EP300_conservative_df$chr) &
    !is.na(EP300_conservative_df$start) &
    !is.na(EP300_conservative_df$end) &
    EP300_conservative_df$end > EP300_conservative_df$start,
]

EP300_conservative_peak <- GRanges(
  seqnames = EP300_conservative_df$chr,
  ranges = IRanges(
    start = EP300_conservative_df$start + 1,
    end = EP300_conservative_df$end
  )
)

seqlevelsStyle(EP300_conservative_peak) <- "UCSC"
EP300_conservative_peak <- keepStandardChromosomes(EP300_conservative_peak, pruning.mode = "coarse")
EP300_conservative_peak <- sort(EP300_conservative_peak)

length(EP300_conservative_peak)

# Clean all peak files
seqlevelsStyle(mcclay_peaks) <- "UCSC"
mcclay_peaks <- keepStandardChromosomes(mcclay_peaks, pruning.mode = "coarse")
mcclay_peaks <- sort(mcclay_peaks)

seqlevelsStyle(npc_peaks) <- "UCSC"
npc_peaks <- keepStandardChromosomes(npc_peaks, pruning.mode = "coarse")
npc_peaks <- sort(npc_peaks)

seqlevelsStyle(forrest_peak) <- "UCSC"
forrest_peak <- keepStandardChromosomes(forrest_peak, pruning.mode = "coarse")
forrest_peak <- sort(forrest_peak)

seqlevelsStyle(EP300_optimal_peak) <- "UCSC"
EP300_optimal_peak <- keepStandardChromosomes(EP300_optimal_peak, pruning.mode = "coarse")
EP300_optimal_peak <- sort(EP300_optimal_peak)

seqlevelsStyle(EP300_conservative_peak) <- "UCSC"
EP300_conservative_peak <- keepStandardChromosomes(EP300_conservative_peak, pruning.mode = "coarse")
EP300_conservative_peak <- sort(EP300_conservative_peak)

length(mcclay_peaks)
length(npc_peaks)
length(forrest_peak)
length(EP300_optimal_peak)
length(EP300_conservative_peak)

# Peak-level overlap: TCF4 datasets vs EP300 optimal
npc_ep300_optimal_hits <- findOverlaps(npc_peaks, EP300_optimal_peak)
npc_ep300_optimal_query_peaks <- unique(queryHits(npc_ep300_optimal_hits))
npc_ep300_optimal_subject_peaks <- unique(subjectHits(npc_ep300_optimal_hits))

mcclay_ep300_optimal_hits <- findOverlaps(mcclay_peaks, EP300_optimal_peak)
mcclay_ep300_optimal_query_peaks <- unique(queryHits(mcclay_ep300_optimal_hits))
mcclay_ep300_optimal_subject_peaks <- unique(subjectHits(mcclay_ep300_optimal_hits))

forrest_ep300_optimal_hits <- findOverlaps(forrest_peak, EP300_optimal_peak)
forrest_ep300_optimal_query_peaks <- unique(queryHits(forrest_ep300_optimal_hits))
forrest_ep300_optimal_subject_peaks <- unique(subjectHits(forrest_ep300_optimal_hits))

length(npc_ep300_optimal_query_peaks)
length(mcclay_ep300_optimal_query_peaks)
length(forrest_ep300_optimal_query_peaks)

# Peak-level overlap: TCF4 datasets vs EP300 conservative
npc_ep300_conservative_hits <- findOverlaps(npc_peaks, EP300_conservative_peak)
npc_ep300_conservative_query_peaks <- unique(queryHits(npc_ep300_conservative_hits))
npc_ep300_conservative_subject_peaks <- unique(subjectHits(npc_ep300_conservative_hits))

mcclay_ep300_conservative_hits <- findOverlaps(mcclay_peaks, EP300_conservative_peak)
mcclay_ep300_conservative_query_peaks <- unique(queryHits(mcclay_ep300_conservative_hits))
mcclay_ep300_conservative_subject_peaks <- unique(subjectHits(mcclay_ep300_conservative_hits))

forrest_ep300_conservative_hits <- findOverlaps(forrest_peak, EP300_conservative_peak)
forrest_ep300_conservative_query_peaks <- unique(queryHits(forrest_ep300_conservative_hits))
forrest_ep300_conservative_subject_peaks <- unique(subjectHits(forrest_ep300_conservative_hits))

length(npc_ep300_conservative_query_peaks)
length(mcclay_ep300_conservative_query_peaks)
length(forrest_ep300_conservative_query_peaks)

# Peak-level overlap summary table
ep300_peak_overlap_summary <- data.frame(
  EP300_Set = c(
    "Optimal",
    "Optimal",
    "Optimal",
    "Conservative",
    "Conservative",
    "Conservative"
  ),
  Comparison = c(
    "NPC_TCF4_vs_EP300",
    "McClay_TCF4_vs_EP300",
    "Forrest_TCF4_vs_EP300",
    "NPC_TCF4_vs_EP300",
    "McClay_TCF4_vs_EP300",
    "Forrest_TCF4_vs_EP300"
  ),
  TCF4_Total_Peaks = c(
    length(npc_peaks),
    length(mcclay_peaks),
    length(forrest_peak),
    length(npc_peaks),
    length(mcclay_peaks),
    length(forrest_peak)
  ),
  EP300_Total_Peaks = c(
    length(EP300_optimal_peak),
    length(EP300_optimal_peak),
    length(EP300_optimal_peak),
    length(EP300_conservative_peak),
    length(EP300_conservative_peak),
    length(EP300_conservative_peak)
  ),
  TCF4_Overlapping_Peaks = c(
    length(npc_ep300_optimal_query_peaks),
    length(mcclay_ep300_optimal_query_peaks),
    length(forrest_ep300_optimal_query_peaks),
    length(npc_ep300_conservative_query_peaks),
    length(mcclay_ep300_conservative_query_peaks),
    length(forrest_ep300_conservative_query_peaks)
  ),
  EP300_Overlapping_Peaks = c(
    length(npc_ep300_optimal_subject_peaks),
    length(mcclay_ep300_optimal_subject_peaks),
    length(forrest_ep300_optimal_subject_peaks),
    length(npc_ep300_conservative_subject_peaks),
    length(mcclay_ep300_conservative_subject_peaks),
    length(forrest_ep300_conservative_subject_peaks)
  )
)

ep300_peak_overlap_summary$Percent_TCF4_Peaks_Overlapping_EP300 <- round(
  ep300_peak_overlap_summary$TCF4_Overlapping_Peaks /
    ep300_peak_overlap_summary$TCF4_Total_Peaks * 100,
  3
)

ep300_peak_overlap_summary$Percent_EP300_Peaks_Overlapping_TCF4 <- round(
  ep300_peak_overlap_summary$EP300_Overlapping_Peaks /
    ep300_peak_overlap_summary$EP300_Total_Peaks * 100,
  3
)

ep300_peak_overlap_summary

# Create EP300-overlapping TCF4 peak objects
npc_TCF4_EP300_optimal_peaks <- npc_peaks[npc_ep300_optimal_query_peaks]
mcclay_TCF4_EP300_optimal_peaks <- mcclay_peaks[mcclay_ep300_optimal_query_peaks]
forrest_TCF4_EP300_optimal_peaks <- forrest_peak[forrest_ep300_optimal_query_peaks]

npc_TCF4_EP300_conservative_peaks <- npc_peaks[npc_ep300_conservative_query_peaks]
mcclay_TCF4_EP300_conservative_peaks <- mcclay_peaks[mcclay_ep300_conservative_query_peaks]
forrest_TCF4_EP300_conservative_peaks <- forrest_peak[forrest_ep300_conservative_query_peaks]

length(npc_TCF4_EP300_optimal_peaks)
length(mcclay_TCF4_EP300_optimal_peaks)
length(forrest_TCF4_EP300_optimal_peaks)

length(npc_TCF4_EP300_conservative_peaks)
length(mcclay_TCF4_EP300_conservative_peaks)
length(forrest_TCF4_EP300_conservative_peaks)

# Gene-level analysis: EP300-overlapping TCF4 genes
npc_ep300_optimal_anno <- annotatePeak(
  npc_TCF4_EP300_optimal_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_ep300_optimal_anno_df <- as.data.frame(npc_ep300_optimal_anno)
npc_ep300_optimal_genes <- unique(npc_ep300_optimal_anno_df$SYMBOL)
npc_ep300_optimal_genes <- npc_ep300_optimal_genes[
  !is.na(npc_ep300_optimal_genes) &
    npc_ep300_optimal_genes != ""
]

mcclay_ep300_optimal_anno <- annotatePeak(
  mcclay_TCF4_EP300_optimal_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_ep300_optimal_anno_df <- as.data.frame(mcclay_ep300_optimal_anno)
mcclay_ep300_optimal_genes <- unique(mcclay_ep300_optimal_anno_df$SYMBOL)
mcclay_ep300_optimal_genes <- mcclay_ep300_optimal_genes[
  !is.na(mcclay_ep300_optimal_genes) &
    mcclay_ep300_optimal_genes != ""
]

forrest_ep300_optimal_anno <- annotatePeak(
  forrest_TCF4_EP300_optimal_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_ep300_optimal_anno_df <- as.data.frame(forrest_ep300_optimal_anno)
forrest_ep300_optimal_genes <- unique(forrest_ep300_optimal_anno_df$SYMBOL)
forrest_ep300_optimal_genes <- forrest_ep300_optimal_genes[
  !is.na(forrest_ep300_optimal_genes) &
    forrest_ep300_optimal_genes != ""
]

npc_ep300_conservative_anno <- annotatePeak(
  npc_TCF4_EP300_conservative_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_ep300_conservative_anno_df <- as.data.frame(npc_ep300_conservative_anno)
npc_ep300_conservative_genes <- unique(npc_ep300_conservative_anno_df$SYMBOL)
npc_ep300_conservative_genes <- npc_ep300_conservative_genes[
  !is.na(npc_ep300_conservative_genes) &
    npc_ep300_conservative_genes != ""
]

mcclay_ep300_conservative_anno <- annotatePeak(
  mcclay_TCF4_EP300_conservative_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_ep300_conservative_anno_df <- as.data.frame(mcclay_ep300_conservative_anno)
mcclay_ep300_conservative_genes <- unique(mcclay_ep300_conservative_anno_df$SYMBOL)
mcclay_ep300_conservative_genes <- mcclay_ep300_conservative_genes[
  !is.na(mcclay_ep300_conservative_genes) &
    mcclay_ep300_conservative_genes != ""
]

forrest_ep300_conservative_anno <- annotatePeak(
  forrest_TCF4_EP300_conservative_peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_ep300_conservative_anno_df <- as.data.frame(forrest_ep300_conservative_anno)
forrest_ep300_conservative_genes <- unique(forrest_ep300_conservative_anno_df$SYMBOL)
forrest_ep300_conservative_genes <- forrest_ep300_conservative_genes[
  !is.na(forrest_ep300_conservative_genes) &
    forrest_ep300_conservative_genes != ""
]

# Gene count summary
ep300_gene_summary <- data.frame(
  EP300_Set = c(
    "Optimal",
    "Optimal",
    "Optimal",
    "Conservative",
    "Conservative",
    "Conservative"
  ),
  Dataset = c(
    "NPC_TCF4_EP300",
    "McClay_TCF4_EP300",
    "Forrest_TCF4_EP300",
    "NPC_TCF4_EP300",
    "McClay_TCF4_EP300",
    "Forrest_TCF4_EP300"
  ),
  Gene_Count = c(
    length(npc_ep300_optimal_genes),
    length(mcclay_ep300_optimal_genes),
    length(forrest_ep300_optimal_genes),
    length(npc_ep300_conservative_genes),
    length(mcclay_ep300_conservative_genes),
    length(forrest_ep300_conservative_genes)
  )
)

ep300_gene_summary

# Pairwise gene overlaps: EP300 optimal
npc_mcclay_ep300_optimal_genes <- intersect(
  npc_ep300_optimal_genes,
  mcclay_ep300_optimal_genes
)

npc_forrest_ep300_optimal_genes <- intersect(
  npc_ep300_optimal_genes,
  forrest_ep300_optimal_genes
)

mcclay_forrest_ep300_optimal_genes <- intersect(
  mcclay_ep300_optimal_genes,
  forrest_ep300_optimal_genes
)

all_three_ep300_optimal_genes <- Reduce(
  intersect,
  list(
    npc_ep300_optimal_genes,
    mcclay_ep300_optimal_genes,
    forrest_ep300_optimal_genes
  )
)

# Pairwise gene overlaps: EP300 conservative
npc_mcclay_ep300_conservative_genes <- intersect(
  npc_ep300_conservative_genes,
  mcclay_ep300_conservative_genes
)

npc_forrest_ep300_conservative_genes <- intersect(
  npc_ep300_conservative_genes,
  forrest_ep300_conservative_genes
)

mcclay_forrest_ep300_conservative_genes <- intersect(
  mcclay_ep300_conservative_genes,
  forrest_ep300_conservative_genes
)

all_three_ep300_conservative_genes <- Reduce(
  intersect,
  list(
    npc_ep300_conservative_genes,
    mcclay_ep300_conservative_genes,
    forrest_ep300_conservative_genes
  )
)

# Gene overlap summary table
ep300_gene_overlap_simple <- data.frame(
  EP300_Set = c(
    "Optimal",
    "Optimal",
    "Optimal",
    "Optimal",
    "Conservative",
    "Conservative",
    "Conservative",
    "Conservative"
  ),
  Comparison = c(
    "NPC vs McClay",
    "NPC vs Forrest",
    "McClay vs Forrest",
    "NPC vs McClay vs Forrest",
    "NPC vs McClay",
    "NPC vs Forrest",
    "McClay vs Forrest",
    "NPC vs McClay vs Forrest"
  ),
  Overlap_Genes = c(
    length(npc_mcclay_ep300_optimal_genes),
    length(npc_forrest_ep300_optimal_genes),
    length(mcclay_forrest_ep300_optimal_genes),
    length(all_three_ep300_optimal_genes),
    length(npc_mcclay_ep300_conservative_genes),
    length(npc_forrest_ep300_conservative_genes),
    length(mcclay_forrest_ep300_conservative_genes),
    length(all_three_ep300_conservative_genes)
  )
)

ep300_gene_overlap_simple


############################################################
# H3K27ac and CHD8 vs NPC TCF4 vs Dr. McClay's vs Forrest
# Genome build: hg38 / GRCh38
############################################################
# Import and clean H3K27ac
k3k27ac <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "27AC_peaks",
  skip = 2,
  col_names = c("chr", "start", "end")
)

k3k27ac <- k3k27ac %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  filter(!is.na(start), !is.na(end), end > start)

k3k27ac_peak <- GRanges(
  seqnames = k3k27ac$chr,
  ranges = IRanges(start = k3k27ac$start + 1, end = k3k27ac$end)
)

seqlevelsStyle(k3k27ac_peak) <- "UCSC"
k3k27ac_peak <- keepStandardChromosomes(k3k27ac_peak, pruning.mode = "coarse")
k3k27ac_peak <- sort(k3k27ac_peak)

length(k3k27ac_peak)

# Import and clean Neuron CHD8
chd8_neurons <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "neuronCHD8_peaks",
  skip = 2,
  col_names = c("chr", "start", "end")
)

chd8_neurons <- chd8_neurons %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start), !is.na(end), end > start)

chd8_neurons_peak <- GRanges(
  seqnames = chd8_neurons$chr,
  ranges = IRanges(start = chd8_neurons$start + 1, end = chd8_neurons$end)
)

seqlevelsStyle(chd8_neurons_peak) <- "UCSC"
chd8_neurons_peak <- keepStandardChromosomes(chd8_neurons_peak, pruning.mode = "coarse")
chd8_neurons_peak <- sort(chd8_neurons_peak)

length(chd8_neurons_peak)

# Import and clean Brain CHD8
chd8_brain <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "brainCHD8_peaks",
  skip = 2,
  col_names = c("chr", "start", "end")
)

chd8_brain <- chd8_brain %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start), !is.na(end), end > start)

chd8_brain_peak <- GRanges(
  seqnames = chd8_brain$chr,
  ranges = IRanges(start = chd8_brain$start + 1, end = chd8_brain$end)
)

seqlevelsStyle(chd8_brain_peak) <- "UCSC"
chd8_brain_peak <- keepStandardChromosomes(chd8_brain_peak, pruning.mode = "coarse")
chd8_brain_peak <- sort(chd8_brain_peak)

length(chd8_brain_peak)

# Import and clean ATAC-seq
atac <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "ATAC_Seq_peaks",
  skip = 3,
  col_names = c("chr", "start", "end")
)

atac <- atac %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start), !is.na(end), end > start)

atac_peak <- GRanges(
  seqnames = atac$chr,
  ranges = IRanges(start = atac$start + 1, end = atac$end)
)

seqlevelsStyle(atac_peak) <- "UCSC"
atac_peak <- keepStandardChromosomes(atac_peak, pruning.mode = "coarse")
atac_peak <- sort(atac_peak)

length(atac_peak)

# Import TCF4 peak files
mcclay_peaks <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
npc_peaks <- import("NPC_hg38_summit_500bp.bed")
forrest_peak <- import("Forrest_hg38.bed")

# Clean TCF4 peak files
seqlevelsStyle(mcclay_peaks) <- "UCSC"
mcclay_peaks <- keepStandardChromosomes(mcclay_peaks, pruning.mode = "coarse")
mcclay_peaks <- sort(mcclay_peaks)

seqlevelsStyle(npc_peaks) <- "UCSC"
npc_peaks <- keepStandardChromosomes(npc_peaks, pruning.mode = "coarse")
npc_peaks <- sort(npc_peaks)

seqlevelsStyle(forrest_peak) <- "UCSC"
forrest_peak <- keepStandardChromosomes(forrest_peak, pruning.mode = "coarse")
forrest_peak <- sort(forrest_peak)

length(mcclay_peaks)
length(npc_peaks)
length(forrest_peak)

# Peak-level overlap: NPC TCF4
npc_gse79965_h3k27ac_hits <- findOverlaps(npc_peaks, k3k27ac_peak)
npc_gse79965_h3k27ac_query <- unique(queryHits(npc_gse79965_h3k27ac_hits))
npc_gse79965_h3k27ac_subject <- unique(subjectHits(npc_gse79965_h3k27ac_hits))

npc_chd8_neuron_hits <- findOverlaps(npc_peaks, chd8_neurons_peak)
npc_chd8_neuron_query <- unique(queryHits(npc_chd8_neuron_hits))
npc_chd8_neuron_subject <- unique(subjectHits(npc_chd8_neuron_hits))

npc_chd8_brain_hits <- findOverlaps(npc_peaks, chd8_brain_peak)
npc_chd8_brain_query <- unique(queryHits(npc_chd8_brain_hits))
npc_chd8_brain_subject <- unique(subjectHits(npc_chd8_brain_hits))

npc_atac_hits <- findOverlaps(npc_peaks, atac_peak)
npc_atac_query <- unique(queryHits(npc_atac_hits))
npc_atac_subject <- unique(subjectHits(npc_atac_hits))

# Peak-level overlap: McClay TCF4
mcclay_gse79965_h3k27ac_hits <- findOverlaps(mcclay_peaks, k3k27ac_peak)
mcclay_gse79965_h3k27ac_query <- unique(queryHits(mcclay_gse79965_h3k27ac_hits))
mcclay_gse79965_h3k27ac_subject <- unique(subjectHits(mcclay_gse79965_h3k27ac_hits))

mcclay_chd8_neuron_hits <- findOverlaps(mcclay_peaks, chd8_neurons_peak)
mcclay_chd8_neuron_query <- unique(queryHits(mcclay_chd8_neuron_hits))
mcclay_chd8_neuron_subject <- unique(subjectHits(mcclay_chd8_neuron_hits))

mcclay_chd8_brain_hits <- findOverlaps(mcclay_peaks, chd8_brain_peak)
mcclay_chd8_brain_query <- unique(queryHits(mcclay_chd8_brain_hits))
mcclay_chd8_brain_subject <- unique(subjectHits(mcclay_chd8_brain_hits))

mcclay_atac_hits <- findOverlaps(mcclay_peaks, atac_peak)
mcclay_atac_query <- unique(queryHits(mcclay_atac_hits))
mcclay_atac_subject <- unique(subjectHits(mcclay_atac_hits))

# Peak-level overlap: Forrest TCF4
forrest_gse79965_h3k27ac_hits <- findOverlaps(forrest_peak, k3k27ac_peak)
forrest_gse79965_h3k27ac_query <- unique(queryHits(forrest_gse79965_h3k27ac_hits))
forrest_gse79965_h3k27ac_subject <- unique(subjectHits(forrest_gse79965_h3k27ac_hits))

forrest_chd8_neuron_hits <- findOverlaps(forrest_peak, chd8_neurons_peak)
forrest_chd8_neuron_query <- unique(queryHits(forrest_chd8_neuron_hits))
forrest_chd8_neuron_subject <- unique(subjectHits(forrest_chd8_neuron_hits))

forrest_chd8_brain_hits <- findOverlaps(forrest_peak, chd8_brain_peak)
forrest_chd8_brain_query <- unique(queryHits(forrest_chd8_brain_hits))
forrest_chd8_brain_subject <- unique(subjectHits(forrest_chd8_brain_hits))

forrest_atac_hits <- findOverlaps(forrest_peak, atac_peak)
forrest_atac_query <- unique(queryHits(forrest_atac_hits))
forrest_atac_subject <- unique(subjectHits(forrest_atac_hits))

# Summary table
gse79965_tcf4_overlap_summary <- data.frame(
  Feature = c(
    "Neuron_H3K27ac",
    "Neuron_CHD8",
    "Brain_CHD8",
    "Neuron_ATAC",
    "Neuron_H3K27ac",
    "Neuron_CHD8",
    "Brain_CHD8",
    "Neuron_ATAC",
    "Neuron_H3K27ac",
    "Neuron_CHD8",
    "Brain_CHD8",
    "Neuron_ATAC"
  ),
  TCF4_Dataset = c(
    "NPC_TCF4",
    "NPC_TCF4",
    "NPC_TCF4",
    "NPC_TCF4",
    "McClay_TCF4",
    "McClay_TCF4",
    "McClay_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4",
    "Forrest_TCF4",
    "Forrest_TCF4",
    "Forrest_TCF4"
  ),
  TCF4_Total_Peaks = c(
    length(npc_peaks),
    length(npc_peaks),
    length(npc_peaks),
    length(npc_peaks),
    length(mcclay_peaks),
    length(mcclay_peaks),
    length(mcclay_peaks),
    length(forrest_peak),
    length(forrest_peak),
    length(forrest_peak),
    length(forrest_peak),
    length(forrest_peak)
  ),
  Feature_Total_Peaks = c(
    length(k3k27ac_peak),
    length(chd8_neurons_peak),
    length(chd8_brain_peak),
    length(atac_peak),
    length(k3k27ac_peak),
    length(chd8_neurons_peak),
    length(chd8_brain_peak),
    length(atac_peak),
    length(k3k27ac_peak),
    length(chd8_neurons_peak),
    length(chd8_brain_peak),
    length(atac_peak)
  ),
  TCF4_Overlapping_Peaks = c(
    length(npc_gse79965_h3k27ac_query),
    length(npc_chd8_neuron_query),
    length(npc_chd8_brain_query),
    length(npc_atac_query),
    length(mcclay_gse79965_h3k27ac_query),
    length(mcclay_chd8_neuron_query),
    length(mcclay_chd8_brain_query),
    length(mcclay_atac_query),
    length(forrest_gse79965_h3k27ac_query),
    length(forrest_chd8_neuron_query),
    length(forrest_chd8_brain_query),
    length(forrest_atac_query)
  ),
  Feature_Overlapping_Peaks = c(
    length(npc_gse79965_h3k27ac_subject),
    length(npc_chd8_neuron_subject),
    length(npc_chd8_brain_subject),
    length(npc_atac_subject),
    length(mcclay_gse79965_h3k27ac_subject),
    length(mcclay_chd8_neuron_subject),
    length(mcclay_chd8_brain_subject),
    length(mcclay_atac_subject),
    length(forrest_gse79965_h3k27ac_subject),
    length(forrest_chd8_neuron_subject),
    length(forrest_chd8_brain_subject),
    length(forrest_atac_subject)
  )
)

gse79965_tcf4_overlap_summary$Percent_TCF4_Peaks_Overlapping_Feature <- round(
  gse79965_tcf4_overlap_summary$TCF4_Overlapping_Peaks /
    gse79965_tcf4_overlap_summary$TCF4_Total_Peaks * 100,
  3
)

gse79965_tcf4_overlap_summary$Percent_Feature_Peaks_Overlapping_TCF4 <- round(
  gse79965_tcf4_overlap_summary$Feature_Overlapping_Peaks /
    gse79965_tcf4_overlap_summary$Feature_Total_Peaks * 100,
  3
)

gse79965_tcf4_overlap_summary

############################################################
# GSE79965 H3K27ac / CHD8 / ATAC-seq vs hg19 TCF4 peak sets
# Testing whether GSE79965 coordinates behave more like hg19
############################################################
# Import hg19 TCF4 peak files
npc_hg19 <- import("NPC_ab21_idr_500bp_summit.bed")
forrest_hg19 <- import("Forrest_TCF4_IDR_sorted_hg19.bed")
mcclay_hg19_df <- read.csv("McClay_TCF4_11322_consensus_hg19_annotated.csv")

colnames(mcclay_hg19_df)

# If columns are seqnames/start/end
mcclay_hg19 <- GRanges(
  seqnames = mcclay_hg19_df$chrom,
  ranges = IRanges(
    start = mcclay_hg19_df$start,
    end = mcclay_hg19_df$end
  )
)

# Clean hg19 TCF4 peak files
seqlevelsStyle(npc_hg19) <- "UCSC"
npc_hg19 <- keepStandardChromosomes(npc_hg19, pruning.mode = "coarse")
npc_hg19 <- sort(npc_hg19)

seqlevelsStyle(mcclay_hg19) <- "UCSC"
mcclay_hg19 <- keepStandardChromosomes(mcclay_hg19, pruning.mode = "coarse")
mcclay_hg19 <- sort(mcclay_hg19)

seqlevelsStyle(forrest_hg19) <- "UCSC"
forrest_hg19 <- keepStandardChromosomes(forrest_hg19, pruning.mode = "coarse")
forrest_hg19 <- sort(forrest_hg19)

length(npc_hg19)
length(mcclay_hg19)
length(forrest_hg19)

# Import GSE79965 H3K27ac peaks
gse_h3k27ac <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "27AC_peaks",
  skip = 2,
  col_names = c("chr", "start", "end")
)

gse_h3k27ac <- gse_h3k27ac %>%
  dplyr::filter(!is.na(chr), chr != "Chromosome") %>%
  dplyr::mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  dplyr::filter(!is.na(start), !is.na(end), end > start)

gse_h3k27ac_peak <- GRanges(
  seqnames = gse_h3k27ac$chr,
  ranges = IRanges(
    start = gse_h3k27ac$start + 1,
    end = gse_h3k27ac$end
  )
)

seqlevelsStyle(gse_h3k27ac_peak) <- "UCSC"
gse_h3k27ac_peak <- keepStandardChromosomes(gse_h3k27ac_peak, pruning.mode = "coarse")
gse_h3k27ac_peak <- sort(gse_h3k27ac_peak)

length(gse_h3k27ac_peak)

# Import GSE79965 neuron CHD8 peaks
gse_chd8_neuron <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "neuronCHD8_peaks",
  skip = 2,
  col_names = c("chr", "start", "end")
)

gse_chd8_neuron <- gse_chd8_neuron %>%
  dplyr::filter(!is.na(chr), chr != "Chromosome") %>%
  dplyr::mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  dplyr::filter(!is.na(start), !is.na(end), end > start)

gse_chd8_neuron_peak <- GRanges(
  seqnames = gse_chd8_neuron$chr,
  ranges = IRanges(
    start = gse_chd8_neuron$start + 1,
    end = gse_chd8_neuron$end
  )
)

seqlevelsStyle(gse_chd8_neuron_peak) <- "UCSC"
gse_chd8_neuron_peak <- keepStandardChromosomes(gse_chd8_neuron_peak, pruning.mode = "coarse")
gse_chd8_neuron_peak <- sort(gse_chd8_neuron_peak)

length(gse_chd8_neuron_peak)

# Import GSE79965 brain CHD8 peaks
gse_chd8_brain <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "brainCHD8_peaks",
  skip = 2,
  col_names = c("chr", "start", "end")
)

gse_chd8_brain <- gse_chd8_brain %>%
  dplyr::filter(!is.na(chr), chr != "Chromosome") %>%
  dplyr::mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  dplyr::filter(!is.na(start), !is.na(end), end > start)

gse_chd8_brain_peak <- GRanges(
  seqnames = gse_chd8_brain$chr,
  ranges = IRanges(
    start = gse_chd8_brain$start + 1,
    end = gse_chd8_brain$end
  )
)

seqlevelsStyle(gse_chd8_brain_peak) <- "UCSC"
gse_chd8_brain_peak <- keepStandardChromosomes(gse_chd8_brain_peak, pruning.mode = "coarse")
gse_chd8_brain_peak <- sort(gse_chd8_brain_peak)

length(gse_chd8_brain_peak)

# Import GSE79965 ATAC-seq peaks
gse_atac <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "ATAC_Seq_peaks",
  skip = 3,
  col_names = c("chr", "start", "end")
)

gse_atac <- gse_atac %>%
  dplyr::filter(!is.na(chr), chr != "Chromosome") %>%
  dplyr::mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  dplyr::filter(!is.na(start), !is.na(end), end > start)

gse_atac_peak <- GRanges(
  seqnames = gse_atac$chr,
  ranges = IRanges(
    start = gse_atac$start + 1,
    end = gse_atac$end
  )
)

seqlevelsStyle(gse_atac_peak) <- "UCSC"
gse_atac_peak <- keepStandardChromosomes(gse_atac_peak, pruning.mode = "coarse")
gse_atac_peak <- sort(gse_atac_peak)

length(gse_atac_peak)

# NPC hg19 overlaps
npc_hg19_h3k27ac_hits <- findOverlaps(npc_hg19, gse_h3k27ac_peak)
npc_hg19_chd8_neuron_hits <- findOverlaps(npc_hg19, gse_chd8_neuron_peak)
npc_hg19_chd8_brain_hits <- findOverlaps(npc_hg19, gse_chd8_brain_peak)
npc_hg19_atac_hits <- findOverlaps(npc_hg19, gse_atac_peak)

# McClay hg19 overlaps
mcclay_hg19_h3k27ac_hits <- findOverlaps(mcclay_hg19, gse_h3k27ac_peak)
mcclay_hg19_chd8_neuron_hits <- findOverlaps(mcclay_hg19, gse_chd8_neuron_peak)
mcclay_hg19_chd8_brain_hits <- findOverlaps(mcclay_hg19, gse_chd8_brain_peak)
mcclay_hg19_atac_hits <- findOverlaps(mcclay_hg19, gse_atac_peak)

# Forrest hg19 overlaps
forrest_hg19_h3k27ac_hits <- findOverlaps(forrest_hg19, gse_h3k27ac_peak)
forrest_hg19_chd8_neuron_hits <- findOverlaps(forrest_hg19, gse_chd8_neuron_peak)
forrest_hg19_chd8_brain_hits <- findOverlaps(forrest_hg19, gse_chd8_brain_peak)
forrest_hg19_atac_hits <- findOverlaps(forrest_hg19, gse_atac_peak)

# Summary table
gse79965_hg19_overlap_summary <- data.frame(
  Feature = c(
    "GSE79965_Neuron_H3K27ac",
    "GSE79965_Neuron_CHD8",
    "GSE79965_Brain_CHD8",
    "GSE79965_Neuron_ATAC",
    "GSE79965_Neuron_H3K27ac",
    "GSE79965_Neuron_CHD8",
    "GSE79965_Brain_CHD8",
    "GSE79965_Neuron_ATAC",
    "GSE79965_Neuron_H3K27ac",
    "GSE79965_Neuron_CHD8",
    "GSE79965_Brain_CHD8",
    "GSE79965_Neuron_ATAC"
  ),
  TCF4_Dataset = c(
    "NPC_TCF4_hg19",
    "NPC_TCF4_hg19",
    "NPC_TCF4_hg19",
    "NPC_TCF4_hg19",
    "McClay_TCF4_hg19",
    "McClay_TCF4_hg19",
    "McClay_TCF4_hg19",
    "McClay_TCF4_hg19",
    "Forrest_TCF4_hg19",
    "Forrest_TCF4_hg19",
    "Forrest_TCF4_hg19",
    "Forrest_TCF4_hg19"
  ),
  TCF4_Total_Peaks = c(
    length(npc_hg19),
    length(npc_hg19),
    length(npc_hg19),
    length(npc_hg19),
    length(mcclay_hg19),
    length(mcclay_hg19),
    length(mcclay_hg19),
    length(mcclay_hg19),
    length(forrest_hg19),
    length(forrest_hg19),
    length(forrest_hg19),
    length(forrest_hg19)
  ),
  Feature_Total_Peaks = c(
    length(gse_h3k27ac_peak),
    length(gse_chd8_neuron_peak),
    length(gse_chd8_brain_peak),
    length(gse_atac_peak),
    length(gse_h3k27ac_peak),
    length(gse_chd8_neuron_peak),
    length(gse_chd8_brain_peak),
    length(gse_atac_peak),
    length(gse_h3k27ac_peak),
    length(gse_chd8_neuron_peak),
    length(gse_chd8_brain_peak),
    length(gse_atac_peak)
  ),
  TCF4_Overlapping_Peaks = c(
    length(unique(queryHits(npc_hg19_h3k27ac_hits))),
    length(unique(queryHits(npc_hg19_chd8_neuron_hits))),
    length(unique(queryHits(npc_hg19_chd8_brain_hits))),
    length(unique(queryHits(npc_hg19_atac_hits))),
    length(unique(queryHits(mcclay_hg19_h3k27ac_hits))),
    length(unique(queryHits(mcclay_hg19_chd8_neuron_hits))),
    length(unique(queryHits(mcclay_hg19_chd8_brain_hits))),
    length(unique(queryHits(mcclay_hg19_atac_hits))),
    length(unique(queryHits(forrest_hg19_h3k27ac_hits))),
    length(unique(queryHits(forrest_hg19_chd8_neuron_hits))),
    length(unique(queryHits(forrest_hg19_chd8_brain_hits))),
    length(unique(queryHits(forrest_hg19_atac_hits)))
  ),
  Feature_Overlapping_Peaks = c(
    length(unique(subjectHits(npc_hg19_h3k27ac_hits))),
    length(unique(subjectHits(npc_hg19_chd8_neuron_hits))),
    length(unique(subjectHits(npc_hg19_chd8_brain_hits))),
    length(unique(subjectHits(npc_hg19_atac_hits))),
    length(unique(subjectHits(mcclay_hg19_h3k27ac_hits))),
    length(unique(subjectHits(mcclay_hg19_chd8_neuron_hits))),
    length(unique(subjectHits(mcclay_hg19_chd8_brain_hits))),
    length(unique(subjectHits(mcclay_hg19_atac_hits))),
    length(unique(subjectHits(forrest_hg19_h3k27ac_hits))),
    length(unique(subjectHits(forrest_hg19_chd8_neuron_hits))),
    length(unique(subjectHits(forrest_hg19_chd8_brain_hits))),
    length(unique(subjectHits(forrest_hg19_atac_hits)))
  )
)

gse79965_hg19_overlap_summary$Percent_TCF4_Peaks_Overlapping_Feature <- round(
  gse79965_hg19_overlap_summary$TCF4_Overlapping_Peaks /
    gse79965_hg19_overlap_summary$TCF4_Total_Peaks * 100,
  3
)

gse79965_hg19_overlap_summary$Percent_Feature_Peaks_Overlapping_TCF4 <- round(
  gse79965_hg19_overlap_summary$Feature_Overlapping_Peaks /
    gse79965_hg19_overlap_summary$Feature_Total_Peaks * 100,
  3
)

gse79965_hg19_overlap_summary

############################################################
# GSE79965 Hi-C TAD Analysis with TCF4 and Epigenomic Features
# Genome build: hg19 / GRCh37
############################################################

# The goal of this analysis is not simple peak overlap.
# TADs are large chromatin domains.
# Here we use TADs as regulatory neighborhoods.

# Input objects expected from previous steps:
# npc_hg19
# mcclay_hg19
# forrest_hg19
# gse_h3k27ac_peak
# gse_chd8_neuron_peak
# gse_chd8_brain_peak
# gse_atac_peak


# Import Hi-C TADs from GSE79965
hic_df <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "Hi-C",
  skip = 2,
  col_names = c("chr", "start", "end")
)

hic_df <- hic_df %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  filter(!is.na(start), !is.na(end), end > start)

hic_tads <- GRanges(
  seqnames = hic_df$chr,
  ranges = IRanges(
    start = hic_df$start + 1,
    end = hic_df$end
  )
)

seqlevelsStyle(hic_tads) <- "UCSC"
hic_tads <- keepStandardChromosomes(hic_tads, pruning.mode = "coarse")
hic_tads <- sort(hic_tads)

length(hic_tads)


# Count how many peaks from each dataset fall inside each TAD
tad_summary <- data.frame(
  TAD_ID = paste0("TAD_", seq_along(hic_tads)),
  Chromosome = as.character(seqnames(hic_tads)),
  Start = start(hic_tads),
  End = end(hic_tads),
  Width_bp = width(hic_tads),

  NPC_TCF4_Peaks = countOverlaps(hic_tads, npc_hg19),
  McClay_TCF4_Peaks = countOverlaps(hic_tads, mcclay_hg19),
  Forrest_TCF4_Peaks = countOverlaps(hic_tads, forrest_hg19),

  H3K27ac_Peaks = countOverlaps(hic_tads, gse_h3k27ac_peak),
  Neuron_CHD8_Peaks = countOverlaps(hic_tads, gse_chd8_neuron_peak),
  Brain_CHD8_Peaks = countOverlaps(hic_tads, gse_chd8_brain_peak),
  ATAC_Peaks = countOverlaps(hic_tads, gse_atac_peak)
)

head(tad_summary)


# Convert peak counts to TRUE/FALSE presence columns
tad_summary$Has_NPC_TCF4 <- tad_summary$NPC_TCF4_Peaks > 0
tad_summary$Has_McClay_TCF4 <- tad_summary$McClay_TCF4_Peaks > 0
tad_summary$Has_Forrest_TCF4 <- tad_summary$Forrest_TCF4_Peaks > 0

tad_summary$Has_H3K27ac <- tad_summary$H3K27ac_Peaks > 0
tad_summary$Has_ATAC <- tad_summary$ATAC_Peaks > 0
tad_summary$Has_Neuron_CHD8 <- tad_summary$Neuron_CHD8_Peaks > 0
tad_summary$Has_Brain_CHD8 <- tad_summary$Brain_CHD8_Peaks > 0

tad_summary$Has_Any_CHD8 <-
  tad_summary$Has_Neuron_CHD8 | tad_summary$Has_Brain_CHD8


# Define biologically meaningful TAD classes
tad_summary$NPC_TCF4_TAD <-
  tad_summary$Has_NPC_TCF4

tad_summary$NPC_Active_TAD <-
  tad_summary$Has_NPC_TCF4 &
  tad_summary$Has_H3K27ac &
  tad_summary$Has_ATAC

tad_summary$NPC_CHD8_TAD <-
  tad_summary$Has_NPC_TCF4 &
  tad_summary$Has_Any_CHD8

tad_summary$NPC_Active_CHD8_TAD <-
  tad_summary$Has_NPC_TCF4 &
  tad_summary$Has_H3K27ac &
  tad_summary$Has_ATAC &
  tad_summary$Has_Any_CHD8

tad_summary$Shared_TCF4_TAD <-
  tad_summary$Has_NPC_TCF4 &
  tad_summary$Has_McClay_TCF4 &
  tad_summary$Has_Forrest_TCF4


# Summary of TAD categories
tad_category_summary <- data.frame(
  Category = c(
    "Total Hi-C TADs",
    "TADs with NPC TCF4",
    "TADs with McClay TCF4",
    "TADs with Forrest TCF4",
    "TADs with H3K27ac",
    "TADs with ATAC",
    "TADs with Neuron CHD8",
    "TADs with Brain CHD8",
    "TADs with Any CHD8",
    "NPC TCF4 + H3K27ac + ATAC TADs",
    "NPC TCF4 + CHD8 TADs",
    "NPC TCF4 + H3K27ac + ATAC + CHD8 TADs",
    "TADs shared by NPC, McClay, and Forrest TCF4"
  ),
  TAD_Count = c(
    nrow(tad_summary),
    sum(tad_summary$Has_NPC_TCF4),
    sum(tad_summary$Has_McClay_TCF4),
    sum(tad_summary$Has_Forrest_TCF4),
    sum(tad_summary$Has_H3K27ac),
    sum(tad_summary$Has_ATAC),
    sum(tad_summary$Has_Neuron_CHD8),
    sum(tad_summary$Has_Brain_CHD8),
    sum(tad_summary$Has_Any_CHD8),
    sum(tad_summary$NPC_Active_TAD),
    sum(tad_summary$NPC_CHD8_TAD),
    sum(tad_summary$NPC_Active_CHD8_TAD),
    sum(tad_summary$Shared_TCF4_TAD)
  )
)

tad_category_summary$Percent_of_All_TADs <- round(
  tad_category_summary$TAD_Count / nrow(tad_summary) * 100,
  3
)

tad_category_summary


# Compare TCF4 datasets at TAD level
npc_tads <- tad_summary$TAD_ID[tad_summary$Has_NPC_TCF4]
mcclay_tads <- tad_summary$TAD_ID[tad_summary$Has_McClay_TCF4]
forrest_tads <- tad_summary$TAD_ID[tad_summary$Has_Forrest_TCF4]

npc_mcclay_tads <- intersect(npc_tads, mcclay_tads)
npc_forrest_tads <- intersect(npc_tads, forrest_tads)
mcclay_forrest_tads <- intersect(mcclay_tads, forrest_tads)

all_three_tcf4_tads <- Reduce(
  intersect,
  list(npc_tads, mcclay_tads, forrest_tads)
)

tcf4_tad_overlap_summary <- data.frame(
  Comparison = c(
    "NPC TCF4 TADs vs McClay TCF4 TADs",
    "NPC TCF4 TADs vs Forrest TCF4 TADs",
    "McClay TCF4 TADs vs Forrest TCF4 TADs",
    "NPC vs McClay vs Forrest TCF4 TADs"
  ),
  Shared_TADs = c(
    length(npc_mcclay_tads),
    length(npc_forrest_tads),
    length(mcclay_forrest_tads),
    length(all_three_tcf4_tads)
  ),
  Percent_of_First_TAD_Set = c(
    round(length(npc_mcclay_tads) / length(npc_tads) * 100, 3),
    round(length(npc_forrest_tads) / length(npc_tads) * 100, 3),
    round(length(mcclay_forrest_tads) / length(mcclay_tads) * 100, 3),
    round(length(all_three_tcf4_tads) / length(npc_tads) * 100, 3)
  ),
  Percent_of_Second_TAD_Set = c(
    round(length(npc_mcclay_tads) / length(mcclay_tads) * 100, 3),
    round(length(npc_forrest_tads) / length(forrest_tads) * 100, 3),
    round(length(mcclay_forrest_tads) / length(forrest_tads) * 100, 3),
    NA
  )
)

tcf4_tad_overlap_summary


# Define regulatory scores for each TCF4 dataset
tad_summary$NPC_Regulatory_Score <-
  tad_summary$NPC_TCF4_Peaks +
  tad_summary$H3K27ac_Peaks +
  tad_summary$ATAC_Peaks +
  tad_summary$Neuron_CHD8_Peaks +
  tad_summary$Brain_CHD8_Peaks

tad_summary$McClay_Regulatory_Score <-
  tad_summary$McClay_TCF4_Peaks +
  tad_summary$H3K27ac_Peaks +
  tad_summary$ATAC_Peaks +
  tad_summary$Neuron_CHD8_Peaks +
  tad_summary$Brain_CHD8_Peaks

tad_summary$Forrest_Regulatory_Score <-
  tad_summary$Forrest_TCF4_Peaks +
  tad_summary$H3K27ac_Peaks +
  tad_summary$ATAC_Peaks +
  tad_summary$Neuron_CHD8_Peaks +
  tad_summary$Brain_CHD8_Peaks


# Identify strongest NPC active CHD8 TADs
npc_active_chd8_tads <- tad_summary %>%
  filter(NPC_Active_CHD8_TAD == TRUE) %>%
  arrange(
    desc(NPC_TCF4_Peaks),
    desc(H3K27ac_Peaks),
    desc(ATAC_Peaks),
    desc(Neuron_CHD8_Peaks + Brain_CHD8_Peaks)
  )

head(npc_active_chd8_tads, 20)


# Rank all TADs by combined regulatory score
top_npc_regulatory_tads <- tad_summary %>%
  arrange(desc(NPC_Regulatory_Score))

top_mcclay_regulatory_tads <- tad_summary %>%
  arrange(desc(McClay_Regulatory_Score))

top_forrest_regulatory_tads <- tad_summary %>%
  arrange(desc(Forrest_Regulatory_Score))

head(top_npc_regulatory_tads, 20)
head(top_mcclay_regulatory_tads, 20)
head(top_forrest_regulatory_tads, 20)


# Rank TADs separately by each feature
top_npc_tcf4_tads <- tad_summary %>%
  arrange(desc(NPC_TCF4_Peaks)) %>%
  slice_head(n = 100)

top_mcclay_tcf4_tads <- tad_summary %>%
  arrange(desc(McClay_TCF4_Peaks)) %>%
  slice_head(n = 100)

top_forrest_tcf4_tads <- tad_summary %>%
  arrange(desc(Forrest_TCF4_Peaks)) %>%
  slice_head(n = 100)

top_h3k27ac_tads <- tad_summary %>%
  arrange(desc(H3K27ac_Peaks)) %>%
  slice_head(n = 100)

top_atac_tads <- tad_summary %>%
  arrange(desc(ATAC_Peaks)) %>%
  slice_head(n = 100)

top_neuron_chd8_tads <- tad_summary %>%
  arrange(desc(Neuron_CHD8_Peaks)) %>%
  slice_head(n = 100)

top_brain_chd8_tads <- tad_summary %>%
  arrange(desc(Brain_CHD8_Peaks)) %>%
  slice_head(n = 100)


# Compare top 100 TAD overlaps for all TCF4 datasets
top100_tcf4_feature_overlap_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC_TCF4", "NPC_TCF4", "NPC_TCF4", "NPC_TCF4",
    "McClay_TCF4", "McClay_TCF4", "McClay_TCF4", "McClay_TCF4",
    "Forrest_TCF4", "Forrest_TCF4", "Forrest_TCF4", "Forrest_TCF4"
  ),
  Feature = c(
    "H3K27ac", "ATAC", "Neuron_CHD8", "Brain_CHD8",
    "H3K27ac", "ATAC", "Neuron_CHD8", "Brain_CHD8",
    "H3K27ac", "ATAC", "Neuron_CHD8", "Brain_CHD8"
  ),
  Shared_Top100_TADs = c(
    length(intersect(top_npc_tcf4_tads$TAD_ID, top_h3k27ac_tads$TAD_ID)),
    length(intersect(top_npc_tcf4_tads$TAD_ID, top_atac_tads$TAD_ID)),
    length(intersect(top_npc_tcf4_tads$TAD_ID, top_neuron_chd8_tads$TAD_ID)),
    length(intersect(top_npc_tcf4_tads$TAD_ID, top_brain_chd8_tads$TAD_ID)),

    length(intersect(top_mcclay_tcf4_tads$TAD_ID, top_h3k27ac_tads$TAD_ID)),
    length(intersect(top_mcclay_tcf4_tads$TAD_ID, top_atac_tads$TAD_ID)),
    length(intersect(top_mcclay_tcf4_tads$TAD_ID, top_neuron_chd8_tads$TAD_ID)),
    length(intersect(top_mcclay_tcf4_tads$TAD_ID, top_brain_chd8_tads$TAD_ID)),

    length(intersect(top_forrest_tcf4_tads$TAD_ID, top_h3k27ac_tads$TAD_ID)),
    length(intersect(top_forrest_tcf4_tads$TAD_ID, top_atac_tads$TAD_ID)),
    length(intersect(top_forrest_tcf4_tads$TAD_ID, top_neuron_chd8_tads$TAD_ID)),
    length(intersect(top_forrest_tcf4_tads$TAD_ID, top_brain_chd8_tads$TAD_ID))
  )
)

top100_tcf4_feature_overlap_summary$Percent_Shared <-
  top100_tcf4_feature_overlap_summary$Shared_Top100_TADs

top100_tcf4_feature_overlap_summary


# Compare top 100 TAD overlaps among epigenomic features
top100_feature_feature_overlap_summary <- data.frame(
  Comparison = c(
    "H3K27ac vs ATAC",
    "H3K27ac vs Neuron CHD8",
    "H3K27ac vs Brain CHD8",
    "ATAC vs Neuron CHD8",
    "ATAC vs Brain CHD8",
    "Neuron CHD8 vs Brain CHD8"
  ),
  Shared_Top100_TADs = c(
    length(intersect(top_h3k27ac_tads$TAD_ID, top_atac_tads$TAD_ID)),
    length(intersect(top_h3k27ac_tads$TAD_ID, top_neuron_chd8_tads$TAD_ID)),
    length(intersect(top_h3k27ac_tads$TAD_ID, top_brain_chd8_tads$TAD_ID)),
    length(intersect(top_atac_tads$TAD_ID, top_neuron_chd8_tads$TAD_ID)),
    length(intersect(top_atac_tads$TAD_ID, top_brain_chd8_tads$TAD_ID)),
    length(intersect(top_neuron_chd8_tads$TAD_ID, top_brain_chd8_tads$TAD_ID))
  )
)

top100_feature_feature_overlap_summary$Percent_Shared <-
  top100_feature_feature_overlap_summary$Shared_Top100_TADs

top100_feature_feature_overlap_summary


# TADs that are top 100 for all active/regulatory features with NPC TCF4
top_all_npc_active_tads <- Reduce(
  intersect,
  list(
    top_npc_tcf4_tads$TAD_ID,
    top_h3k27ac_tads$TAD_ID,
    top_atac_tads$TAD_ID,
    top_neuron_chd8_tads$TAD_ID,
    top_brain_chd8_tads$TAD_ID
  )
)

length(top_all_npc_active_tads)

top_all_npc_active_tad_table <- tad_summary %>%
  filter(TAD_ID %in% top_all_npc_active_tads) %>%
  arrange(desc(NPC_Regulatory_Score))

top_all_npc_active_tad_table


# TADs that are top 100 for all active/regulatory features with McClay TCF4
top_all_mcclay_active_tads <- Reduce(
  intersect,
  list(
    top_mcclay_tcf4_tads$TAD_ID,
    top_h3k27ac_tads$TAD_ID,
    top_atac_tads$TAD_ID,
    top_neuron_chd8_tads$TAD_ID,
    top_brain_chd8_tads$TAD_ID
  )
)

length(top_all_mcclay_active_tads)

top_all_mcclay_active_tad_table <- tad_summary %>%
  filter(TAD_ID %in% top_all_mcclay_active_tads) %>%
  arrange(desc(McClay_Regulatory_Score))

top_all_mcclay_active_tad_table


# TADs that are top 100 for all active/regulatory features with Forrest TCF4
top_all_forrest_active_tads <- Reduce(
  intersect,
  list(
    top_forrest_tcf4_tads$TAD_ID,
    top_h3k27ac_tads$TAD_ID,
    top_atac_tads$TAD_ID,
    top_neuron_chd8_tads$TAD_ID,
    top_brain_chd8_tads$TAD_ID
  )
)

length(top_all_forrest_active_tads)

top_all_forrest_active_tad_table <- tad_summary %>%
  filter(TAD_ID %in% top_all_forrest_active_tads) %>%
  arrange(desc(Forrest_Regulatory_Score))

top_all_forrest_active_tad_table


# Spearman correlation between peak counts across TADs
tad_count_correlation <- cor(
  tad_summary[, c(
    "NPC_TCF4_Peaks",
    "McClay_TCF4_Peaks",
    "Forrest_TCF4_Peaks",
    "H3K27ac_Peaks",
    "ATAC_Peaks",
    "Neuron_CHD8_Peaks",
    "Brain_CHD8_Peaks"
  )],
  method = "spearman"
)

round(tad_count_correlation, 3)


# TCF4-feature correlation summary
tcf4_feature_correlation_summary <- data.frame(
  TCF4_Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  H3K27ac = c(
    tad_count_correlation["NPC_TCF4_Peaks", "H3K27ac_Peaks"],
    tad_count_correlation["McClay_TCF4_Peaks", "H3K27ac_Peaks"],
    tad_count_correlation["Forrest_TCF4_Peaks", "H3K27ac_Peaks"]
  ),
  ATAC = c(
    tad_count_correlation["NPC_TCF4_Peaks", "ATAC_Peaks"],
    tad_count_correlation["McClay_TCF4_Peaks", "ATAC_Peaks"],
    tad_count_correlation["Forrest_TCF4_Peaks", "ATAC_Peaks"]
  ),
  Neuron_CHD8 = c(
    tad_count_correlation["NPC_TCF4_Peaks", "Neuron_CHD8_Peaks"],
    tad_count_correlation["McClay_TCF4_Peaks", "Neuron_CHD8_Peaks"],
    tad_count_correlation["Forrest_TCF4_Peaks", "Neuron_CHD8_Peaks"]
  ),
  Brain_CHD8 = c(
    tad_count_correlation["NPC_TCF4_Peaks", "Brain_CHD8_Peaks"],
    tad_count_correlation["McClay_TCF4_Peaks", "Brain_CHD8_Peaks"],
    tad_count_correlation["Forrest_TCF4_Peaks", "Brain_CHD8_Peaks"]
  )
)

tcf4_feature_correlation_summary[, 2:5] <- round(
  tcf4_feature_correlation_summary[, 2:5],
  3
)

tcf4_feature_correlation_summary

# Check the TAD_1375 (chr19:57.81 Mb - 58.98 Mb) out to see any important gene
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Define TAD_1375
tad_1375 <- GRanges(
  seqnames = "chr19",
  ranges = IRanges(
    start = 57810001,
    end = 58980000
  )
)

# Get hg19 gene coordinates
hg19_genes <- genes(txdb)

seqlevelsStyle(hg19_genes) <- "UCSC"
hg19_genes <- keepStandardChromosomes(hg19_genes, pruning.mode = "coarse")

# Find genes inside or overlapping this TAD
tad_1375_gene_hits <- findOverlaps(hg19_genes, tad_1375)

tad_1375_genes <- hg19_genes[unique(queryHits(tad_1375_gene_hits))]

# Convert Entrez IDs to gene symbols
tad_1375_entrez <- names(tad_1375_genes)

tad_1375_symbols <- mapIds(
  org.Hs.eg.db,
  keys = tad_1375_entrez,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Make clean table
tad_1375_gene_table <- data.frame(
  Entrez_ID = tad_1375_entrez,
  Gene_Symbol = as.character(tad_1375_symbols),
  Chromosome = as.character(seqnames(tad_1375_genes)),
  Start = start(tad_1375_genes),
  End = end(tad_1375_genes)
)

tad_1375_gene_table <- tad_1375_gene_table[
  !is.na(tad_1375_gene_table$Gene_Symbol),
]

tad_1375_gene_table


############################################################
# RNA-seq integration with TCF4 peaks and GSE79965 Hi-C TADs
# Genome build: hg19 / GRCh37
############################################################
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Import RNA-seq data
rnaseq_raw <- read_xlsx(
  "GSE79965_processedData.xlsx",
  sheet = "RNA-SEQ",
  col_names = FALSE
)

colnames(rnaseq_raw) <- c("Gene_Symbol", "FPKM")

rnaseq_df <- rnaseq_raw %>%
  filter(!is.na(Gene_Symbol), Gene_Symbol != "Symbols") %>%
  mutate(FPKM = as.numeric(FPKM)) %>%
  filter(!is.na(FPKM))

head(rnaseq_df)
summary(rnaseq_df$FPKM)

# Define expressed genes
expressed_genes_fpkm_0 <- rnaseq_df$Gene_Symbol[rnaseq_df$FPKM > 0]
expressed_genes_fpkm_0_5 <- rnaseq_df$Gene_Symbol[rnaseq_df$FPKM > 0.5]
expressed_genes_fpkm_1 <- rnaseq_df$Gene_Symbol[rnaseq_df$FPKM > 1]

length(expressed_genes_fpkm_0)
length(expressed_genes_fpkm_0_5)
length(expressed_genes_fpkm_1)

# Annotate TCF4 peaks to genes
npc_anno <- annotatePeak(
  npc_hg19,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

npc_anno_df <- as.data.frame(npc_anno)

npc_genes <- unique(npc_anno_df$SYMBOL)
npc_genes <- npc_genes[!is.na(npc_genes) & npc_genes != ""]

mcclay_anno <- annotatePeak(
  mcclay_hg19,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_anno_df <- as.data.frame(mcclay_anno)

mcclay_genes <- unique(mcclay_anno_df$SYMBOL)
mcclay_genes <- mcclay_genes[!is.na(mcclay_genes) & mcclay_genes != ""]

forrest_anno <- annotatePeak(
  forrest_hg19,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_anno_df <- as.data.frame(forrest_anno)

forrest_genes <- unique(forrest_anno_df$SYMBOL)
forrest_genes <- forrest_genes[!is.na(forrest_genes) & forrest_genes != ""]

length(npc_genes)
length(mcclay_genes)
length(forrest_genes)

# Expression overlap summary for TCF4-associated genes
npc_expressed_0 <- intersect(npc_genes, expressed_genes_fpkm_0)
npc_expressed_0_5 <- intersect(npc_genes, expressed_genes_fpkm_0_5)
npc_expressed_1 <- intersect(npc_genes, expressed_genes_fpkm_1)

mcclay_expressed_0 <- intersect(mcclay_genes, expressed_genes_fpkm_0)
mcclay_expressed_0_5 <- intersect(mcclay_genes, expressed_genes_fpkm_0_5)
mcclay_expressed_1 <- intersect(mcclay_genes, expressed_genes_fpkm_1)

forrest_expressed_0 <- intersect(forrest_genes, expressed_genes_fpkm_0)
forrest_expressed_0_5 <- intersect(forrest_genes, expressed_genes_fpkm_0_5)
forrest_expressed_1 <- intersect(forrest_genes, expressed_genes_fpkm_1)

tcf4_gene_expression_summary <- data.frame(
  Dataset = c("NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"),
  Total_TCF4_Genes = c(
    length(npc_genes),
    length(mcclay_genes),
    length(forrest_genes)
  ),
  Expressed_Genes_FPKM_gt_0 = c(
    length(npc_expressed_0),
    length(mcclay_expressed_0),
    length(forrest_expressed_0)
  ),
  Expressed_Genes_FPKM_gt_0_5 = c(
    length(npc_expressed_0_5),
    length(mcclay_expressed_0_5),
    length(forrest_expressed_0_5)
  ),
  Expressed_Genes_FPKM_gt_1 = c(
    length(npc_expressed_1),
    length(mcclay_expressed_1),
    length(forrest_expressed_1)
  )
)

tcf4_gene_expression_summary$Percent_Expressed_FPKM_gt_1 <- round(
  tcf4_gene_expression_summary$Expressed_Genes_FPKM_gt_1 /
    tcf4_gene_expression_summary$Total_TCF4_Genes * 100,
  3
)

tcf4_gene_expression_summary

# Add TCF4 target labels to RNA-seq table
rnaseq_df$NPC_TCF4_Target <- rnaseq_df$Gene_Symbol %in% npc_genes
rnaseq_df$McClay_TCF4_Target <- rnaseq_df$Gene_Symbol %in% mcclay_genes
rnaseq_df$Forrest_TCF4_Target <- rnaseq_df$Gene_Symbol %in% forrest_genes

rnaseq_df$log2_FPKM_plus1 <- log2(rnaseq_df$FPKM + 1)

# Compare expression of TCF4 target genes vs non-target genes
npc_expression_test <- wilcox.test(
  log2_FPKM_plus1 ~ NPC_TCF4_Target,
  data = rnaseq_df
)

mcclay_expression_test <- wilcox.test(
  log2_FPKM_plus1 ~ McClay_TCF4_Target,
  data = rnaseq_df
)

forrest_expression_test <- wilcox.test(
  log2_FPKM_plus1 ~ Forrest_TCF4_Target,
  data = rnaseq_df
)

npc_expression_test
mcclay_expression_test
forrest_expression_test

expression_test_summary <- data.frame(
  Dataset = c("NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"),
  Median_Expression_Target_Genes = c(
    median(rnaseq_df$log2_FPKM_plus1[rnaseq_df$NPC_TCF4_Target], na.rm = TRUE),
    median(rnaseq_df$log2_FPKM_plus1[rnaseq_df$McClay_TCF4_Target], na.rm = TRUE),
    median(rnaseq_df$log2_FPKM_plus1[rnaseq_df$Forrest_TCF4_Target], na.rm = TRUE)
  ),
  Median_Expression_NonTarget_Genes = c(
    median(rnaseq_df$log2_FPKM_plus1[!rnaseq_df$NPC_TCF4_Target], na.rm = TRUE),
    median(rnaseq_df$log2_FPKM_plus1[!rnaseq_df$McClay_TCF4_Target], na.rm = TRUE),
    median(rnaseq_df$log2_FPKM_plus1[!rnaseq_df$Forrest_TCF4_Target], na.rm = TRUE)
  ),
  Wilcoxon_P_Value = c(
    npc_expression_test$p.value,
    mcclay_expression_test$p.value,
    forrest_expression_test$p.value
  )
)

expression_test_summary

# Optional boxplots
boxplot(
  log2_FPKM_plus1 ~ NPC_TCF4_Target,
  data = rnaseq_df,
  main = "Expression of NPC TCF4 Target Genes",
  xlab = "NPC TCF4 Target",
  ylab = "log2(FPKM + 1)"
)

boxplot(
  log2_FPKM_plus1 ~ McClay_TCF4_Target,
  data = rnaseq_df,
  main = "Expression of McClay TCF4 Target Genes",
  xlab = "McClay TCF4 Target",
  ylab = "log2(FPKM + 1)"
)

boxplot(
  log2_FPKM_plus1 ~ Forrest_TCF4_Target,
  data = rnaseq_df,
  main = "Expression of Forrest TCF4 Target Genes",
  xlab = "Forrest TCF4 Target",
  ylab = "log2(FPKM + 1)"
)

# RNA-seq integration with TADs
hg19_genes_gr <- genes(txdb)

seqlevelsStyle(hg19_genes_gr) <- "UCSC"
hg19_genes_gr <- keepStandardChromosomes(
  hg19_genes_gr,
  pruning.mode = "coarse"
)

hg19_gene_entrez <- names(hg19_genes_gr)

hg19_gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = hg19_gene_entrez,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

mcols(hg19_genes_gr)$ENTREZID <- hg19_gene_entrez
mcols(hg19_genes_gr)$SYMBOL <- hg19_gene_symbols

hg19_genes_gr <- hg19_genes_gr[
  !is.na(mcols(hg19_genes_gr)$SYMBOL)
]

# Map genes to TADs
gene_tad_hits <- findOverlaps(hg19_genes_gr, hic_tads)

gene_tad_table <- data.frame(
  Gene_Symbol = mcols(hg19_genes_gr)$SYMBOL[queryHits(gene_tad_hits)],
  TAD_ID = tad_summary$TAD_ID[subjectHits(gene_tad_hits)]
)

gene_tad_table <- unique(gene_tad_table)

# Add RNA-seq expression to gene-TAD table
gene_tad_expression <- gene_tad_table %>%
  left_join(rnaseq_df, by = "Gene_Symbol")

head(gene_tad_expression)

# Calculate expression summary per TAD
tad_expression_summary <- gene_tad_expression %>%
  group_by(TAD_ID) %>%
  summarise(
    Genes_in_TAD = n_distinct(Gene_Symbol),
    Genes_with_RNAseq_Data = sum(!is.na(FPKM)),
    Expressed_Genes_FPKM_gt_1 = sum(FPKM > 1, na.rm = TRUE),
    Mean_FPKM = mean(FPKM, na.rm = TRUE),
    Median_FPKM = median(FPKM, na.rm = TRUE),
    Mean_log2_FPKM_plus1 = mean(log2_FPKM_plus1, na.rm = TRUE),
    Median_log2_FPKM_plus1 = median(log2_FPKM_plus1, na.rm = TRUE)
  )

# Add expression summary to TAD summary
tad_summary_expression <- tad_summary %>%
  left_join(tad_expression_summary, by = "TAD_ID")

head(tad_summary_expression)

# Correlate TAD regulatory features with expression
tad_expression_correlation <- cor(
  tad_summary_expression[, c(
    "NPC_TCF4_Peaks",
    "McClay_TCF4_Peaks",
    "Forrest_TCF4_Peaks",
    "H3K27ac_Peaks",
    "ATAC_Peaks",
    "Neuron_CHD8_Peaks",
    "Brain_CHD8_Peaks",
    "NPC_Regulatory_Score",
    "McClay_Regulatory_Score",
    "Forrest_Regulatory_Score",
    "Mean_log2_FPKM_plus1",
    "Median_log2_FPKM_plus1",
    "Expressed_Genes_FPKM_gt_1"
  )],
  use = "pairwise.complete.obs",
  method = "spearman"
)

round(tad_expression_correlation, 3)

# Extract most relevant correlations
tad_expression_correlation_summary <- data.frame(
  Feature = c(
    "NPC_TCF4_Peaks",
    "McClay_TCF4_Peaks",
    "Forrest_TCF4_Peaks",
    "H3K27ac_Peaks",
    "ATAC_Peaks",
    "Neuron_CHD8_Peaks",
    "Brain_CHD8_Peaks",
    "NPC_Regulatory_Score",
    "McClay_Regulatory_Score",
    "Forrest_Regulatory_Score"
  ),
  Correlation_with_Mean_log2_FPKM = c(
    tad_expression_correlation["NPC_TCF4_Peaks", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["McClay_TCF4_Peaks", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["Forrest_TCF4_Peaks", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["H3K27ac_Peaks", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["ATAC_Peaks", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["Neuron_CHD8_Peaks", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["Brain_CHD8_Peaks", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["NPC_Regulatory_Score", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["McClay_Regulatory_Score", "Mean_log2_FPKM_plus1"],
    tad_expression_correlation["Forrest_Regulatory_Score", "Mean_log2_FPKM_plus1"]
  ),
  Correlation_with_Expressed_Genes_Count = c(
    tad_expression_correlation["NPC_TCF4_Peaks", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["McClay_TCF4_Peaks", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["Forrest_TCF4_Peaks", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["H3K27ac_Peaks", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["ATAC_Peaks", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["Neuron_CHD8_Peaks", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["Brain_CHD8_Peaks", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["NPC_Regulatory_Score", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["McClay_Regulatory_Score", "Expressed_Genes_FPKM_gt_1"],
    tad_expression_correlation["Forrest_Regulatory_Score", "Expressed_Genes_FPKM_gt_1"]
  )
)

tad_expression_correlation_summary[, 2:3] <- round(
  tad_expression_correlation_summary[, 2:3],
  3
)

tad_expression_correlation_summary

# Analyze TAD_1375 expression
tad_1375_expression <- tad_1375_gene_table %>%
  left_join(rnaseq_df, by = "Gene_Symbol") %>%
  arrange(desc(FPKM))

tad_1375_expression

tad_1375_expression_summary <- data.frame(
  TAD_ID = "TAD_1375",
  Region = "chr19:57810001-58980000",
  Total_Genes = nrow(tad_1375_expression),
  Genes_with_RNAseq_Data = sum(!is.na(tad_1375_expression$FPKM)),
  Expressed_Genes_FPKM_gt_0 = sum(tad_1375_expression$FPKM > 0, na.rm = TRUE),
  Expressed_Genes_FPKM_gt_0_5 = sum(tad_1375_expression$FPKM > 0.5, na.rm = TRUE),
  Expressed_Genes_FPKM_gt_1 = sum(tad_1375_expression$FPKM > 1, na.rm = TRUE),
  Mean_FPKM = mean(tad_1375_expression$FPKM, na.rm = TRUE),
  Median_FPKM = median(tad_1375_expression$FPKM, na.rm = TRUE)
)

tad_1375_expression_summary

# Top expressed genes in TAD_1375
top_expressed_tad_1375_genes <- tad_1375_expression %>%
  filter(!is.na(FPKM)) %>%
  arrange(desc(FPKM))

head(top_expressed_tad_1375_genes, 20)

# Compare expression of genes inside NPC active CHD8 TADs
npc_active_chd8_tad_ids <- tad_summary$TAD_ID[
  tad_summary$NPC_Active_CHD8_TAD == TRUE
]

npc_active_chd8_tad_genes <- gene_tad_table$Gene_Symbol[
  gene_tad_table$TAD_ID %in% npc_active_chd8_tad_ids
]

npc_active_chd8_tad_genes <- unique(npc_active_chd8_tad_genes)

npc_active_chd8_tad_expression <- rnaseq_df %>%
  filter(Gene_Symbol %in% npc_active_chd8_tad_genes) %>%
  arrange(desc(FPKM))

npc_active_chd8_tad_expression_summary <- data.frame(
  Gene_Set = "Genes in NPC TCF4 + H3K27ac + ATAC + CHD8 TADs",
  Total_Genes = length(npc_active_chd8_tad_genes),
  Genes_with_RNAseq_Data = nrow(npc_active_chd8_tad_expression),
  Expressed_Genes_FPKM_gt_1 = sum(npc_active_chd8_tad_expression$FPKM > 1, na.rm = TRUE),
  Median_FPKM = median(npc_active_chd8_tad_expression$FPKM, na.rm = TRUE),
  Mean_FPKM = mean(npc_active_chd8_tad_expression$FPKM, na.rm = TRUE)
)

npc_active_chd8_tad_expression_summary$Percent_Expressed_FPKM_gt_1 <- round(
  npc_active_chd8_tad_expression_summary$Expressed_Genes_FPKM_gt_1 /
    npc_active_chd8_tad_expression_summary$Genes_with_RNAseq_Data * 100,
  3
)

npc_active_chd8_tad_expression_summary
