############################################################
# Peak-level and gene-level overlap analysis
# Genome builds used:
#   - Most peak overlaps: hg38 / GRCh38
#   - GSE79965 TAD/RNA-seq integration: hg19 / GRCh37
############################################################

############################################################
# 0. Libraries and shared setup
############################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(readxl)
  library(dplyr)
  library(R.utils)
})

txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

standard_chr <- paste0("chr", c(1:22, "X", "Y"))

############################################################
# 1. Core TCF4 and CTCF overlap analysis, hg38
############################################################

# Inputs:
#   McClay_TCF4_11322_consensus_hg38_sorted.bed
#   NPC_hg38_summit_500bp.bed
#   Forrest_hg38.bed
#   GSE123202_ENCFF102XIH_conservative_idr_thresholded_peaks_GRCh38.bed

mcclay_peaks <- import("McClay_TCF4_11322_consensus_hg38_sorted.bed")
npc_peaks <- import("NPC_hg38_summit_500bp.bed")
forrest_peaks <- import("Forrest_hg38.bed")

ctcf_df <- read.table(
  "GSE123202_ENCFF102XIH_conservative_idr_thresholded_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

ctcf_df <- ctcf_df[, 1:3]
colnames(ctcf_df) <- c("chr", "start", "end")
ctcf_df$start <- as.numeric(ctcf_df$start)
ctcf_df$end <- as.numeric(ctcf_df$end)
ctcf_df <- ctcf_df[!is.na(ctcf_df$start) & !is.na(ctcf_df$end) & ctcf_df$end > ctcf_df$start, ]

ctcf_peaks <- GRanges(
  seqnames = ctcf_df$chr,
  ranges = IRanges(start = ctcf_df$start + 1, end = ctcf_df$end)
)

seqlevelsStyle(mcclay_peaks) <- "UCSC"
seqlevelsStyle(npc_peaks) <- "UCSC"
seqlevelsStyle(forrest_peaks) <- "UCSC"
seqlevelsStyle(ctcf_peaks) <- "UCSC"

mcclay_peaks <- keepSeqlevels(mcclay_peaks, intersect(seqlevels(mcclay_peaks), standard_chr), pruning.mode = "coarse")
npc_peaks <- keepSeqlevels(npc_peaks, intersect(seqlevels(npc_peaks), standard_chr), pruning.mode = "coarse")
forrest_peaks <- keepSeqlevels(forrest_peaks, intersect(seqlevels(forrest_peaks), standard_chr), pruning.mode = "coarse")
ctcf_peaks <- keepSeqlevels(ctcf_peaks, intersect(seqlevels(ctcf_peaks), standard_chr), pruning.mode = "coarse")

mcclay_peaks <- sort(mcclay_peaks)
npc_peaks <- sort(npc_peaks)
forrest_peaks <- sort(forrest_peaks)
ctcf_peaks <- sort(ctcf_peaks)

# Peak-level overlaps.
mcclay_ctcf_peaks <- subsetByOverlaps(mcclay_peaks, ctcf_peaks)
npc_ctcf_peaks <- subsetByOverlaps(npc_peaks, ctcf_peaks)
mcclay_npc_peaks <- subsetByOverlaps(mcclay_peaks, npc_peaks)
forrest_ctcf_peaks <- subsetByOverlaps(forrest_peaks, ctcf_peaks)
forrest_npc_peaks <- subsetByOverlaps(forrest_peaks, npc_peaks)
forrest_mcclay_peaks <- subsetByOverlaps(forrest_peaks, mcclay_peaks)
mcclay_npc_ctcf_peaks <- subsetByOverlaps(mcclay_npc_peaks, ctcf_peaks)
forrest_mcclay_npc_ctcf_peaks <- subsetByOverlaps(forrest_peaks, mcclay_npc_ctcf_peaks)

peak_overlap_summary <- data.frame(
  Comparison = c(
    "McClay_TCF4_vs_CTCF",
    "NPC_TCF4_vs_CTCF",
    "McClay_TCF4_vs_NPC_TCF4",
    "McClay_TCF4_vs_NPC_TCF4_vs_CTCF",
    "Forrest_TCF4_vs_CTCF",
    "Forrest_TCF4_vs_McClay_TCF4",
    "Forrest_TCF4_vs_NPC_TCF4",
    "Forrest_TCF4_vs_McClay_TCF4_vs_NPC_TCF4_vs_CTCF"
  ),
  Peak_Count = c(
    length(mcclay_ctcf_peaks),
    length(npc_ctcf_peaks),
    length(mcclay_npc_peaks),
    length(mcclay_npc_ctcf_peaks),
    length(forrest_ctcf_peaks),
    length(forrest_mcclay_peaks),
    length(forrest_npc_peaks),
    length(forrest_mcclay_npc_ctcf_peaks)
  ),
  First_Dataset_Size = c(
    length(mcclay_peaks),
    length(npc_peaks),
    length(mcclay_peaks),
    length(mcclay_peaks),
    length(forrest_peaks),
    length(forrest_peaks),
    length(forrest_peaks),
    length(forrest_peaks)
  )
)

peak_overlap_summary$Percent_of_First_Dataset <- round(
  peak_overlap_summary$Peak_Count / peak_overlap_summary$First_Dataset_Size * 100,
  3
)

peak_overlap_summary
write.csv(peak_overlap_summary, "peak_overlap_summary_tcf4_ctcf_hg38.csv", row.names = FALSE)

# Gene annotation for complete and overlap peak sets.
mcclay_anno_df <- as.data.frame(annotatePeak(mcclay_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))
npc_anno_df <- as.data.frame(annotatePeak(npc_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))
forrest_anno_df <- as.data.frame(annotatePeak(forrest_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))
ctcf_anno_df <- as.data.frame(annotatePeak(ctcf_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))
mcclay_ctcf_anno_df <- as.data.frame(annotatePeak(mcclay_ctcf_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))
npc_ctcf_anno_df <- as.data.frame(annotatePeak(npc_ctcf_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))
mcclay_npc_ctcf_anno_df <- as.data.frame(annotatePeak(mcclay_npc_ctcf_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))

mcclay_genes <- unique(na.omit(mcclay_anno_df$SYMBOL))
npc_genes <- unique(na.omit(npc_anno_df$SYMBOL))
forrest_genes <- unique(na.omit(forrest_anno_df$SYMBOL))
ctcf_genes <- unique(na.omit(ctcf_anno_df$SYMBOL))
mcclay_ctcf_genes <- unique(na.omit(mcclay_ctcf_anno_df$SYMBOL))
npc_ctcf_genes <- unique(na.omit(npc_ctcf_anno_df$SYMBOL))
mcclay_npc_ctcf_genes <- unique(na.omit(mcclay_npc_ctcf_anno_df$SYMBOL))

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
    length(Reduce(intersect, list(mcclay_genes, npc_genes, ctcf_genes))),
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
  Second_Gene_Set_Size = c(
    length(ctcf_genes),
    length(ctcf_genes),
    length(npc_genes),
    length(ctcf_genes),
    length(ctcf_genes),
    length(npc_genes),
    length(mcclay_genes)
  )
)

gene_overlap_summary$Percent_of_First_Gene_Set <- round(
  gene_overlap_summary$Overlap_Genes / gene_overlap_summary$First_Gene_Set_Size * 100,
  2
)

gene_overlap_summary$Percent_of_Second_Gene_Set <- round(
  gene_overlap_summary$Overlap_Genes / gene_overlap_summary$Second_Gene_Set_Size * 100,
  2
)

gene_overlap_summary
write.csv(gene_overlap_summary, "gene_overlap_summary_tcf4_ctcf_hg38.csv", row.names = FALSE)

############################################################
# 2. FGFR1 and NOTCH overlap analysis, hg19 lifted to hg38
############################################################

# Inputs:
#   FGFR1_NOTCH/GSM2439179_C4-R1.bed.gz
#   FGFR1_NOTCH/GSM2439181_P3-R1.bed.gz
#   FGFR1_NOTCH/GSM2439182_P3-N.bed.gz

control_fgfr1_hg19 <- import("FGFR1_NOTCH/GSM2439179_C4-R1.bed.gz")
patient_fgfr1_hg19 <- import("FGFR1_NOTCH/GSM2439181_P3-R1.bed.gz")
patient_notch_hg19 <- import("FGFR1_NOTCH/GSM2439182_P3-N.bed.gz")

download.file(
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
  destfile = "hg19ToHg38.over.chain.gz",
  mode = "wb"
)

R.utils::gunzip("hg19ToHg38.over.chain.gz", remove = FALSE, overwrite = TRUE)
hg19_to_hg38_chain <- import.chain("hg19ToHg38.over.chain")

control_fgfr1_hg38_list <- liftOver(control_fgfr1_hg19, hg19_to_hg38_chain)
patient_fgfr1_hg38_list <- liftOver(patient_fgfr1_hg19, hg19_to_hg38_chain)
patient_notch_hg38_list <- liftOver(patient_notch_hg19, hg19_to_hg38_chain)

control_fgfr1_hg38 <- unlist(control_fgfr1_hg38_list[elementNROWS(control_fgfr1_hg38_list) == 1])
patient_fgfr1_hg38 <- unlist(patient_fgfr1_hg38_list[elementNROWS(patient_fgfr1_hg38_list) == 1])
patient_notch_hg38 <- unlist(patient_notch_hg38_list[elementNROWS(patient_notch_hg38_list) == 1])

seqlevelsStyle(control_fgfr1_hg38) <- "UCSC"
seqlevelsStyle(patient_fgfr1_hg38) <- "UCSC"
seqlevelsStyle(patient_notch_hg38) <- "UCSC"

control_fgfr1_hg38 <- sort(keepStandardChromosomes(control_fgfr1_hg38, pruning.mode = "coarse"))
patient_fgfr1_hg38 <- sort(keepStandardChromosomes(patient_fgfr1_hg38, pruning.mode = "coarse"))
patient_notch_hg38 <- sort(keepStandardChromosomes(patient_notch_hg38, pruning.mode = "coarse"))

export(control_fgfr1_hg38, "FGFR1_NOTCH/Control_FGFR1_hg38_clean.bed")
export(patient_fgfr1_hg38, "FGFR1_NOTCH/Patient_FGFR1_hg38_clean.bed")
export(patient_notch_hg38, "FGFR1_NOTCH/Patient_NOTCH_hg38_clean.bed")

# Peak-level pairwise overlaps against McClay, NPC, and Forrest TCF4.
mcclay_control_fgfr1_hits <- findOverlaps(mcclay_peaks, control_fgfr1_hg38)
mcclay_patient_fgfr1_hits <- findOverlaps(mcclay_peaks, patient_fgfr1_hg38)
mcclay_patient_notch_hits <- findOverlaps(mcclay_peaks, patient_notch_hg38)

npc_control_fgfr1_hits <- findOverlaps(npc_peaks, control_fgfr1_hg38)
npc_patient_fgfr1_hits <- findOverlaps(npc_peaks, patient_fgfr1_hg38)
npc_patient_notch_hits <- findOverlaps(npc_peaks, patient_notch_hg38)

forrest_control_fgfr1_hits <- findOverlaps(forrest_peaks, control_fgfr1_hg38)
forrest_patient_fgfr1_hits <- findOverlaps(forrest_peaks, patient_fgfr1_hg38)
forrest_patient_notch_hits <- findOverlaps(forrest_peaks, patient_notch_hg38)

fgfr1_notch_peak_summary <- data.frame(
  Comparison = c(
    "McClay_TCF4_vs_Control_FGFR1",
    "McClay_TCF4_vs_Patient_FGFR1",
    "McClay_TCF4_vs_Patient_NOTCH",
    "NPC_TCF4_vs_Control_FGFR1",
    "NPC_TCF4_vs_Patient_FGFR1",
    "NPC_TCF4_vs_Patient_NOTCH",
    "Forrest_TCF4_vs_Control_FGFR1",
    "Forrest_TCF4_vs_Patient_FGFR1",
    "Forrest_TCF4_vs_Patient_NOTCH"
  ),
  Query_Peaks = c(
    length(mcclay_peaks), length(mcclay_peaks), length(mcclay_peaks),
    length(npc_peaks), length(npc_peaks), length(npc_peaks),
    length(forrest_peaks), length(forrest_peaks), length(forrest_peaks)
  ),
  Subject_Peaks = c(
    length(control_fgfr1_hg38), length(patient_fgfr1_hg38), length(patient_notch_hg38),
    length(control_fgfr1_hg38), length(patient_fgfr1_hg38), length(patient_notch_hg38),
    length(control_fgfr1_hg38), length(patient_fgfr1_hg38), length(patient_notch_hg38)
  ),
  Overlap_Query_Peaks = c(
    length(unique(queryHits(mcclay_control_fgfr1_hits))),
    length(unique(queryHits(mcclay_patient_fgfr1_hits))),
    length(unique(queryHits(mcclay_patient_notch_hits))),
    length(unique(queryHits(npc_control_fgfr1_hits))),
    length(unique(queryHits(npc_patient_fgfr1_hits))),
    length(unique(queryHits(npc_patient_notch_hits))),
    length(unique(queryHits(forrest_control_fgfr1_hits))),
    length(unique(queryHits(forrest_patient_fgfr1_hits))),
    length(unique(queryHits(forrest_patient_notch_hits)))
  ),
  Overlap_Subject_Peaks = c(
    length(unique(subjectHits(mcclay_control_fgfr1_hits))),
    length(unique(subjectHits(mcclay_patient_fgfr1_hits))),
    length(unique(subjectHits(mcclay_patient_notch_hits))),
    length(unique(subjectHits(npc_control_fgfr1_hits))),
    length(unique(subjectHits(npc_patient_fgfr1_hits))),
    length(unique(subjectHits(npc_patient_notch_hits))),
    length(unique(subjectHits(forrest_control_fgfr1_hits))),
    length(unique(subjectHits(forrest_patient_fgfr1_hits))),
    length(unique(subjectHits(forrest_patient_notch_hits)))
  )
)

fgfr1_notch_peak_summary$Percent_Query_Overlap <- round(
  fgfr1_notch_peak_summary$Overlap_Query_Peaks / fgfr1_notch_peak_summary$Query_Peaks * 100,
  3
)

fgfr1_notch_peak_summary$Percent_Subject_Overlap <- round(
  fgfr1_notch_peak_summary$Overlap_Subject_Peaks / fgfr1_notch_peak_summary$Subject_Peaks * 100,
  3
)

fgfr1_notch_peak_summary
write.csv(fgfr1_notch_peak_summary, "peak_overlap_summary_fgfr1_notch_hg38.csv", row.names = FALSE)

# Gene-level FGFR1/NOTCH overlaps.
control_fgfr1_anno_df <- as.data.frame(annotatePeak(control_fgfr1_hg38, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))
patient_fgfr1_anno_df <- as.data.frame(annotatePeak(patient_fgfr1_hg38, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))
patient_notch_anno_df <- as.data.frame(annotatePeak(patient_notch_hg38, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))

control_fgfr1_genes <- unique(na.omit(control_fgfr1_anno_df$SYMBOL))
patient_fgfr1_genes <- unique(na.omit(patient_fgfr1_anno_df$SYMBOL))
patient_notch_genes <- unique(na.omit(patient_notch_anno_df$SYMBOL))

fgfr1_notch_gene_summary <- data.frame(
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
    length(mcclay_genes), length(mcclay_genes), length(mcclay_genes),
    length(npc_genes), length(npc_genes), length(npc_genes),
    length(forrest_genes), length(forrest_genes), length(forrest_genes)
  ),
  Second_Gene_Set_Size = c(
    length(control_fgfr1_genes), length(patient_fgfr1_genes), length(patient_notch_genes),
    length(control_fgfr1_genes), length(patient_fgfr1_genes), length(patient_notch_genes),
    length(control_fgfr1_genes), length(patient_fgfr1_genes), length(patient_notch_genes)
  )
)

fgfr1_notch_gene_summary$Percent_of_First_Gene_Set <- round(
  fgfr1_notch_gene_summary$Overlap_Genes / fgfr1_notch_gene_summary$First_Gene_Set_Size * 100,
  2
)

fgfr1_notch_gene_summary$Percent_of_Second_Gene_Set <- round(
  fgfr1_notch_gene_summary$Overlap_Genes / fgfr1_notch_gene_summary$Second_Gene_Set_Size * 100,
  2
)

fgfr1_notch_gene_summary
write.csv(fgfr1_notch_gene_summary, "gene_overlap_summary_fgfr1_notch_hg38.csv", row.names = FALSE)

############################################################
# 3. Gene-level Fisher tests for FGFR1 and NOTCH overlaps
############################################################

background_gene <- read.csv("final_ref_gene_2026_filtered.csv", stringsAsFactors = FALSE)
all_genes <- unique(toupper(na.omit(background_gene$gene)))

mcclay_genes_bg <- intersect(unique(toupper(mcclay_genes)), all_genes)
npc_genes_bg <- intersect(unique(toupper(npc_genes)), all_genes)
forrest_genes_bg <- intersect(unique(toupper(forrest_genes)), all_genes)
control_fgfr1_genes_bg <- intersect(unique(toupper(control_fgfr1_genes)), all_genes)
patient_fgfr1_genes_bg <- intersect(unique(toupper(patient_fgfr1_genes)), all_genes)
patient_notch_genes_bg <- intersect(unique(toupper(patient_notch_genes)), all_genes)

mcclay_control_fgfr1_fisher <- fisher.test(matrix(c(
  length(intersect(mcclay_genes_bg, control_fgfr1_genes_bg)),
  length(setdiff(mcclay_genes_bg, control_fgfr1_genes_bg)),
  length(setdiff(control_fgfr1_genes_bg, mcclay_genes_bg)),
  length(setdiff(all_genes, union(mcclay_genes_bg, control_fgfr1_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

mcclay_patient_fgfr1_fisher <- fisher.test(matrix(c(
  length(intersect(mcclay_genes_bg, patient_fgfr1_genes_bg)),
  length(setdiff(mcclay_genes_bg, patient_fgfr1_genes_bg)),
  length(setdiff(patient_fgfr1_genes_bg, mcclay_genes_bg)),
  length(setdiff(all_genes, union(mcclay_genes_bg, patient_fgfr1_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

mcclay_patient_notch_fisher <- fisher.test(matrix(c(
  length(intersect(mcclay_genes_bg, patient_notch_genes_bg)),
  length(setdiff(mcclay_genes_bg, patient_notch_genes_bg)),
  length(setdiff(patient_notch_genes_bg, mcclay_genes_bg)),
  length(setdiff(all_genes, union(mcclay_genes_bg, patient_notch_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

npc_control_fgfr1_fisher <- fisher.test(matrix(c(
  length(intersect(npc_genes_bg, control_fgfr1_genes_bg)),
  length(setdiff(npc_genes_bg, control_fgfr1_genes_bg)),
  length(setdiff(control_fgfr1_genes_bg, npc_genes_bg)),
  length(setdiff(all_genes, union(npc_genes_bg, control_fgfr1_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

npc_patient_fgfr1_fisher <- fisher.test(matrix(c(
  length(intersect(npc_genes_bg, patient_fgfr1_genes_bg)),
  length(setdiff(npc_genes_bg, patient_fgfr1_genes_bg)),
  length(setdiff(patient_fgfr1_genes_bg, npc_genes_bg)),
  length(setdiff(all_genes, union(npc_genes_bg, patient_fgfr1_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

npc_patient_notch_fisher <- fisher.test(matrix(c(
  length(intersect(npc_genes_bg, patient_notch_genes_bg)),
  length(setdiff(npc_genes_bg, patient_notch_genes_bg)),
  length(setdiff(patient_notch_genes_bg, npc_genes_bg)),
  length(setdiff(all_genes, union(npc_genes_bg, patient_notch_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

forrest_control_fgfr1_fisher <- fisher.test(matrix(c(
  length(intersect(forrest_genes_bg, control_fgfr1_genes_bg)),
  length(setdiff(forrest_genes_bg, control_fgfr1_genes_bg)),
  length(setdiff(control_fgfr1_genes_bg, forrest_genes_bg)),
  length(setdiff(all_genes, union(forrest_genes_bg, control_fgfr1_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

forrest_patient_fgfr1_fisher <- fisher.test(matrix(c(
  length(intersect(forrest_genes_bg, patient_fgfr1_genes_bg)),
  length(setdiff(forrest_genes_bg, patient_fgfr1_genes_bg)),
  length(setdiff(patient_fgfr1_genes_bg, forrest_genes_bg)),
  length(setdiff(all_genes, union(forrest_genes_bg, patient_fgfr1_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

forrest_patient_notch_fisher <- fisher.test(matrix(c(
  length(intersect(forrest_genes_bg, patient_notch_genes_bg)),
  length(setdiff(forrest_genes_bg, patient_notch_genes_bg)),
  length(setdiff(patient_notch_genes_bg, forrest_genes_bg)),
  length(setdiff(all_genes, union(forrest_genes_bg, patient_notch_genes_bg)))
), nrow = 2, byrow = TRUE), alternative = "greater")

fgfr1_notch_fisher_summary <- data.frame(
  Comparison = fgfr1_notch_gene_summary$Comparison,
  Set1_Size = c(
    length(mcclay_genes_bg), length(mcclay_genes_bg), length(mcclay_genes_bg),
    length(npc_genes_bg), length(npc_genes_bg), length(npc_genes_bg),
    length(forrest_genes_bg), length(forrest_genes_bg), length(forrest_genes_bg)
  ),
  Set2_Size = c(
    length(control_fgfr1_genes_bg), length(patient_fgfr1_genes_bg), length(patient_notch_genes_bg),
    length(control_fgfr1_genes_bg), length(patient_fgfr1_genes_bg), length(patient_notch_genes_bg),
    length(control_fgfr1_genes_bg), length(patient_fgfr1_genes_bg), length(patient_notch_genes_bg)
  ),
  Overlap = c(
    length(intersect(mcclay_genes_bg, control_fgfr1_genes_bg)),
    length(intersect(mcclay_genes_bg, patient_fgfr1_genes_bg)),
    length(intersect(mcclay_genes_bg, patient_notch_genes_bg)),
    length(intersect(npc_genes_bg, control_fgfr1_genes_bg)),
    length(intersect(npc_genes_bg, patient_fgfr1_genes_bg)),
    length(intersect(npc_genes_bg, patient_notch_genes_bg)),
    length(intersect(forrest_genes_bg, control_fgfr1_genes_bg)),
    length(intersect(forrest_genes_bg, patient_fgfr1_genes_bg)),
    length(intersect(forrest_genes_bg, patient_notch_genes_bg))
  ),
  Odds_Ratio = c(
    unname(mcclay_control_fgfr1_fisher$estimate),
    unname(mcclay_patient_fgfr1_fisher$estimate),
    unname(mcclay_patient_notch_fisher$estimate),
    unname(npc_control_fgfr1_fisher$estimate),
    unname(npc_patient_fgfr1_fisher$estimate),
    unname(npc_patient_notch_fisher$estimate),
    unname(forrest_control_fgfr1_fisher$estimate),
    unname(forrest_patient_fgfr1_fisher$estimate),
    unname(forrest_patient_notch_fisher$estimate)
  ),
  P_value = c(
    mcclay_control_fgfr1_fisher$p.value,
    mcclay_patient_fgfr1_fisher$p.value,
    mcclay_patient_notch_fisher$p.value,
    npc_control_fgfr1_fisher$p.value,
    npc_patient_fgfr1_fisher$p.value,
    npc_patient_notch_fisher$p.value,
    forrest_control_fgfr1_fisher$p.value,
    forrest_patient_fgfr1_fisher$p.value,
    forrest_patient_notch_fisher$p.value
  )
)

fgfr1_notch_fisher_summary
write.csv(fgfr1_notch_fisher_summary, "fisher_gene_overlap_fgfr1_notch_hg38.csv", row.names = FALSE)

############################################################
# 4. NPC TCF4 with CTCF, FGFR1, and NOTCH
############################################################

npc_control_fgfr1_peaks <- subsetByOverlaps(npc_peaks, control_fgfr1_hg38)
npc_patient_fgfr1_peaks <- subsetByOverlaps(npc_peaks, patient_fgfr1_hg38)
npc_patient_notch_peaks <- subsetByOverlaps(npc_peaks, patient_notch_hg38)
npc_ctcf_control_fgfr1_peaks <- subsetByOverlaps(npc_ctcf_peaks, control_fgfr1_hg38)
npc_ctcf_patient_fgfr1_peaks <- subsetByOverlaps(npc_ctcf_peaks, patient_fgfr1_hg38)
npc_ctcf_patient_notch_peaks <- subsetByOverlaps(npc_ctcf_peaks, patient_notch_hg38)
npc_ctcf_patient_fgfr1_notch_peaks <- subsetByOverlaps(npc_ctcf_patient_fgfr1_peaks, patient_notch_hg38)

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
    length(npc_ctcf_peaks),
    length(npc_control_fgfr1_peaks),
    length(npc_patient_fgfr1_peaks),
    length(npc_patient_notch_peaks),
    length(npc_ctcf_control_fgfr1_peaks),
    length(npc_ctcf_patient_fgfr1_peaks),
    length(npc_ctcf_patient_notch_peaks),
    length(npc_ctcf_patient_fgfr1_notch_peaks)
  ) / length(npc_peaks) * 100, 2),
  Percent_of_NPC_CTCF = round(c(
    NA, NA, NA, NA,
    length(npc_ctcf_control_fgfr1_peaks),
    length(npc_ctcf_patient_fgfr1_peaks),
    length(npc_ctcf_patient_notch_peaks),
    length(npc_ctcf_patient_fgfr1_notch_peaks)
  ) / length(npc_ctcf_peaks) * 100, 2)
)

npc_ctcf_fgfr1_notch_peak_summary
write.csv(npc_ctcf_fgfr1_notch_peak_summary, "npc_ctcf_fgfr1_notch_peak_summary.csv", row.names = FALSE)

npc_ctcf_shared_genes <- intersect(npc_genes, ctcf_genes)
npc_control_fgfr1_shared_genes <- intersect(npc_genes, control_fgfr1_genes)
npc_patient_fgfr1_shared_genes <- intersect(npc_genes, patient_fgfr1_genes)
npc_patient_notch_shared_genes <- intersect(npc_genes, patient_notch_genes)
npc_ctcf_control_fgfr1_shared_genes <- Reduce(intersect, list(npc_genes, ctcf_genes, control_fgfr1_genes))
npc_ctcf_patient_fgfr1_shared_genes <- Reduce(intersect, list(npc_genes, ctcf_genes, patient_fgfr1_genes))
npc_ctcf_patient_notch_shared_genes <- Reduce(intersect, list(npc_genes, ctcf_genes, patient_notch_genes))
npc_ctcf_patient_fgfr1_notch_shared_genes <- Reduce(intersect, list(npc_genes, ctcf_genes, patient_fgfr1_genes, patient_notch_genes))

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
    length(npc_ctcf_shared_genes),
    length(npc_control_fgfr1_shared_genes),
    length(npc_patient_fgfr1_shared_genes),
    length(npc_patient_notch_shared_genes),
    length(npc_ctcf_control_fgfr1_shared_genes),
    length(npc_ctcf_patient_fgfr1_shared_genes),
    length(npc_ctcf_patient_notch_shared_genes),
    length(npc_ctcf_patient_fgfr1_notch_shared_genes)
  ) / length(npc_genes) * 100, 2),
  Percent_of_NPC_CTCF_Genes = round(c(
    NA, NA, NA, NA,
    length(npc_ctcf_control_fgfr1_shared_genes),
    length(npc_ctcf_patient_fgfr1_shared_genes),
    length(npc_ctcf_patient_notch_shared_genes),
    length(npc_ctcf_patient_fgfr1_notch_shared_genes)
  ) / length(npc_ctcf_shared_genes) * 100, 2)
)

npc_ctcf_fgfr1_notch_gene_summary
write.csv(npc_ctcf_fgfr1_notch_gene_summary, "npc_ctcf_fgfr1_notch_gene_summary.csv", row.names = FALSE)

############################################################
# 5. H3K27ac, H3K4me3, and H3K27me3 peak overlaps, hg38
############################################################

background_genes <- unique(na.omit(read.csv("final_ref_gene_2026_filtered.csv")$gene))
background_genes <- background_genes[background_genes != ""]

# H3K27ac.
h3k27ac_df <- read.table(
  "GSE96178_ENCFF407DWP_replicated_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

h3k27ac_df <- h3k27ac_df[, 1:3]
colnames(h3k27ac_df) <- c("chr", "start", "end")
h3k27ac_df$start <- as.numeric(h3k27ac_df$start)
h3k27ac_df$end <- as.numeric(h3k27ac_df$end)
h3k27ac_df <- h3k27ac_df[!is.na(h3k27ac_df$start) & !is.na(h3k27ac_df$end) & h3k27ac_df$end > h3k27ac_df$start, ]
h3k27ac_peaks <- GRanges(seqnames = h3k27ac_df$chr, ranges = IRanges(start = h3k27ac_df$start + 1, end = h3k27ac_df$end))
seqlevelsStyle(h3k27ac_peaks) <- "UCSC"
h3k27ac_peaks <- sort(keepStandardChromosomes(h3k27ac_peaks, pruning.mode = "coarse"))

# H3K4me3.
h3k4me3_df <- read.table(
  "GSE176914_ENCFF907ZRF_replicated_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

h3k4me3_df <- h3k4me3_df[, 1:3]
colnames(h3k4me3_df) <- c("chr", "start", "end")
h3k4me3_df$start <- as.numeric(h3k4me3_df$start)
h3k4me3_df$end <- as.numeric(h3k4me3_df$end)
h3k4me3_df <- h3k4me3_df[!is.na(h3k4me3_df$start) & !is.na(h3k4me3_df$end) & h3k4me3_df$end > h3k4me3_df$start, ]
h3k4me3_peaks <- GRanges(seqnames = h3k4me3_df$chr, ranges = IRanges(start = h3k4me3_df$start + 1, end = h3k4me3_df$end))
seqlevelsStyle(h3k4me3_peaks) <- "UCSC"
h3k4me3_peaks <- sort(keepStandardChromosomes(h3k4me3_peaks, pruning.mode = "coarse"))

# H3K27me3.
h3k27me3_df <- read.table(
  gzfile("ENCFF056AFA.bed.gz"),
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

h3k27me3_df <- h3k27me3_df[, 1:3]
colnames(h3k27me3_df) <- c("chr", "start", "end")
h3k27me3_df$start <- as.numeric(h3k27me3_df$start)
h3k27me3_df$end <- as.numeric(h3k27me3_df$end)
h3k27me3_df <- h3k27me3_df[!is.na(h3k27me3_df$start) & !is.na(h3k27me3_df$end) & h3k27me3_df$end > h3k27me3_df$start, ]
h3k27me3_peaks <- GRanges(seqnames = h3k27me3_df$chr, ranges = IRanges(start = h3k27me3_df$start + 1, end = h3k27me3_df$end))
seqlevelsStyle(h3k27me3_peaks) <- "UCSC"
h3k27me3_peaks <- sort(keepStandardChromosomes(h3k27me3_peaks, pruning.mode = "coarse"))

# Peak-level histone overlaps.
npc_h3k27ac_hits <- findOverlaps(npc_peaks, h3k27ac_peaks)
mcclay_h3k27ac_hits <- findOverlaps(mcclay_peaks, h3k27ac_peaks)
forrest_h3k27ac_hits <- findOverlaps(forrest_peaks, h3k27ac_peaks)

npc_h3k4me3_hits <- findOverlaps(npc_peaks, h3k4me3_peaks)
mcclay_h3k4me3_hits <- findOverlaps(mcclay_peaks, h3k4me3_peaks)
forrest_h3k4me3_hits <- findOverlaps(forrest_peaks, h3k4me3_peaks)

npc_h3k27me3_hits <- findOverlaps(npc_peaks, h3k27me3_peaks)
mcclay_h3k27me3_hits <- findOverlaps(mcclay_peaks, h3k27me3_peaks)
forrest_h3k27me3_hits <- findOverlaps(forrest_peaks, h3k27me3_peaks)

histone_peak_overlap_summary <- data.frame(
  Mark = c(
    "H3K27ac", "H3K27ac", "H3K27ac",
    "H3K4me3", "H3K4me3", "H3K4me3",
    "H3K27me3", "H3K27me3", "H3K27me3"
  ),
  TCF4_Dataset = c(
    "NPC_TCF4", "McClay_TCF4", "Forrest_TCF4",
    "NPC_TCF4", "McClay_TCF4", "Forrest_TCF4",
    "NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"
  ),
  TCF4_Total_Peaks = c(
    length(npc_peaks), length(mcclay_peaks), length(forrest_peaks),
    length(npc_peaks), length(mcclay_peaks), length(forrest_peaks),
    length(npc_peaks), length(mcclay_peaks), length(forrest_peaks)
  ),
  Mark_Total_Peaks = c(
    length(h3k27ac_peaks), length(h3k27ac_peaks), length(h3k27ac_peaks),
    length(h3k4me3_peaks), length(h3k4me3_peaks), length(h3k4me3_peaks),
    length(h3k27me3_peaks), length(h3k27me3_peaks), length(h3k27me3_peaks)
  ),
  TCF4_Overlapping_Peaks = c(
    length(unique(queryHits(npc_h3k27ac_hits))),
    length(unique(queryHits(mcclay_h3k27ac_hits))),
    length(unique(queryHits(forrest_h3k27ac_hits))),
    length(unique(queryHits(npc_h3k4me3_hits))),
    length(unique(queryHits(mcclay_h3k4me3_hits))),
    length(unique(queryHits(forrest_h3k4me3_hits))),
    length(unique(queryHits(npc_h3k27me3_hits))),
    length(unique(queryHits(mcclay_h3k27me3_hits))),
    length(unique(queryHits(forrest_h3k27me3_hits)))
  ),
  Mark_Overlapping_Peaks = c(
    length(unique(subjectHits(npc_h3k27ac_hits))),
    length(unique(subjectHits(mcclay_h3k27ac_hits))),
    length(unique(subjectHits(forrest_h3k27ac_hits))),
    length(unique(subjectHits(npc_h3k4me3_hits))),
    length(unique(subjectHits(mcclay_h3k4me3_hits))),
    length(unique(subjectHits(forrest_h3k4me3_hits))),
    length(unique(subjectHits(npc_h3k27me3_hits))),
    length(unique(subjectHits(mcclay_h3k27me3_hits))),
    length(unique(subjectHits(forrest_h3k27me3_hits)))
  )
)

histone_peak_overlap_summary$Percent_TCF4_Peaks_Overlapping_Mark <- round(
  histone_peak_overlap_summary$TCF4_Overlapping_Peaks / histone_peak_overlap_summary$TCF4_Total_Peaks * 100,
  3
)

histone_peak_overlap_summary$Percent_Mark_Peaks_Overlapping_TCF4 <- round(
  histone_peak_overlap_summary$Mark_Overlapping_Peaks / histone_peak_overlap_summary$Mark_Total_Peaks * 100,
  3
)

histone_peak_overlap_summary
write.csv(histone_peak_overlap_summary, "histone_peak_overlap_summary_hg38.csv", row.names = FALSE)

# Gene-level histone-overlapping TCF4 targets.
npc_h3k27ac_peaks <- npc_peaks[unique(queryHits(npc_h3k27ac_hits))]
mcclay_h3k27ac_peaks <- mcclay_peaks[unique(queryHits(mcclay_h3k27ac_hits))]
forrest_h3k27ac_peaks <- forrest_peaks[unique(queryHits(forrest_h3k27ac_hits))]

npc_h3k4me3_peaks <- npc_peaks[unique(queryHits(npc_h3k4me3_hits))]
mcclay_h3k4me3_peaks <- mcclay_peaks[unique(queryHits(mcclay_h3k4me3_hits))]
forrest_h3k4me3_peaks <- forrest_peaks[unique(queryHits(forrest_h3k4me3_hits))]

npc_h3k27me3_peaks <- npc_peaks[unique(queryHits(npc_h3k27me3_hits))]
mcclay_h3k27me3_peaks <- mcclay_peaks[unique(queryHits(mcclay_h3k27me3_hits))]
forrest_h3k27me3_peaks <- forrest_peaks[unique(queryHits(forrest_h3k27me3_hits))]

npc_h3k27ac_genes <- unique(na.omit(as.data.frame(annotatePeak(npc_h3k27ac_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
mcclay_h3k27ac_genes <- unique(na.omit(as.data.frame(annotatePeak(mcclay_h3k27ac_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
forrest_h3k27ac_genes <- unique(na.omit(as.data.frame(annotatePeak(forrest_h3k27ac_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))

npc_h3k4me3_genes <- unique(na.omit(as.data.frame(annotatePeak(npc_h3k4me3_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
mcclay_h3k4me3_genes <- unique(na.omit(as.data.frame(annotatePeak(mcclay_h3k4me3_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
forrest_h3k4me3_genes <- unique(na.omit(as.data.frame(annotatePeak(forrest_h3k4me3_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))

npc_h3k27me3_genes <- unique(na.omit(as.data.frame(annotatePeak(npc_h3k27me3_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
mcclay_h3k27me3_genes <- unique(na.omit(as.data.frame(annotatePeak(mcclay_h3k27me3_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
forrest_h3k27me3_genes <- unique(na.omit(as.data.frame(annotatePeak(forrest_h3k27me3_peaks, TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))

histone_gene_summary <- data.frame(
  Mark = c(
    "H3K27ac", "H3K27ac", "H3K27ac",
    "H3K4me3", "H3K4me3", "H3K4me3",
    "H3K27me3", "H3K27me3", "H3K27me3"
  ),
  TCF4_Dataset = c(
    "NPC_TCF4", "McClay_TCF4", "Forrest_TCF4",
    "NPC_TCF4", "McClay_TCF4", "Forrest_TCF4",
    "NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"
  ),
  Total_TCF4_Genes = c(
    length(npc_genes), length(mcclay_genes), length(forrest_genes),
    length(npc_genes), length(mcclay_genes), length(forrest_genes),
    length(npc_genes), length(mcclay_genes), length(forrest_genes)
  ),
  Mark_Overlapping_TCF4_Genes = c(
    length(npc_h3k27ac_genes), length(mcclay_h3k27ac_genes), length(forrest_h3k27ac_genes),
    length(npc_h3k4me3_genes), length(mcclay_h3k4me3_genes), length(forrest_h3k4me3_genes),
    length(npc_h3k27me3_genes), length(mcclay_h3k27me3_genes), length(forrest_h3k27me3_genes)
  )
)

histone_gene_summary$Percent_TCF4_Genes_With_Mark <- round(
  histone_gene_summary$Mark_Overlapping_TCF4_Genes / histone_gene_summary$Total_TCF4_Genes * 100,
  3
)

histone_gene_summary
write.csv(histone_gene_summary, "histone_gene_overlap_summary_hg38.csv", row.names = FALSE)

############################################################
# 6. EP300 overlap analysis, hg38
############################################################

ep300_optimal_df <- read.table(
  "GSE127584_ENCFF459ARL_optimal_idr_thresholded_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

ep300_optimal_df <- ep300_optimal_df[, 1:3]
colnames(ep300_optimal_df) <- c("chr", "start", "end")
ep300_optimal_df$start <- as.numeric(ep300_optimal_df$start)
ep300_optimal_df$end <- as.numeric(ep300_optimal_df$end)
ep300_optimal_df <- ep300_optimal_df[!is.na(ep300_optimal_df$start) & !is.na(ep300_optimal_df$end) & ep300_optimal_df$end > ep300_optimal_df$start, ]
ep300_optimal_peaks <- GRanges(seqnames = ep300_optimal_df$chr, ranges = IRanges(start = ep300_optimal_df$start + 1, end = ep300_optimal_df$end))
seqlevelsStyle(ep300_optimal_peaks) <- "UCSC"
ep300_optimal_peaks <- sort(keepStandardChromosomes(ep300_optimal_peaks, pruning.mode = "coarse"))

ep300_conservative_df <- read.table(
  "GSE127584_ENCFF747GDL_conservative_idr_thresholded_peaks_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

ep300_conservative_df <- ep300_conservative_df[, 1:3]
colnames(ep300_conservative_df) <- c("chr", "start", "end")
ep300_conservative_df$start <- as.numeric(ep300_conservative_df$start)
ep300_conservative_df$end <- as.numeric(ep300_conservative_df$end)
ep300_conservative_df <- ep300_conservative_df[!is.na(ep300_conservative_df$start) & !is.na(ep300_conservative_df$end) & ep300_conservative_df$end > ep300_conservative_df$start, ]
ep300_conservative_peaks <- GRanges(seqnames = ep300_conservative_df$chr, ranges = IRanges(start = ep300_conservative_df$start + 1, end = ep300_conservative_df$end))
seqlevelsStyle(ep300_conservative_peaks) <- "UCSC"
ep300_conservative_peaks <- sort(keepStandardChromosomes(ep300_conservative_peaks, pruning.mode = "coarse"))

npc_ep300_optimal_hits <- findOverlaps(npc_peaks, ep300_optimal_peaks)
mcclay_ep300_optimal_hits <- findOverlaps(mcclay_peaks, ep300_optimal_peaks)
forrest_ep300_optimal_hits <- findOverlaps(forrest_peaks, ep300_optimal_peaks)
npc_ep300_conservative_hits <- findOverlaps(npc_peaks, ep300_conservative_peaks)
mcclay_ep300_conservative_hits <- findOverlaps(mcclay_peaks, ep300_conservative_peaks)
forrest_ep300_conservative_hits <- findOverlaps(forrest_peaks, ep300_conservative_peaks)

ep300_peak_overlap_summary <- data.frame(
  EP300_Set = c("Optimal", "Optimal", "Optimal", "Conservative", "Conservative", "Conservative"),
  TCF4_Dataset = c("NPC_TCF4", "McClay_TCF4", "Forrest_TCF4", "NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"),
  TCF4_Total_Peaks = c(
    length(npc_peaks), length(mcclay_peaks), length(forrest_peaks),
    length(npc_peaks), length(mcclay_peaks), length(forrest_peaks)
  ),
  EP300_Total_Peaks = c(
    length(ep300_optimal_peaks), length(ep300_optimal_peaks), length(ep300_optimal_peaks),
    length(ep300_conservative_peaks), length(ep300_conservative_peaks), length(ep300_conservative_peaks)
  ),
  TCF4_Overlapping_Peaks = c(
    length(unique(queryHits(npc_ep300_optimal_hits))),
    length(unique(queryHits(mcclay_ep300_optimal_hits))),
    length(unique(queryHits(forrest_ep300_optimal_hits))),
    length(unique(queryHits(npc_ep300_conservative_hits))),
    length(unique(queryHits(mcclay_ep300_conservative_hits))),
    length(unique(queryHits(forrest_ep300_conservative_hits)))
  ),
  EP300_Overlapping_Peaks = c(
    length(unique(subjectHits(npc_ep300_optimal_hits))),
    length(unique(subjectHits(mcclay_ep300_optimal_hits))),
    length(unique(subjectHits(forrest_ep300_optimal_hits))),
    length(unique(subjectHits(npc_ep300_conservative_hits))),
    length(unique(subjectHits(mcclay_ep300_conservative_hits))),
    length(unique(subjectHits(forrest_ep300_conservative_hits)))
  )
)

ep300_peak_overlap_summary$Percent_TCF4_Peaks_Overlapping_EP300 <- round(
  ep300_peak_overlap_summary$TCF4_Overlapping_Peaks / ep300_peak_overlap_summary$TCF4_Total_Peaks * 100,
  3
)

ep300_peak_overlap_summary$Percent_EP300_Peaks_Overlapping_TCF4 <- round(
  ep300_peak_overlap_summary$EP300_Overlapping_Peaks / ep300_peak_overlap_summary$EP300_Total_Peaks * 100,
  3
)

ep300_peak_overlap_summary
write.csv(ep300_peak_overlap_summary, "ep300_peak_overlap_summary_hg38.csv", row.names = FALSE)

npc_ep300_optimal_genes <- unique(na.omit(as.data.frame(annotatePeak(npc_peaks[unique(queryHits(npc_ep300_optimal_hits))], TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
mcclay_ep300_optimal_genes <- unique(na.omit(as.data.frame(annotatePeak(mcclay_peaks[unique(queryHits(mcclay_ep300_optimal_hits))], TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
forrest_ep300_optimal_genes <- unique(na.omit(as.data.frame(annotatePeak(forrest_peaks[unique(queryHits(forrest_ep300_optimal_hits))], TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))

npc_ep300_conservative_genes <- unique(na.omit(as.data.frame(annotatePeak(npc_peaks[unique(queryHits(npc_ep300_conservative_hits))], TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
mcclay_ep300_conservative_genes <- unique(na.omit(as.data.frame(annotatePeak(mcclay_peaks[unique(queryHits(mcclay_ep300_conservative_hits))], TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))
forrest_ep300_conservative_genes <- unique(na.omit(as.data.frame(annotatePeak(forrest_peaks[unique(queryHits(forrest_ep300_conservative_hits))], TxDb = txdb_hg38, annoDb = "org.Hs.eg.db"))$SYMBOL))

ep300_gene_overlap_summary <- data.frame(
  EP300_Set = c("Optimal", "Optimal", "Optimal", "Optimal", "Conservative", "Conservative", "Conservative", "Conservative"),
  Comparison = c(
    "NPC_vs_McClay",
    "NPC_vs_Forrest",
    "McClay_vs_Forrest",
    "NPC_vs_McClay_vs_Forrest",
    "NPC_vs_McClay",
    "NPC_vs_Forrest",
    "McClay_vs_Forrest",
    "NPC_vs_McClay_vs_Forrest"
  ),
  Overlap_Genes = c(
    length(intersect(npc_ep300_optimal_genes, mcclay_ep300_optimal_genes)),
    length(intersect(npc_ep300_optimal_genes, forrest_ep300_optimal_genes)),
    length(intersect(mcclay_ep300_optimal_genes, forrest_ep300_optimal_genes)),
    length(Reduce(intersect, list(npc_ep300_optimal_genes, mcclay_ep300_optimal_genes, forrest_ep300_optimal_genes))),
    length(intersect(npc_ep300_conservative_genes, mcclay_ep300_conservative_genes)),
    length(intersect(npc_ep300_conservative_genes, forrest_ep300_conservative_genes)),
    length(intersect(mcclay_ep300_conservative_genes, forrest_ep300_conservative_genes)),
    length(Reduce(intersect, list(npc_ep300_conservative_genes, mcclay_ep300_conservative_genes, forrest_ep300_conservative_genes)))
  )
)

ep300_gene_overlap_summary
write.csv(ep300_gene_overlap_summary, "ep300_gene_overlap_summary_hg38.csv", row.names = FALSE)

############################################################
# 7. GSE79965 H3K27ac, CHD8, and ATAC overlap analysis
############################################################

# GSE79965 coordinates are tested against both hg38 and hg19 TCF4 peak sets.
# The hg19 section below is the preferred coordinate system for the TAD and RNA-seq work.

gse_h3k27ac <- read_xlsx("GSE79965_processedData.xlsx", sheet = "27AC_peaks", skip = 2, col_names = c("chr", "start", "end")) %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start), !is.na(end), end > start)

gse_chd8_neuron <- read_xlsx("GSE79965_processedData.xlsx", sheet = "neuronCHD8_peaks", skip = 2, col_names = c("chr", "start", "end")) %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start), !is.na(end), end > start)

gse_chd8_brain <- read_xlsx("GSE79965_processedData.xlsx", sheet = "brainCHD8_peaks", skip = 2, col_names = c("chr", "start", "end")) %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start), !is.na(end), end > start)

gse_atac <- read_xlsx("GSE79965_processedData.xlsx", sheet = "ATAC_Seq_peaks", skip = 3, col_names = c("chr", "start", "end")) %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start), !is.na(end), end > start)

gse_h3k27ac_peaks <- GRanges(seqnames = gse_h3k27ac$chr, ranges = IRanges(start = gse_h3k27ac$start + 1, end = gse_h3k27ac$end))
gse_chd8_neuron_peaks <- GRanges(seqnames = gse_chd8_neuron$chr, ranges = IRanges(start = gse_chd8_neuron$start + 1, end = gse_chd8_neuron$end))
gse_chd8_brain_peaks <- GRanges(seqnames = gse_chd8_brain$chr, ranges = IRanges(start = gse_chd8_brain$start + 1, end = gse_chd8_brain$end))
gse_atac_peaks <- GRanges(seqnames = gse_atac$chr, ranges = IRanges(start = gse_atac$start + 1, end = gse_atac$end))

seqlevelsStyle(gse_h3k27ac_peaks) <- "UCSC"
seqlevelsStyle(gse_chd8_neuron_peaks) <- "UCSC"
seqlevelsStyle(gse_chd8_brain_peaks) <- "UCSC"
seqlevelsStyle(gse_atac_peaks) <- "UCSC"

gse_h3k27ac_peaks <- sort(keepStandardChromosomes(gse_h3k27ac_peaks, pruning.mode = "coarse"))
gse_chd8_neuron_peaks <- sort(keepStandardChromosomes(gse_chd8_neuron_peaks, pruning.mode = "coarse"))
gse_chd8_brain_peaks <- sort(keepStandardChromosomes(gse_chd8_brain_peaks, pruning.mode = "coarse"))
gse_atac_peaks <- sort(keepStandardChromosomes(gse_atac_peaks, pruning.mode = "coarse"))

# hg38 TCF4 comparison.
npc_gse_h3k27ac_hg38_hits <- findOverlaps(npc_peaks, gse_h3k27ac_peaks)
npc_gse_chd8_neuron_hg38_hits <- findOverlaps(npc_peaks, gse_chd8_neuron_peaks)
npc_gse_chd8_brain_hg38_hits <- findOverlaps(npc_peaks, gse_chd8_brain_peaks)
npc_gse_atac_hg38_hits <- findOverlaps(npc_peaks, gse_atac_peaks)

mcclay_gse_h3k27ac_hg38_hits <- findOverlaps(mcclay_peaks, gse_h3k27ac_peaks)
mcclay_gse_chd8_neuron_hg38_hits <- findOverlaps(mcclay_peaks, gse_chd8_neuron_peaks)
mcclay_gse_chd8_brain_hg38_hits <- findOverlaps(mcclay_peaks, gse_chd8_brain_peaks)
mcclay_gse_atac_hg38_hits <- findOverlaps(mcclay_peaks, gse_atac_peaks)

forrest_gse_h3k27ac_hg38_hits <- findOverlaps(forrest_peaks, gse_h3k27ac_peaks)
forrest_gse_chd8_neuron_hg38_hits <- findOverlaps(forrest_peaks, gse_chd8_neuron_peaks)
forrest_gse_chd8_brain_hg38_hits <- findOverlaps(forrest_peaks, gse_chd8_brain_peaks)
forrest_gse_atac_hg38_hits <- findOverlaps(forrest_peaks, gse_atac_peaks)

gse79965_hg38_overlap_summary <- data.frame(
  Feature = rep(c("Neuron_H3K27ac", "Neuron_CHD8", "Brain_CHD8", "Neuron_ATAC"), 3),
  TCF4_Dataset = c(rep("NPC_TCF4_hg38", 4), rep("McClay_TCF4_hg38", 4), rep("Forrest_TCF4_hg38", 4)),
  TCF4_Total_Peaks = c(rep(length(npc_peaks), 4), rep(length(mcclay_peaks), 4), rep(length(forrest_peaks), 4)),
  Feature_Total_Peaks = rep(c(length(gse_h3k27ac_peaks), length(gse_chd8_neuron_peaks), length(gse_chd8_brain_peaks), length(gse_atac_peaks)), 3),
  TCF4_Overlapping_Peaks = c(
    length(unique(queryHits(npc_gse_h3k27ac_hg38_hits))),
    length(unique(queryHits(npc_gse_chd8_neuron_hg38_hits))),
    length(unique(queryHits(npc_gse_chd8_brain_hg38_hits))),
    length(unique(queryHits(npc_gse_atac_hg38_hits))),
    length(unique(queryHits(mcclay_gse_h3k27ac_hg38_hits))),
    length(unique(queryHits(mcclay_gse_chd8_neuron_hg38_hits))),
    length(unique(queryHits(mcclay_gse_chd8_brain_hg38_hits))),
    length(unique(queryHits(mcclay_gse_atac_hg38_hits))),
    length(unique(queryHits(forrest_gse_h3k27ac_hg38_hits))),
    length(unique(queryHits(forrest_gse_chd8_neuron_hg38_hits))),
    length(unique(queryHits(forrest_gse_chd8_brain_hg38_hits))),
    length(unique(queryHits(forrest_gse_atac_hg38_hits)))
  ),
  Feature_Overlapping_Peaks = c(
    length(unique(subjectHits(npc_gse_h3k27ac_hg38_hits))),
    length(unique(subjectHits(npc_gse_chd8_neuron_hg38_hits))),
    length(unique(subjectHits(npc_gse_chd8_brain_hg38_hits))),
    length(unique(subjectHits(npc_gse_atac_hg38_hits))),
    length(unique(subjectHits(mcclay_gse_h3k27ac_hg38_hits))),
    length(unique(subjectHits(mcclay_gse_chd8_neuron_hg38_hits))),
    length(unique(subjectHits(mcclay_gse_chd8_brain_hg38_hits))),
    length(unique(subjectHits(mcclay_gse_atac_hg38_hits))),
    length(unique(subjectHits(forrest_gse_h3k27ac_hg38_hits))),
    length(unique(subjectHits(forrest_gse_chd8_neuron_hg38_hits))),
    length(unique(subjectHits(forrest_gse_chd8_brain_hg38_hits))),
    length(unique(subjectHits(forrest_gse_atac_hg38_hits)))
  )
)

gse79965_hg38_overlap_summary$Percent_TCF4_Peaks_Overlapping_Feature <- round(
  gse79965_hg38_overlap_summary$TCF4_Overlapping_Peaks / gse79965_hg38_overlap_summary$TCF4_Total_Peaks * 100,
  3
)

gse79965_hg38_overlap_summary$Percent_Feature_Peaks_Overlapping_TCF4 <- round(
  gse79965_hg38_overlap_summary$Feature_Overlapping_Peaks / gse79965_hg38_overlap_summary$Feature_Total_Peaks * 100,
  3
)

gse79965_hg38_overlap_summary
write.csv(gse79965_hg38_overlap_summary, "gse79965_hg38_overlap_summary.csv", row.names = FALSE)

# hg19 TCF4 comparison.
npc_hg19 <- import("NPC_ab21_idr_500bp_summit.bed")
forrest_hg19 <- import("Forrest_TCF4_IDR_sorted_hg19.bed")
mcclay_hg19_df <- read.csv("McClay_TCF4_11322_consensus_hg19_annotated.csv")

mcclay_hg19 <- GRanges(
  seqnames = mcclay_hg19_df$chrom,
  ranges = IRanges(start = mcclay_hg19_df$start, end = mcclay_hg19_df$end)
)

seqlevelsStyle(npc_hg19) <- "UCSC"
seqlevelsStyle(mcclay_hg19) <- "UCSC"
seqlevelsStyle(forrest_hg19) <- "UCSC"

npc_hg19 <- sort(keepStandardChromosomes(npc_hg19, pruning.mode = "coarse"))
mcclay_hg19 <- sort(keepStandardChromosomes(mcclay_hg19, pruning.mode = "coarse"))
forrest_hg19 <- sort(keepStandardChromosomes(forrest_hg19, pruning.mode = "coarse"))

npc_gse_h3k27ac_hg19_hits <- findOverlaps(npc_hg19, gse_h3k27ac_peaks)
npc_gse_chd8_neuron_hg19_hits <- findOverlaps(npc_hg19, gse_chd8_neuron_peaks)
npc_gse_chd8_brain_hg19_hits <- findOverlaps(npc_hg19, gse_chd8_brain_peaks)
npc_gse_atac_hg19_hits <- findOverlaps(npc_hg19, gse_atac_peaks)

mcclay_gse_h3k27ac_hg19_hits <- findOverlaps(mcclay_hg19, gse_h3k27ac_peaks)
mcclay_gse_chd8_neuron_hg19_hits <- findOverlaps(mcclay_hg19, gse_chd8_neuron_peaks)
mcclay_gse_chd8_brain_hg19_hits <- findOverlaps(mcclay_hg19, gse_chd8_brain_peaks)
mcclay_gse_atac_hg19_hits <- findOverlaps(mcclay_hg19, gse_atac_peaks)

forrest_gse_h3k27ac_hg19_hits <- findOverlaps(forrest_hg19, gse_h3k27ac_peaks)
forrest_gse_chd8_neuron_hg19_hits <- findOverlaps(forrest_hg19, gse_chd8_neuron_peaks)
forrest_gse_chd8_brain_hg19_hits <- findOverlaps(forrest_hg19, gse_chd8_brain_peaks)
forrest_gse_atac_hg19_hits <- findOverlaps(forrest_hg19, gse_atac_peaks)

gse79965_hg19_overlap_summary <- data.frame(
  Feature = rep(c("Neuron_H3K27ac", "Neuron_CHD8", "Brain_CHD8", "Neuron_ATAC"), 3),
  TCF4_Dataset = c(rep("NPC_TCF4_hg19", 4), rep("McClay_TCF4_hg19", 4), rep("Forrest_TCF4_hg19", 4)),
  TCF4_Total_Peaks = c(rep(length(npc_hg19), 4), rep(length(mcclay_hg19), 4), rep(length(forrest_hg19), 4)),
  Feature_Total_Peaks = rep(c(length(gse_h3k27ac_peaks), length(gse_chd8_neuron_peaks), length(gse_chd8_brain_peaks), length(gse_atac_peaks)), 3),
  TCF4_Overlapping_Peaks = c(
    length(unique(queryHits(npc_gse_h3k27ac_hg19_hits))),
    length(unique(queryHits(npc_gse_chd8_neuron_hg19_hits))),
    length(unique(queryHits(npc_gse_chd8_brain_hg19_hits))),
    length(unique(queryHits(npc_gse_atac_hg19_hits))),
    length(unique(queryHits(mcclay_gse_h3k27ac_hg19_hits))),
    length(unique(queryHits(mcclay_gse_chd8_neuron_hg19_hits))),
    length(unique(queryHits(mcclay_gse_chd8_brain_hg19_hits))),
    length(unique(queryHits(mcclay_gse_atac_hg19_hits))),
    length(unique(queryHits(forrest_gse_h3k27ac_hg19_hits))),
    length(unique(queryHits(forrest_gse_chd8_neuron_hg19_hits))),
    length(unique(queryHits(forrest_gse_chd8_brain_hg19_hits))),
    length(unique(queryHits(forrest_gse_atac_hg19_hits)))
  ),
  Feature_Overlapping_Peaks = c(
    length(unique(subjectHits(npc_gse_h3k27ac_hg19_hits))),
    length(unique(subjectHits(npc_gse_chd8_neuron_hg19_hits))),
    length(unique(subjectHits(npc_gse_chd8_brain_hg19_hits))),
    length(unique(subjectHits(npc_gse_atac_hg19_hits))),
    length(unique(subjectHits(mcclay_gse_h3k27ac_hg19_hits))),
    length(unique(subjectHits(mcclay_gse_chd8_neuron_hg19_hits))),
    length(unique(subjectHits(mcclay_gse_chd8_brain_hg19_hits))),
    length(unique(subjectHits(mcclay_gse_atac_hg19_hits))),
    length(unique(subjectHits(forrest_gse_h3k27ac_hg19_hits))),
    length(unique(subjectHits(forrest_gse_chd8_neuron_hg19_hits))),
    length(unique(subjectHits(forrest_gse_chd8_brain_hg19_hits))),
    length(unique(subjectHits(forrest_gse_atac_hg19_hits)))
  )
)

gse79965_hg19_overlap_summary$Percent_TCF4_Peaks_Overlapping_Feature <- round(
  gse79965_hg19_overlap_summary$TCF4_Overlapping_Peaks / gse79965_hg19_overlap_summary$TCF4_Total_Peaks * 100,
  3
)

gse79965_hg19_overlap_summary$Percent_Feature_Peaks_Overlapping_TCF4 <- round(
  gse79965_hg19_overlap_summary$Feature_Overlapping_Peaks / gse79965_hg19_overlap_summary$Feature_Total_Peaks * 100,
  3
)

gse79965_hg19_overlap_summary
write.csv(gse79965_hg19_overlap_summary, "gse79965_hg19_overlap_summary.csv", row.names = FALSE)

############################################################
# 8. GSE79965 Hi-C TAD analysis, hg19
############################################################

hic_df <- read_xlsx("GSE79965_processedData.xlsx", sheet = "Hi-C", skip = 2, col_names = c("chr", "start", "end")) %>%
  filter(!is.na(chr), chr != "Chromosome") %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  filter(!is.na(start), !is.na(end), end > start)

hic_tads <- GRanges(
  seqnames = hic_df$chr,
  ranges = IRanges(start = hic_df$start + 1, end = hic_df$end)
)

seqlevelsStyle(hic_tads) <- "UCSC"
hic_tads <- sort(keepStandardChromosomes(hic_tads, pruning.mode = "coarse"))

tad_summary <- data.frame(
  TAD_ID = paste0("TAD_", seq_along(hic_tads)),
  Chromosome = as.character(seqnames(hic_tads)),
  Start = start(hic_tads),
  End = end(hic_tads),
  Width_bp = width(hic_tads),
  NPC_TCF4_Peaks = countOverlaps(hic_tads, npc_hg19),
  McClay_TCF4_Peaks = countOverlaps(hic_tads, mcclay_hg19),
  Forrest_TCF4_Peaks = countOverlaps(hic_tads, forrest_hg19),
  H3K27ac_Peaks = countOverlaps(hic_tads, gse_h3k27ac_peaks),
  Neuron_CHD8_Peaks = countOverlaps(hic_tads, gse_chd8_neuron_peaks),
  Brain_CHD8_Peaks = countOverlaps(hic_tads, gse_chd8_brain_peaks),
  ATAC_Peaks = countOverlaps(hic_tads, gse_atac_peaks)
)

tad_summary$Has_NPC_TCF4 <- tad_summary$NPC_TCF4_Peaks > 0
tad_summary$Has_McClay_TCF4 <- tad_summary$McClay_TCF4_Peaks > 0
tad_summary$Has_Forrest_TCF4 <- tad_summary$Forrest_TCF4_Peaks > 0
tad_summary$Has_H3K27ac <- tad_summary$H3K27ac_Peaks > 0
tad_summary$Has_ATAC <- tad_summary$ATAC_Peaks > 0
tad_summary$Has_Neuron_CHD8 <- tad_summary$Neuron_CHD8_Peaks > 0
tad_summary$Has_Brain_CHD8 <- tad_summary$Brain_CHD8_Peaks > 0
tad_summary$Has_Any_CHD8 <- tad_summary$Has_Neuron_CHD8 | tad_summary$Has_Brain_CHD8

tad_summary$NPC_Active_TAD <- tad_summary$Has_NPC_TCF4 & tad_summary$Has_H3K27ac & tad_summary$Has_ATAC
tad_summary$NPC_CHD8_TAD <- tad_summary$Has_NPC_TCF4 & tad_summary$Has_Any_CHD8
tad_summary$NPC_Active_CHD8_TAD <- tad_summary$Has_NPC_TCF4 & tad_summary$Has_H3K27ac & tad_summary$Has_ATAC & tad_summary$Has_Any_CHD8
tad_summary$Shared_TCF4_TAD <- tad_summary$Has_NPC_TCF4 & tad_summary$Has_McClay_TCF4 & tad_summary$Has_Forrest_TCF4

tad_category_summary <- data.frame(
  Category = c(
    "Total_HiC_TADs",
    "TADs_with_NPC_TCF4",
    "TADs_with_McClay_TCF4",
    "TADs_with_Forrest_TCF4",
    "TADs_with_H3K27ac",
    "TADs_with_ATAC",
    "TADs_with_Neuron_CHD8",
    "TADs_with_Brain_CHD8",
    "TADs_with_Any_CHD8",
    "NPC_TCF4_H3K27ac_ATAC_TADs",
    "NPC_TCF4_CHD8_TADs",
    "NPC_TCF4_H3K27ac_ATAC_CHD8_TADs",
    "Shared_NPC_McClay_Forrest_TCF4_TADs"
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

tad_category_summary$Percent_of_All_TADs <- round(tad_category_summary$TAD_Count / nrow(tad_summary) * 100, 3)
tad_category_summary
write.csv(tad_category_summary, "tad_category_summary_gse79965_hg19.csv", row.names = FALSE)

npc_tads <- tad_summary$TAD_ID[tad_summary$Has_NPC_TCF4]
mcclay_tads <- tad_summary$TAD_ID[tad_summary$Has_McClay_TCF4]
forrest_tads <- tad_summary$TAD_ID[tad_summary$Has_Forrest_TCF4]

tcf4_tad_overlap_summary <- data.frame(
  Comparison = c(
    "NPC_TCF4_TADs_vs_McClay_TCF4_TADs",
    "NPC_TCF4_TADs_vs_Forrest_TCF4_TADs",
    "McClay_TCF4_TADs_vs_Forrest_TCF4_TADs",
    "NPC_vs_McClay_vs_Forrest_TCF4_TADs"
  ),
  Shared_TADs = c(
    length(intersect(npc_tads, mcclay_tads)),
    length(intersect(npc_tads, forrest_tads)),
    length(intersect(mcclay_tads, forrest_tads)),
    length(Reduce(intersect, list(npc_tads, mcclay_tads, forrest_tads)))
  ),
  Percent_of_First_TAD_Set = c(
    round(length(intersect(npc_tads, mcclay_tads)) / length(npc_tads) * 100, 3),
    round(length(intersect(npc_tads, forrest_tads)) / length(npc_tads) * 100, 3),
    round(length(intersect(mcclay_tads, forrest_tads)) / length(mcclay_tads) * 100, 3),
    round(length(Reduce(intersect, list(npc_tads, mcclay_tads, forrest_tads))) / length(npc_tads) * 100, 3)
  ),
  Percent_of_Second_TAD_Set = c(
    round(length(intersect(npc_tads, mcclay_tads)) / length(mcclay_tads) * 100, 3),
    round(length(intersect(npc_tads, forrest_tads)) / length(forrest_tads) * 100, 3),
    round(length(intersect(mcclay_tads, forrest_tads)) / length(forrest_tads) * 100, 3),
    NA
  )
)

tcf4_tad_overlap_summary
write.csv(tcf4_tad_overlap_summary, "tcf4_tad_overlap_summary_gse79965_hg19.csv", row.names = FALSE)

tad_summary$NPC_Regulatory_Score <- tad_summary$NPC_TCF4_Peaks + tad_summary$H3K27ac_Peaks + tad_summary$ATAC_Peaks + tad_summary$Neuron_CHD8_Peaks + tad_summary$Brain_CHD8_Peaks
tad_summary$McClay_Regulatory_Score <- tad_summary$McClay_TCF4_Peaks + tad_summary$H3K27ac_Peaks + tad_summary$ATAC_Peaks + tad_summary$Neuron_CHD8_Peaks + tad_summary$Brain_CHD8_Peaks
tad_summary$Forrest_Regulatory_Score <- tad_summary$Forrest_TCF4_Peaks + tad_summary$H3K27ac_Peaks + tad_summary$ATAC_Peaks + tad_summary$Neuron_CHD8_Peaks + tad_summary$Brain_CHD8_Peaks

top_npc_regulatory_tads <- tad_summary %>% arrange(desc(NPC_Regulatory_Score))
top_mcclay_regulatory_tads <- tad_summary %>% arrange(desc(McClay_Regulatory_Score))
top_forrest_regulatory_tads <- tad_summary %>% arrange(desc(Forrest_Regulatory_Score))
npc_active_chd8_tads <- tad_summary %>% filter(NPC_Active_CHD8_TAD) %>% arrange(desc(NPC_TCF4_Peaks), desc(H3K27ac_Peaks), desc(ATAC_Peaks), desc(Neuron_CHD8_Peaks + Brain_CHD8_Peaks))

write.csv(top_npc_regulatory_tads, "top_npc_regulatory_tads_gse79965_hg19.csv", row.names = FALSE)
write.csv(top_mcclay_regulatory_tads, "top_mcclay_regulatory_tads_gse79965_hg19.csv", row.names = FALSE)
write.csv(top_forrest_regulatory_tads, "top_forrest_regulatory_tads_gse79965_hg19.csv", row.names = FALSE)
write.csv(npc_active_chd8_tads, "npc_active_chd8_tads_gse79965_hg19.csv", row.names = FALSE)

top_npc_tcf4_tads <- tad_summary %>% arrange(desc(NPC_TCF4_Peaks)) %>% slice_head(n = 100)
top_mcclay_tcf4_tads <- tad_summary %>% arrange(desc(McClay_TCF4_Peaks)) %>% slice_head(n = 100)
top_forrest_tcf4_tads <- tad_summary %>% arrange(desc(Forrest_TCF4_Peaks)) %>% slice_head(n = 100)
top_h3k27ac_tads <- tad_summary %>% arrange(desc(H3K27ac_Peaks)) %>% slice_head(n = 100)
top_atac_tads <- tad_summary %>% arrange(desc(ATAC_Peaks)) %>% slice_head(n = 100)
top_neuron_chd8_tads <- tad_summary %>% arrange(desc(Neuron_CHD8_Peaks)) %>% slice_head(n = 100)
top_brain_chd8_tads <- tad_summary %>% arrange(desc(Brain_CHD8_Peaks)) %>% slice_head(n = 100)

top100_tcf4_feature_overlap_summary <- data.frame(
  TCF4_Dataset = c(rep("NPC_TCF4", 4), rep("McClay_TCF4", 4), rep("Forrest_TCF4", 4)),
  Feature = rep(c("H3K27ac", "ATAC", "Neuron_CHD8", "Brain_CHD8"), 3),
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

top100_tcf4_feature_overlap_summary$Percent_Shared <- top100_tcf4_feature_overlap_summary$Shared_Top100_TADs
top100_tcf4_feature_overlap_summary
write.csv(top100_tcf4_feature_overlap_summary, "top100_tcf4_feature_tad_overlap_summary.csv", row.names = FALSE)

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

tcf4_feature_correlation_summary <- data.frame(
  TCF4_Dataset = c("NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"),
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

tcf4_feature_correlation_summary[, 2:5] <- round(tcf4_feature_correlation_summary[, 2:5], 3)
tcf4_feature_correlation_summary
write.csv(tcf4_feature_correlation_summary, "tcf4_feature_correlation_summary_tads.csv", row.names = FALSE)

############################################################
# 9. TAD_1375 gene lookup, hg19
############################################################

tad_1375 <- GRanges(seqnames = "chr19", ranges = IRanges(start = 57810001, end = 58980000))
hg19_genes <- genes(txdb_hg19)
seqlevelsStyle(hg19_genes) <- "UCSC"
hg19_genes <- keepStandardChromosomes(hg19_genes, pruning.mode = "coarse")

tad_1375_gene_hits <- findOverlaps(hg19_genes, tad_1375)
tad_1375_genes <- hg19_genes[unique(queryHits(tad_1375_gene_hits))]
tad_1375_entrez <- names(tad_1375_genes)
tad_1375_symbols <- mapIds(org.Hs.eg.db, keys = tad_1375_entrez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

tad_1375_gene_table <- data.frame(
  Entrez_ID = tad_1375_entrez,
  Gene_Symbol = as.character(tad_1375_symbols),
  Chromosome = as.character(seqnames(tad_1375_genes)),
  Start = start(tad_1375_genes),
  End = end(tad_1375_genes)
)

tad_1375_gene_table <- tad_1375_gene_table[!is.na(tad_1375_gene_table$Gene_Symbol), ]
tad_1375_gene_table
write.csv(tad_1375_gene_table, "tad_1375_gene_table_hg19.csv", row.names = FALSE)

############################################################
# 10. RNA-seq integration with TCF4 peaks and TADs, hg19
############################################################

rnaseq_raw <- read_xlsx("GSE79965_processedData.xlsx", sheet = "RNA-SEQ", col_names = FALSE)
colnames(rnaseq_raw) <- c("Gene_Symbol", "FPKM")

rnaseq_df <- rnaseq_raw %>%
  filter(!is.na(Gene_Symbol), Gene_Symbol != "Symbols") %>%
  mutate(FPKM = as.numeric(FPKM)) %>%
  filter(!is.na(FPKM))

expressed_genes_fpkm_0 <- rnaseq_df$Gene_Symbol[rnaseq_df$FPKM > 0]
expressed_genes_fpkm_0_5 <- rnaseq_df$Gene_Symbol[rnaseq_df$FPKM > 0.5]
expressed_genes_fpkm_1 <- rnaseq_df$Gene_Symbol[rnaseq_df$FPKM > 1]

npc_hg19_anno_df <- as.data.frame(annotatePeak(npc_hg19, TxDb = txdb_hg19, annoDb = "org.Hs.eg.db"))
mcclay_hg19_anno_df <- as.data.frame(annotatePeak(mcclay_hg19, TxDb = txdb_hg19, annoDb = "org.Hs.eg.db"))
forrest_hg19_anno_df <- as.data.frame(annotatePeak(forrest_hg19, TxDb = txdb_hg19, annoDb = "org.Hs.eg.db"))

npc_hg19_genes <- unique(na.omit(npc_hg19_anno_df$SYMBOL))
mcclay_hg19_genes <- unique(na.omit(mcclay_hg19_anno_df$SYMBOL))
forrest_hg19_genes <- unique(na.omit(forrest_hg19_anno_df$SYMBOL))

tcf4_gene_expression_summary <- data.frame(
  Dataset = c("NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"),
  Total_TCF4_Genes = c(length(npc_hg19_genes), length(mcclay_hg19_genes), length(forrest_hg19_genes)),
  Expressed_Genes_FPKM_gt_0 = c(
    length(intersect(npc_hg19_genes, expressed_genes_fpkm_0)),
    length(intersect(mcclay_hg19_genes, expressed_genes_fpkm_0)),
    length(intersect(forrest_hg19_genes, expressed_genes_fpkm_0))
  ),
  Expressed_Genes_FPKM_gt_0_5 = c(
    length(intersect(npc_hg19_genes, expressed_genes_fpkm_0_5)),
    length(intersect(mcclay_hg19_genes, expressed_genes_fpkm_0_5)),
    length(intersect(forrest_hg19_genes, expressed_genes_fpkm_0_5))
  ),
  Expressed_Genes_FPKM_gt_1 = c(
    length(intersect(npc_hg19_genes, expressed_genes_fpkm_1)),
    length(intersect(mcclay_hg19_genes, expressed_genes_fpkm_1)),
    length(intersect(forrest_hg19_genes, expressed_genes_fpkm_1))
  )
)

tcf4_gene_expression_summary$Percent_Expressed_FPKM_gt_1 <- round(
  tcf4_gene_expression_summary$Expressed_Genes_FPKM_gt_1 / tcf4_gene_expression_summary$Total_TCF4_Genes * 100,
  3
)

tcf4_gene_expression_summary
write.csv(tcf4_gene_expression_summary, "tcf4_gene_expression_summary_hg19.csv", row.names = FALSE)

rnaseq_df$NPC_TCF4_Target <- rnaseq_df$Gene_Symbol %in% npc_hg19_genes
rnaseq_df$McClay_TCF4_Target <- rnaseq_df$Gene_Symbol %in% mcclay_hg19_genes
rnaseq_df$Forrest_TCF4_Target <- rnaseq_df$Gene_Symbol %in% forrest_hg19_genes
rnaseq_df$log2_FPKM_plus1 <- log2(rnaseq_df$FPKM + 1)

npc_expression_test <- wilcox.test(log2_FPKM_plus1 ~ NPC_TCF4_Target, data = rnaseq_df)
mcclay_expression_test <- wilcox.test(log2_FPKM_plus1 ~ McClay_TCF4_Target, data = rnaseq_df)
forrest_expression_test <- wilcox.test(log2_FPKM_plus1 ~ Forrest_TCF4_Target, data = rnaseq_df)

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
  Wilcoxon_P_Value = c(npc_expression_test$p.value, mcclay_expression_test$p.value, forrest_expression_test$p.value)
)

expression_test_summary
write.csv(expression_test_summary, "tcf4_expression_test_summary_hg19.csv", row.names = FALSE)

hg19_genes_gr <- genes(txdb_hg19)
seqlevelsStyle(hg19_genes_gr) <- "UCSC"
hg19_genes_gr <- keepStandardChromosomes(hg19_genes_gr, pruning.mode = "coarse")
hg19_gene_entrez <- names(hg19_genes_gr)
hg19_gene_symbols <- mapIds(org.Hs.eg.db, keys = hg19_gene_entrez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
mcols(hg19_genes_gr)$ENTREZID <- hg19_gene_entrez
mcols(hg19_genes_gr)$SYMBOL <- hg19_gene_symbols
hg19_genes_gr <- hg19_genes_gr[!is.na(mcols(hg19_genes_gr)$SYMBOL)]

gene_tad_hits <- findOverlaps(hg19_genes_gr, hic_tads)
gene_tad_table <- data.frame(
  Gene_Symbol = mcols(hg19_genes_gr)$SYMBOL[queryHits(gene_tad_hits)],
  TAD_ID = tad_summary$TAD_ID[subjectHits(gene_tad_hits)]
)

gene_tad_table <- unique(gene_tad_table)

gene_tad_expression <- gene_tad_table %>%
  left_join(rnaseq_df, by = "Gene_Symbol")

tad_expression_summary <- gene_tad_expression %>%
  group_by(TAD_ID) %>%
  summarise(
    Genes_in_TAD = n_distinct(Gene_Symbol),
    Genes_with_RNAseq_Data = sum(!is.na(FPKM)),
    Expressed_Genes_FPKM_gt_1 = sum(FPKM > 1, na.rm = TRUE),
    Mean_FPKM = mean(FPKM, na.rm = TRUE),
    Median_FPKM = median(FPKM, na.rm = TRUE),
    Mean_log2_FPKM_plus1 = mean(log2_FPKM_plus1, na.rm = TRUE),
    Median_log2_FPKM_plus1 = median(log2_FPKM_plus1, na.rm = TRUE),
    .groups = "drop"
  )

tad_summary_expression <- tad_summary %>%
  left_join(tad_expression_summary, by = "TAD_ID")

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

tad_expression_correlation_summary[, 2:3] <- round(tad_expression_correlation_summary[, 2:3], 3)
tad_expression_correlation_summary
write.csv(tad_expression_correlation_summary, "tad_expression_correlation_summary_hg19.csv", row.names = FALSE)

tad_1375_expression <- tad_1375_gene_table %>%
  left_join(rnaseq_df, by = "Gene_Symbol") %>%
  arrange(desc(FPKM))

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

tad_1375_expression
tad_1375_expression_summary
write.csv(tad_1375_expression, "tad_1375_expression_hg19.csv", row.names = FALSE)
write.csv(tad_1375_expression_summary, "tad_1375_expression_summary_hg19.csv", row.names = FALSE)

npc_active_chd8_tad_ids <- tad_summary$TAD_ID[tad_summary$NPC_Active_CHD8_TAD]
npc_active_chd8_tad_genes <- unique(gene_tad_table$Gene_Symbol[gene_tad_table$TAD_ID %in% npc_active_chd8_tad_ids])

npc_active_chd8_tad_expression <- rnaseq_df %>%
  filter(Gene_Symbol %in% npc_active_chd8_tad_genes) %>%
  arrange(desc(FPKM))

npc_active_chd8_tad_expression_summary <- data.frame(
  Gene_Set = "Genes_in_NPC_TCF4_H3K27ac_ATAC_CHD8_TADs",
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
