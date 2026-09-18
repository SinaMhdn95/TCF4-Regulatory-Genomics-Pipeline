###############################################################
# TCF4 peak overlap with Ripke PGC schizophrenia GWAS loci
# Source: Supplementary_Table_3_Bioinformatic_Summary_Data.xlsx
# Genome build: hg19 / GRCh37
###############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(GenomicRanges)
  library(rtracklayer)
})

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

input_file <- "Supplementary_Table_3_Bioinformatic_Summary_Data.xlsx"
input_sheet <- "Supplementary Table 3"

raw_loci <- read_xlsx(input_file, sheet = input_sheet)

# Clean Ripke loci
locus_clean <- raw_loci %>%
  rename(
    rank = Rank,
    p_value = `P-value`,
    position_hg19 = `Position (hg19)`,
    scz = SCZ,
    genes = `Protein coding genes`,
    omim = OMIM,
    nhgri = `NHGRI GWAS catalog`,
    ko = `KO phenotype`
  ) %>%
  mutate(
    position_hg19 = as.character(position_hg19),
    position_hg19 = str_replace(position_hg19, "^chr", ""),
    chr = str_replace(position_hg19, ":.*$", ""),
    start = as.integer(str_replace(position_hg19, "^.*:(\\d+)-.*$", "\\1")),
    end = as.integer(str_replace(position_hg19, "^.*-(\\d+)$", "\\1")),
    chr = paste0("chr", chr),
    gene_text = str_squish(
      str_replace_all(genes, "\\s*\\(.*?\\)", "")
    ),
    gene_list = str_split(gene_text, "\\s+"),
    corrected_start = pmin(start, end, na.rm = TRUE),
    corrected_end = pmax(start, end, na.rm = TRUE)
  )

ripke_invalid_positions <- locus_clean %>%
  filter(is.na(corrected_start) | is.na(corrected_end))

ripke_clean <- locus_clean %>%
  filter(!is.na(corrected_start) & !is.na(corrected_end)) %>%
  select(
    rank,
    p_value,
    chr,
    start = corrected_start,
    end = corrected_end,
    scz,
    genes,
    gene_text,
    gene_list,
    omim,
    nhgri,
    ko
  )

View(ripke_clean)

# Build GRanges for Ripke loci
ripke_loci <- makeGRangesFromDataFrame(
  ripke_clean,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

# Extract one gene per row
ripke_genes <- ripke_clean %>%
  unnest(gene_list) %>%
  filter(
    gene_list != "",
    !str_detect(gene_list, "^(Locus|too|broad|\\*)$")
  ) %>%
  mutate(
    gene = str_replace_all(gene_list, "[^A-Za-z0-9._-]", "")
  ) %>%
  filter(gene != "") %>%
  distinct(rank, gene, .keep_all = TRUE)

View(ripke_genes)

# Import hg19 TCF4 peak files
npc_hg19 <- import("NPC_ab21_idr_500bp_summit.bed")
forrest_hg19 <- import("Forrest_TCF4_IDR_sorted_hg19.bed")
mcclay_hg19_df <- read.csv("McClay_TCF4_11322_consensus_hg19_annotated.csv")

# Convert McClay hg19 data frame to GRanges
mcclay_hg19 <- makeGRangesFromDataFrame(
  mcclay_hg19_df,
  seqnames.field = "chrom",
  start.field = "start_hg19",
  end.field = "end_hg19",
  keep.extra.columns = TRUE
)

# Check peak counts
length(npc_hg19)
length(mcclay_hg19)
length(forrest_hg19)
length(ripke_loci)

# Find overlaps between TCF4 peaks and Ripke loci
npc_ripke_hits <- findOverlaps(npc_hg19, ripke_loci)
mcclay_ripke_hits <- findOverlaps(mcclay_hg19, ripke_loci)
forrest_ripke_hits <- findOverlaps(forrest_hg19, ripke_loci)

# Peak-level overlap summary
peak_overlap_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Total_Peaks = c(
    length(npc_hg19),
    length(mcclay_hg19),
    length(forrest_hg19)
  ),
  Peaks_in_SCZ_Loci = c(
    length(unique(queryHits(npc_ripke_hits))),
    length(unique(queryHits(mcclay_ripke_hits))),
    length(unique(queryHits(forrest_ripke_hits)))
  )
)

peak_overlap_summary$Percent_Peaks_in_SCZ_Loci <- round(
  peak_overlap_summary$Peaks_in_SCZ_Loci /
    peak_overlap_summary$Total_Peaks * 100,
  3
)

peak_overlap_summary

# Locus-level overlap summary
npc_loci <- unique(subjectHits(npc_ripke_hits))
mcclay_loci <- unique(subjectHits(mcclay_ripke_hits))
forrest_loci <- unique(subjectHits(forrest_ripke_hits))

locus_overlap_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Total_SCZ_Loci = length(ripke_loci),
  SCZ_Loci_With_TCF4_Peaks = c(
    length(npc_loci),
    length(mcclay_loci),
    length(forrest_loci)
  )
)

locus_overlap_summary$Percent_SCZ_Loci_With_TCF4_Peaks <- round(
  locus_overlap_summary$SCZ_Loci_With_TCF4_Peaks /
    locus_overlap_summary$Total_SCZ_Loci * 100,
  3
)

locus_overlap_summary

# Shared SCZ loci across all three TCF4 datasets
shared_scz_loci <- Reduce(
  intersect,
  list(
    npc_loci,
    mcclay_loci,
    forrest_loci
  )
)

shared_scz_loci_summary <- data.frame(
  Total_SCZ_Loci = length(ripke_loci),
  Shared_SCZ_Loci_With_TCF4_Peaks_All_Three = length(shared_scz_loci),
  Percent_Shared_SCZ_Loci_All_Three = round(
    length(shared_scz_loci) / length(ripke_loci) * 100,
    3
  )
)

shared_scz_loci_summary

# Table of shared SCZ loci
shared_scz_loci_table <- ripke_clean[shared_scz_loci, ]

shared_scz_loci_table

View(shared_scz_loci_table)


###############################################################
# TCF4 peak overlap with genomic loci implicates genes and synaptic biology in SCZ
# Source: Mapping genomic loci implicates genes and synaptic biology in schizophrenia
# Genome build: hg19 / GRCh37
###############################################################

# Extract TCF4-associated genes from hg19 peak files
# Already imported:
# npc_hg19 <- import("NPC_ab21_idr_500bp_summit.bed")
# forrest_hg19 <- import("Forrest_TCF4_IDR_sorted_hg19.bed")
# mcclay_hg19_df <- read.csv("McClay_TCF4_11322_consensus_hg19_annotated.csv")

# Convert McClay hg19 data frame to GRanges
mcclay_hg19 <- makeGRangesFromDataFrame(
  mcclay_hg19_df,
  seqnames.field = "chrom",
  start.field = "start_hg19",
  end.field = "end_hg19",
  keep.extra.columns = TRUE
)

# Annotate peaks to nearest genes
npc_anno <- annotatePeak(
  npc_hg19,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

mcclay_anno <- annotatePeak(
  mcclay_hg19,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

forrest_anno <- annotatePeak(
  forrest_hg19,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

# Convert annotation objects to data frames
npc_anno_df <- as.data.frame(npc_anno)
mcclay_anno_df <- as.data.frame(mcclay_anno)
forrest_anno_df <- as.data.frame(forrest_anno)


# Extract unique gene symbols
npc_genes <- npc_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

mcclay_genes <- mcclay_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

forrest_genes <- forrest_anno_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

length(npc_genes)
length(mcclay_genes)
length(forrest_genes)

# PGC3 SCZ Table S1 / LD clumps vs TCF4 peaks
# Read file
st1_raw <- read_xls("NIHMS1824605-supplement-Table_S1.xls")

# Clean column names manually if needed
colnames(st1_raw)

st1 <- st1_raw %>%
  rename(
    SNP = SNP,
    CHR = CHR,
    BP = BP,
    P = P,
    OR = OR,
    SE = SE,
    A1A2 = A1A2,
    INFO = INFO,
    ngt = ngt,
    range_left = `range.left`,
    range_right = `range.right`,
    span_kb = `span(kb)`
  ) %>%
  mutate(
    CHR = paste0("chr", CHR),
    BP = as.numeric(BP),
    P = as.numeric(P),
    OR = as.numeric(OR),
    SE = as.numeric(SE),
    range_left = as.numeric(range_left),
    range_right = as.numeric(range_right)
  ) %>%
  filter(!is.na(CHR), !is.na(BP), !is.na(P))

# Check total clumps
nrow(st1)

# Genome-wide significant index SNPs
st1_gws <- st1 %>%
  filter(P < 5e-8)

nrow(st1_gws)

# Convert exact index SNPs to GRanges
st1_gws_snp_gr <- GRanges(
  seqnames = st1_gws$CHR,
  ranges = IRanges(
    start = st1_gws$BP,
    end = st1_gws$BP
  )
)

mcols(st1_gws_snp_gr) <- st1_gws

# Convert LD clump regions to GRanges
st1_gws_region_gr <- GRanges(
  seqnames = st1_gws$CHR,
  ranges = IRanges(
    start = pmin(st1_gws$range_left, st1_gws$range_right),
    end = pmax(st1_gws$range_left, st1_gws$range_right)
  )
)

mcols(st1_gws_region_gr) <- st1_gws

# Find overlap with your TCF4 peaks
# Exact SNP overlap
npc_snp_hits <- findOverlaps(st1_gws_snp_gr, npc_hg19)
mcclay_snp_hits <- findOverlaps(st1_gws_snp_gr, mcclay_hg19)
forrest_snp_hits <- findOverlaps(st1_gws_snp_gr, forrest_hg19)

# LD region overlap
npc_region_hits <- findOverlaps(st1_gws_region_gr, npc_hg19)
mcclay_region_hits <- findOverlaps(st1_gws_region_gr, mcclay_hg19)
forrest_region_hits <- findOverlaps(st1_gws_region_gr, forrest_hg19)

# Summary table
gws_overlap_summary <- data.frame(
  Dataset = c("NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"),
  Total_GWS_Index_SNPs = length(st1_gws_snp_gr),
  GWS_Index_SNPs_Inside_TCF4_Peaks = c(
    length(unique(queryHits(npc_snp_hits))),
    length(unique(queryHits(mcclay_snp_hits))),
    length(unique(queryHits(forrest_snp_hits)))
  ),
  GWS_LD_Regions_Overlapping_TCF4_Peaks = c(
    length(unique(queryHits(npc_region_hits))),
    length(unique(queryHits(mcclay_region_hits))),
    length(unique(queryHits(forrest_region_hits)))
  ),
  TCF4_Peaks_Overlapping_GWS_LD_Regions = c(
    length(unique(subjectHits(npc_region_hits))),
    length(unique(subjectHits(mcclay_region_hits))),
    length(unique(subjectHits(forrest_region_hits)))
  )
)

gws_overlap_summary$Percent_GWS_Index_SNPs_Inside_TCF4_Peaks <- round(
  gws_overlap_summary$GWS_Index_SNPs_Inside_TCF4_Peaks /
    gws_overlap_summary$Total_GWS_Index_SNPs * 100,
  3
)

gws_overlap_summary$Percent_GWS_LD_Regions_Overlapping_TCF4 <- round(
  gws_overlap_summary$GWS_LD_Regions_Overlapping_TCF4_Peaks /
    gws_overlap_summary$Total_GWS_Index_SNPs * 100,
  3
)

gws_overlap_summary

# PGC3 / Trubetskoy Table S3 vs TCF4 peaks
# 287 distinct genome-wide significant regions
s3_raw <- read_xls(
  "NIHMS1824605-supplement-Table_S3.xls",
  sheet = "ST3 -  discovery + replication"
)

colnames(s3_raw)

# Clean Table S3
s3_clean <- s3_raw %>%
  rename(
    chr = Chromosome,
    top_index = `top-index`,
    top_pos = `top-pos`,
    top_alleles = `top-alleles`,
    top_freq = `top-freq`,
    top_info = `top-info`,
    top_p = `top-P`,
    top_or = `top-OR`,
    top_se = `top-SE`,
    merge_left = `merge-LEFT`,
    merge_right = `merge-RIGHT`,
    ensembl_n_genes = `N-genes`,
    ensembl_top5_genes = genes,
    ensembl_gene_names = genes_all
  ) %>%
  mutate(
    chr = paste0("chr", chr),
    top_pos = as.numeric(top_pos),
    top_p = as.numeric(top_p),
    top_or = as.numeric(top_or),
    top_se = as.numeric(top_se),
    merge_left = as.numeric(merge_left),
    merge_right = as.numeric(merge_right),
    region_start = pmin(merge_left, merge_right),
    region_end = pmax(merge_left, merge_right)
  ) %>%
  filter(
    !is.na(region_start),
    !is.na(region_end),
    region_end > region_start
  )

nrow(s3_clean)


# Convert S3 merged regions to GRanges
s3_regions_gr <- GRanges(
  seqnames = s3_clean$chr,
  ranges = IRanges(
    start = s3_clean$region_start,
    end = s3_clean$region_end
  )
)

mcols(s3_regions_gr) <- s3_clean


# Keep autosomes only
standard_autosomes <- paste0("chr", 1:22)

s3_regions_gr <- keepSeqlevels(
  s3_regions_gr,
  standard_autosomes,
  pruning.mode = "coarse"
)

npc_hg19_clean <- keepSeqlevels(
  npc_hg19,
  standard_autosomes,
  pruning.mode = "coarse"
)

mcclay_hg19_clean <- keepSeqlevels(
  mcclay_hg19,
  standard_autosomes,
  pruning.mode = "coarse"
)

forrest_hg19_clean <- keepSeqlevels(
  forrest_hg19,
  standard_autosomes,
  pruning.mode = "coarse"
)


# Overlap S3 SCZ regions with TCF4 peaks
npc_s3_hits <- findOverlaps(s3_regions_gr, npc_hg19_clean)
mcclay_s3_hits <- findOverlaps(s3_regions_gr, mcclay_hg19_clean)
forrest_s3_hits <- findOverlaps(s3_regions_gr, forrest_hg19_clean)

# Region-level overlap summary
s3_tcf4_overlap_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Total_S3_SCZ_Regions = length(s3_regions_gr),
  SCZ_Regions_Overlapping_TCF4_Peaks = c(
    length(unique(queryHits(npc_s3_hits))),
    length(unique(queryHits(mcclay_s3_hits))),
    length(unique(queryHits(forrest_s3_hits)))
  ),
  TCF4_Peaks_Overlapping_SCZ_Regions = c(
    length(unique(subjectHits(npc_s3_hits))),
    length(unique(subjectHits(mcclay_s3_hits))),
    length(unique(subjectHits(forrest_s3_hits)))
  )
)

s3_tcf4_overlap_summary$Percent_SCZ_Regions_Overlapping_TCF4 <- round(
  s3_tcf4_overlap_summary$SCZ_Regions_Overlapping_TCF4_Peaks /
    s3_tcf4_overlap_summary$Total_S3_SCZ_Regions * 100,
  3
)

s3_tcf4_overlap_summary


# Shared SCZ regions across all three TCF4 datasets
npc_s3_regions <- unique(queryHits(npc_s3_hits))
mcclay_s3_regions <- unique(queryHits(mcclay_s3_hits))
forrest_s3_regions <- unique(queryHits(forrest_s3_hits))

shared_s3_regions <- Reduce(
  intersect,
  list(
    npc_s3_regions,
    mcclay_s3_regions,
    forrest_s3_regions
  )
)

shared_s3_region_summary <- data.frame(
  Total_S3_SCZ_Regions = length(s3_regions_gr),
  Shared_SCZ_Regions_All_Three_TCF4 = length(shared_s3_regions),
  Percent_Shared_SCZ_Regions_All_Three = round(
    length(shared_s3_regions) / length(s3_regions_gr) * 100,
    3
  )
)

shared_s3_region_summary

shared_s3_region_table <- s3_clean[shared_s3_regions, ]

shared_s3_region_table


# Extract regions overlapping each dataset
npc_s3_region_table <- s3_clean[npc_s3_regions, ] %>%
  arrange(top_p)

mcclay_s3_region_table <- s3_clean[mcclay_s3_regions, ] %>%
  arrange(top_p)

forrest_s3_region_table <- s3_clean[forrest_s3_regions, ] %>%
  arrange(top_p)

# Gene-level table from S3
s3_genes <- s3_clean %>%
  select(
    chr,
    region_start,
    region_end,
    top_index,
    top_p,
    ensembl_gene_names
  ) %>%
  separate_rows(ensembl_gene_names, sep = ",|;|\\s+") %>%   #Split genes into separate rows
  mutate(
    gene = str_trim(ensembl_gene_names),                    #Clean whitespace
    gene = str_replace_all(gene, "[^A-Za-z0-9._-]", "")     #Remove weird symbols
  ) %>%
  filter(                 #Remove missing values
    !is.na(gene),
    gene != "",
    gene != "NA"
  ) %>%
  distinct(gene, .keep_all = TRUE)

s3_genes


# Compare S3 genes with TCF4-associated genes
s3_gene_list <- unique(s3_genes$gene)

length(s3_gene_list)

npc_s3_genes <- intersect(npc_genes, s3_gene_list)
mcclay_s3_genes <- intersect(mcclay_genes, s3_gene_list)
forrest_s3_genes <- intersect(forrest_genes, s3_gene_list)

s3_gene_overlap_summary <- data.frame(
  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),
  Total_TCF4_Genes = c(
    length(npc_genes),
    length(mcclay_genes),
    length(forrest_genes)
  ),
  TCF4_Genes_in_S3_SCZ_Regions = c(
    length(npc_s3_genes),
    length(mcclay_s3_genes),
    length(forrest_s3_genes)
  )
)

s3_gene_overlap_summary$Percent_TCF4_Genes_in_S3 <- round(
  s3_gene_overlap_summary$TCF4_Genes_in_S3_SCZ_Regions /
    s3_gene_overlap_summary$Total_TCF4_Genes * 100,
  3
)

s3_gene_overlap_summary


# Shared S3 genes across TCF4 datasets
shared_s3_tcf4_genes <- Reduce(
  intersect,
  list(
    npc_s3_genes,
    mcclay_s3_genes,
    forrest_s3_genes
  )
)

length(shared_s3_tcf4_genes)

shared_s3_tcf4_genes

# PGC3 / Trubetskoy Table S10 vs TCF4 peaks
# ST10b conditional independent SCZ loci vs TCF4 peaks
# Import dataset
st10b_raw <- read_xlsx("NIHMS1824605-supplement-Table_S10.xlsx", sheet = "10b")
colnames(st10b_raw)
dim(st10b_raw)
View(st10b_raw)

# Clean ST10b table
st10b_clean <- st10b_raw %>%
rename(                                            # Clean names based on R best practices of coding
  index_chr = `Index CHR`,                         # Chromosome of the lead/index SNP
  index_snp = `Index SNP`,                         # Lead/index SNP name
  bp_left = `BP Left`,                             # Left/lower boundary of the locus
  bp_right = `BP Right`,                           # Right/upper boundary of the locus
  n_conditional_independent = `# Conditional independent`, # Number of independent signals in this locus

  index_bp = `Index BP`,                           # Base-pair position of lead/index SNP
  index_p = `Index P value`,                       # P-value of lead/index SNP
  merged_indexes = `Merged Indexes`,               # Index SNPs merged into this locus

  second_snp = `2nd SNP`,                          # Secondary independent SNP, if present
  second_bp = `2nd BP`,
  second_p = `2nd P value`,

  third_snp = `3rd SNP`,                           # Third independent SNP, if present
  third_bp = `3rd BP`,
  third_p = `3rd P value`,

  fourth_snp = `4th SNP`,                          # Fourth independent SNP, if present
  fourth_bp = `4th BP`,
  fourth_p = `4th P value`,

  fifth_snp = `5th SNP`,                           # Fifth independent SNP, if present
  fifth_bp = `5th BP`,
  fifth_p = `5th P value`
) %>%
  mutate(                                          # Convert chromosome and numeric columns into usable formats
  chr = paste0("chr", index_chr),                  # Convert chromosome 1 to chr1 format

  bp_left = as.numeric(bp_left),                   # Convert locus boundary to numeric
  bp_right = as.numeric(bp_right),

  index_bp = as.numeric(index_bp),                 # Convert SNP position to numeric
  index_p = as.numeric(index_p),                   # Convert p-value to numeric

  n_conditional_independent =
    as.numeric(n_conditional_independent),         # Convert number of independent signals to numeric
  region_start = pmin(bp_left, bp_right, na.rm = TRUE), # pmin/pmax make sure start is always smaller than end.
  region_end = pmax(bp_left, bp_right, na.rm = TRUE)
) %>%
  filter(                                         # Keep only rows with valid chromosome and genomic region
  !is.na(chr),
  !is.na(region_start),
  !is.na(region_end),
  region_end > region_start
)


nrow(st10b_clean)

# Convert ST10b loci to GRanges
st10b_loci_gr <- GRanges(
  seqnames = st10b_clean$chr,
  ranges = IRanges(
    start = st10b_clean$region_start,
    end = st10b_clean$region_end
  )
)


mcols(st10b_loci_gr) <- st10b_clean

# Keep standard autosomes only
standard_autosomes <- paste0("chr", 1:22)

st10b_loci_gr <- keepSeqlevels(
  st10b_loci_gr,
  standard_autosomes,
  pruning.mode = "coarse"
)

npc_hg19_clean <- keepSeqlevels(
  npc_hg19,
  standard_autosomes,
  pruning.mode = "coarse"
)

mcclay_hg19_clean <- keepSeqlevels(
  mcclay_hg19,
  standard_autosomes,
  pruning.mode = "coarse"
)

forrest_hg19_clean <- keepSeqlevels(
  forrest_hg19,
  standard_autosomes,
  pruning.mode = "coarse"
)

length(st10b_loci_gr)

# Locus-level overlap with TCF4 peak
npc_st10b_hits <- findOverlaps(st10b_loci_gr, npc_hg19_clean)
mcclay_st10b_hits <- findOverlaps(st10b_loci_gr, mcclay_hg19_clean)
forrest_st10b_hits <- findOverlaps(st10b_loci_gr, forrest_hg19_clean)

st10b_tcf4_overlap_summary <- data.frame(
  Dataset = c("NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"),
  Total_ST10b_Loci = length(st10b_loci_gr),
  ST10b_Loci_Overlapping_TCF4_Peaks = c(
    length(unique(queryHits(npc_st10b_hits))),
    length(unique(queryHits(mcclay_st10b_hits))),
    length(unique(queryHits(forrest_st10b_hits)))
  ),
  TCF4_Peaks_Overlapping_ST10b_Loci = c(
    length(unique(subjectHits(npc_st10b_hits))),
    length(unique(subjectHits(mcclay_st10b_hits))),
    length(unique(subjectHits(forrest_st10b_hits)))
  )
)

st10b_tcf4_overlap_summary$Percent_ST10b_Loci_Overlapping_TCF4 <- round(
  st10b_tcf4_overlap_summary$ST10b_Loci_Overlapping_TCF4_Peaks /
    st10b_tcf4_overlap_summary$Total_ST10b_Loci * 100,
  3
)

st10b_tcf4_overlap_summary


# Shared ST10b loci across all three TCF4 datasets
npc_st10b_loci <- unique(queryHits(npc_st10b_hits))
mcclay_st10b_loci <- unique(queryHits(mcclay_st10b_hits))
forrest_st10b_loci <- unique(queryHits(forrest_st10b_hits))

shared_st10b_loci <- Reduce(
  intersect,
  list(
    npc_st10b_loci,
    mcclay_st10b_loci,
    forrest_st10b_loci
  )
)

shared_st10b_locus_summary <- data.frame(
  Total_ST10b_Loci = length(st10b_loci_gr),
  Shared_ST10b_Loci_All_Three_TCF4 = length(shared_st10b_loci),
  Percent_Shared_ST10b_Loci_All_Three = round(
    length(shared_st10b_loci) / length(st10b_loci_gr) * 100,
    3
  )
)

shared_st10b_locus_summary

shared_st10b_locus_table <- st10b_clean[shared_st10b_loci, ] %>%
  arrange(index_p)

shared_st10b_locus_table

# Analyze independent signal SNPs from index, 2nd, 3rd, 4th, and 5th SNP columns
# Build independent-signal SNP table
st10b_signals <- bind_rows(
  st10b_clean %>%
    transmute(
      locus_id = row_number(),
      chr,
      signal_order = "Index",
      snp = index_snp,
      bp = index_bp,
      p_value = index_p
    ),
  st10b_clean %>%
    transmute(
      locus_id = row_number(),
      chr,
      signal_order = "Second",
      snp = second_snp,
      bp = as.numeric(second_bp),
      p_value = as.numeric(second_p)
    ),
  st10b_clean %>%
    transmute(
      locus_id = row_number(),
      chr,
      signal_order = "Third",
      snp = third_snp,
      bp = as.numeric(third_bp),
      p_value = as.numeric(third_p)
    ),
  st10b_clean %>%
    transmute(
      locus_id = row_number(),
      chr,
      signal_order = "Fourth",
      snp = fourth_snp,
      bp = as.numeric(fourth_bp),
      p_value = as.numeric(fourth_p)
    ),
  st10b_clean %>%
    transmute(
      locus_id = row_number(),
      chr,
      signal_order = "Fifth",
      snp = fifth_snp,
      bp = as.numeric(fifth_bp),
      p_value = as.numeric(fifth_p)
    )
) %>%
  filter(
    !is.na(snp),
    snp != "",
    !is.na(bp),
    !is.na(p_value)
  )

nrow(st10b_signals)
head(st10b_signals)


# Convert independent signals to GRanges
st10b_signal_gr <- GRanges(
  seqnames = st10b_signals$chr,
  ranges = IRanges(
    start = st10b_signals$bp,
    end = st10b_signals$bp
  )
)

mcols(st10b_signal_gr) <- st10b_signals

st10b_signal_gr <- keepSeqlevels(
  st10b_signal_gr,
  standard_autosomes,
  pruning.mode = "coarse"
)

length(st10b_signal_gr)

# Exact independent-signal SNP overlap with TCF4 peaks
npc_signal_hits <- findOverlaps(st10b_signal_gr, npc_hg19_clean)
mcclay_signal_hits <- findOverlaps(st10b_signal_gr, mcclay_hg19_clean)
forrest_signal_hits <- findOverlaps(st10b_signal_gr, forrest_hg19_clean)

st10b_signal_overlap_summary <- data.frame(
  Dataset = c("NPC_TCF4", "McClay_TCF4", "Forrest_TCF4"),
  Total_ST10b_Independent_Signals = length(st10b_signal_gr),
  Independent_Signals_Inside_TCF4_Peaks = c(
    length(unique(queryHits(npc_signal_hits))),
    length(unique(queryHits(mcclay_signal_hits))),
    length(unique(queryHits(forrest_signal_hits)))
  ),
  TCF4_Peaks_With_Independent_Signals = c(
    length(unique(subjectHits(npc_signal_hits))),
    length(unique(subjectHits(mcclay_signal_hits))),
    length(unique(subjectHits(forrest_signal_hits)))
  )
)

st10b_signal_overlap_summary$Percent_Signals_Inside_TCF4_Peaks <- round(
  st10b_signal_overlap_summary$Independent_Signals_Inside_TCF4_Peaks /
    st10b_signal_overlap_summary$Total_ST10b_Independent_Signals * 100,
  3
)

st10b_signal_overlap_summary
st10b_tcf4_overlap_summary
shared_st10b_locus_summary
st10b_signal_overlap_summary

####################################################################
# Gene level analysis using S12, S13, S14 gene sets vs TCF datasets
####################################################################
# Import data
s12_raw <- read_xlsx("NIHMS1824605-supplement-Table_S12.xlsx", sheet = "Prioritised")
s12_allcriteria_raw <- read_xlsx("NIHMS1824605-supplement-Table_S12.xlsx", sheet = "ST12 all criteria")
s13_raw <-read_xlsx("NIHMS1824605-supplement-Table_S13.xlsx", sheet = "ST13 Prioritised NS+UTR pp≥.1")
s14_raw <-read_xlsx("NIHMS1824605-supplement-Table_S14.xlsx", sheet = "ST14 Prioritised Single Gene")

colnames(s12_raw)
colnames(s13_raw)
colnames(s14_raw)
colnames(s12_allcriteria_raw)

length(s12_raw$Symbol.ID)
length(s13_raw$gene_symbol)
length(unique(s14_raw$gene_symbol))
length(unique(s12_allcriteria_raw$Symbol.ID))

# Clean S12, S13, and S14 datasets
# S12 - SMR / Hi-C prioritized genes
s12_genes <- unique(
  na.omit(
    trimws(s12_raw$Symbol.ID)
  )
)

length(s12_genes)
head(s12_genes)

# S13 - Coding / UTR prioritized genes
s13_genes <- unique(
  na.omit(
    trimws(s13_raw$gene_symbol)
  )
)

length(s13_genes)
head(s13_genes)

# S14 - Single-gene credible loci
s14_genes <- unique(
  na.omit(
    trimws(s14_raw$gene_symbol)
  )
)

length(s14_genes)
head(s14_genes)

# Calculate the overlaps
# NPC overlaps
npc_s12 <- intersect(
  npc_genes,
  s12_genes
)

npc_s13 <- intersect(
  npc_genes,
  s13_genes
)

npc_s14 <- intersect(
  npc_genes,
  s14_genes
)


# McClay overlaps
mcclay_s12 <- intersect(
  mcclay_genes,
  s12_genes
)

mcclay_s13 <- intersect(
  mcclay_genes,
  s13_genes
)

mcclay_s14 <- intersect(
  mcclay_genes,
  s14_genes
)


# Forrest overlaps
forrest_s12 <- intersect(
  forrest_genes,
  s12_genes
)

forrest_s13 <- intersect(
  forrest_genes,
  s13_genes
)

forrest_s14 <- intersect(
  forrest_genes,
  s14_genes
)


# Summary of prioritized gene overlaps
gene_priority_summary <- data.frame(

  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),

  Total_TCF4_Genes = c(
    length(npc_genes),
    length(mcclay_genes),
    length(forrest_genes)
  ),

  S12_Prioritized_Genes = c(
    length(npc_s12),
    length(mcclay_s12),
    length(forrest_s12)
  ),

  S13_Coding_UTR_Genes = c(
    length(npc_s13),
    length(mcclay_s13),
    length(forrest_s13)
  ),

  S14_Single_Gene_Loci = c(
    length(npc_s14),
    length(mcclay_s14),
    length(forrest_s14)
  )
)

gene_priority_summary


# Gene-level overlap summary with reference set sizes
gene_priority_summary <- data.frame(

  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),

  Total_TCF4_Genes = c(
    length(npc_genes),
    length(mcclay_genes),
    length(forrest_genes)
  ),

  S12_Reference_Genes = length(s12_genes),

  S12_Overlap = c(
    length(npc_s12),
    length(mcclay_s12),
    length(forrest_s12)
  ),

  S13_Reference_Genes = length(s13_genes),

  S13_Overlap = c(
    length(npc_s13),
    length(mcclay_s13),
    length(forrest_s13)
  ),

  S14_Reference_Genes = length(s14_genes),

  S14_Overlap = c(
    length(npc_s14),
    length(mcclay_s14),
    length(forrest_s14)
  )
)

## Calculate percentages relative to each schizophrenia gene set

gene_priority_summary$Percent_S12 <- round(
  gene_priority_summary$S12_Overlap /
    gene_priority_summary$S12_Reference_Genes * 100,
  3
)

gene_priority_summary$Percent_S13 <- round(
  gene_priority_summary$S13_Overlap /
    gene_priority_summary$S13_Reference_Genes * 100,
  3
)

gene_priority_summary$Percent_S14 <- round(
  gene_priority_summary$S14_Overlap /
    gene_priority_summary$S14_Reference_Genes * 100,
  3
)

gene_priority_summary

# Let's do the same with all Criterea on S12
# Extract unique ST12 genes
st12_genes <- unique(
  na.omit(
    trimws(
      s12_allcriteria_raw$Symbol.ID
    )
  )
)

length(st12_genes)


# NPC
npc_st12 <- intersect(
  npc_genes,
  st12_genes
)

# McClay
mcclay_st12 <- intersect(
  mcclay_genes,
  st12_genes
)


# Forrest
forrest_st12 <- intersect(
  forrest_genes,
  st12_genes
)

# ST12 overlap summary
st12_overlap_summary <- data.frame(

  Dataset = c(
    "NPC_TCF4",
    "McClay_TCF4",
    "Forrest_TCF4"
  ),

  Total_TCF4_Genes = length(st12_genes),

  ST12_Genes = c(
    length(npc_st12),
    length(mcclay_st12),
    length(forrest_st12)
  )
)

st12_overlap_summary$Percent_ST12 <- round(
  st12_overlap_summary$ST12_Genes /
    st12_overlap_summary$Total_TCF4_Genes * 100,
  3
)

st12_overlap_summary


# Shared ST12 genes
shared_st12 <- Reduce(
  intersect,
  list(
    npc_st12,
    mcclay_st12,
    forrest_st12
  )
)

length(shared_st12)

shared_st12


# Evidence matrix for shared genes
shared_st12_evidence <- s12_allcriteria_raw %>%
  filter(
    Symbol.ID %in% shared_st12
  )

shared_st12_evidence
