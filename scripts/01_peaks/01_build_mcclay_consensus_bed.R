# Convert the 11,322 high confidence Dr. McClay's peaks to standartt BED file for doing dowstream analysis
# Load libraries
library(tidyverse)
library(tibble)
library(stringr)

# Import dataset
mcclay_overla_peaks <- read.csv("Supplementary Table S1. All 11,322 matching binding sites.csv")

head(mcclay_overla_peaks)
colnames(mcclay_overla_peaks)
view(mcclay_overla_peaks)

# Generate a clean and standard dataset
# We need to define a new start and end points for the peaks
mcclay_consensus <- mcclay_overla_peaks %>%
  mutate(
    chrom = chr.K12,
    start_hg19 = pmin(start.K12, start.N16),
    end_hg19 = pmax(end.K12, end.N16),
    peak_id = paste0("McClay_TCF4_peak_", row_number())
  ) %>%
  select(
    peak_id,
    chrom,
    start_hg19,
    end_hg19,
    score.K12,
    logQ.K12,
    best.score.N16,
    logQ.N16,
    mean.score,
    Ebox.K12,
    Ebox.N16
  )

# Save full annotation table
write_csv(
  mcclay_consensus,
  "McClay_TCF4_11322_consensus_hg19_annotated.csv"
)

# Save BED file for liftOver
mcclay_bed <- mcclay_consensus %>%
  select(
    chrom,
    start_hg19,
    end_hg19,
    peak_id
  )

write_tsv(
  mcclay_bed,
  "McClay_TCF4_11322_consensus_hg19.bed",
  col_names = FALSE
)