#In this script we annotated our full TCF4 ChIP peaks,
#identified the subset of E-box-containing peaks,
#tested whether E-boxes are enriched in promoters,
#defined promoter + E-box genes as candidate direct TCF4 targets,
#and ran pathway enrichment on those genes.


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
    "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "org.Hs.eg.db"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("enrichplot")


#Load packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(base)
library(clusterProfiler)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tibble)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Import data
peakfile <- "NPC_ab21_idr_500bp_summit.bed"


#Peak annotation
peak <- readPeakFile(peakfile)
peakAnno <- annotatePeak(
  peak,
  tssRegion = c(-1000, 1000),
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
  annoDb = "org.Hs.eg.db"
)

anno_df <- as.data.frame(peakAnno)  #full annotation table for all 3947 peaks

if (FALSE) {
#Save output
write.csv(
  anno_df,
  file="NPCs_peak_annotation.csv",
  row.names=FALSE
)
}

#Save plots
pdf("NPC_ab21_peak_annotation_plots.pdf")
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno)
dev.off()

anno_df <- as.data.frame(peakAnno)
#write.csv(anno_df, file="NPC_ab21_EBOX_peak_annotation_hg19.csv", row.names=FALSE)

#Make a reference of E-box motifs
ebox_dir <- "Ebox_predicted"
ebox_list <- list()

# chromosomes 1-22
for (i in 1:22) {

  file <- paste0(ebox_dir, "/CANNTG_", i)  #For each chromosome, read the file and
                                           #put it in the "file"

  if (file.exists(file)) {

    x <- readLines(file, warn = FALSE)  #read the file:used readLines() because those
                                        #files are basically just lines of numbers, one position per line.
    x <- gsub("\r", "", x)              #remove hidden Windows characters
    x <- x[nchar(x) > 0]                #remove empty lines:This keeps only non-empty lines
    pos <- as.integer(x)                #convert positions to integers

    df <- data.frame(                   #Build a proper genomic interval table:This is because
                                        #our files give the motif start position, and E-box length is 6 bp
      chr = paste0("chr", i),
      start = pos,
      end = pos + 5
    )

    ebox_list[[length(ebox_list)+1]] <- df  #store each chromosome table in the list
  }
}

# chrX
file <- paste0(ebox_dir, "/CANNTG_X")
if (file.exists(file)) {
  x <- readLines(file, warn = FALSE)
  x <- gsub("\r", "", x)
  pos <- as.integer(x)

  df <- data.frame(
    chr = "chrX",
    start = pos,
    end = pos + 5
  )

  ebox_list[[length(ebox_list)+1]] <- df
}

# chrY
file <- paste0(ebox_dir, "/CANNTG_Y")
if (file.exists(file)) {
  x <- readLines(file, warn = FALSE)
  x <- gsub("\r", "", x)
  pos <- as.integer(x)

  df <- data.frame(
    chr = "chrY",
    start = pos,
    end = pos + 5
  )

  ebox_list[[length(ebox_list)+1]] <- df
}

ebox_df <- do.call(rbind, ebox_list)

head(ebox_df)
nrow(ebox_df)

if (FALSE) {
# Save file
write.table(
  ebox_df,
  file = "EBOX_all_reference.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
}
#Change the our peak file into a dataframe
peaks_df <- read.table(peakfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#Add column names to peak dataframe
colnames(peaks_df)[1:4] <- c("chr", "start", "end", "peak_id")

#Convert peaks and E-boxes to genomic ranges
peaks_gr <- GRanges(
  seqnames = peaks_df$chr,
  ranges = IRanges(start = peaks_df$start, end = peaks_df$end)
)

ebox_gr <- GRanges(
  seqnames = ebox_df$chr,
  ranges = IRanges(start = ebox_df$start, end = ebox_df$end)
)

#Find overlaps between peaks and E-box motifs
hits <- findOverlaps(peaks_gr, ebox_gr) #Does any peak region overlap any E box region?


# Keep only unique peaks that overlap at least one E-box
peak_hit_index <- unique(queryHits(hits)) #Extracts the row numbers of peaks that overlap E box motifs.
length(peak_hit_index)

peaks_with_ebox <- peaks_df[peak_hit_index, ] #Extract position of peaks that contain E-box

if (FALSE) {
write.table(
  peaks_with_ebox,
  file = "peaks_with_ebox_3040.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
}
#### PART2:Annotate the E-box peaks specifically ########
#Load data
peak2_file <- "peaks_with_ebox_3040.bed"
peak2 <- readPeakFile(peak2_file)

peakAnno2 <- annotatePeak(
  peak2,
  tssRegion = c(-1000, 1000),
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
  annoDb = "org.Hs.eg.db"
)

anno_df2 <- as.data.frame(peakAnno2)

#write.csv(anno_df2, file="NPC_EBOX_peak_annotation_3040.csv", row.names=FALSE)

#### PART 3: Are E-box motifs more common in promoter peaks or in non promoter peaks?
View(anno_df) #full annotation table for all 3947 peaks


#Create binary promoter column (1 = promoter peak , 0 = non promoter peak)
anno_df$region_type <- ifelse(grepl("Promoter", anno_df$annotation),1,0)
table(anno_df$region_type)

#Add E-box binary column
#Load e-box peak
ebox <- read.table("peaks_with_ebox_3040.bed", header=FALSE)
colnames(ebox) <- c("chr","start","end","peak_id")
View(ebox)

#Mark which peaks contain e-box(1 = peak contains E-box, 0 = peak does not contain E-box)
anno_df$ebox_binary <- ifelse(anno_df$V4 %in% ebox$peak_id, 1, 0)

table(anno_df$ebox_binary)

#Build a 2x2 table
table(anno_df$region_type, anno_df$ebox_binary)
prop.table(table(anno_df$region_type, anno_df$ebox_binary), margin = 1)

#Fisher test to see if this enrichment statistically significant?
mat <- table(anno_df$region_type, anno_df$ebox_binary)
fisher.test(mat)

#Data integration and compare our TCF4 gene sets with external gene sets
#Extract all peak genes
gene_all <- unique(anno_df$SYMBOL)
length(gene_all)

if (FALSE) {
write.csv(
  data.frame(Gene = gene_all),
  file = "All_TCF4_peak_genes.csv",
  row.names = FALSE
)
}
#Extract E-box peak genes
gene_motif <- unique(anno_df2$SYMBOL)
length(gene_motif)

if (FALSE) {
write.csv(
  data.frame(Gene = gene_motif),
  file = "EBOX_peak_genes.csv",
  row.names = FALSE
)
}
#Extract Promoter + E-box genes
promoter_ebox_genes <- unique(
  anno_df$SYMBOL[
    anno_df$ebox_binary == 1 &
      anno_df$annotation == "Promoter"])

length(promoter_ebox_genes)

if (FALSE) {
write.csv(
  data.frame(Gene = promoter_ebox_genes),
  file = "Promoter_EBOX_genes.csv",
  row.names = FALSE
)
}
##### MAKE A NEW REFERENCE #######
#Load data
#Background gene list
#New UCSC style reference, transcript level
new_refGene <- read.csv("ref_Gene_2026.csv", stringsAsFactors = FALSE)

#Dr. McClay old gene-level reference
JLM_Ref_gene <- read.csv("Data integration/refGene_autoX_geneposns_redux.csv",
                         stringsAsFactors = FALSE)

#Check column names
colnames(new_refGene)
colnames(JLM_Ref_gene)

#Count number of genes in new_refGene and copmare with Dr.McClay's refGene list
length(unique(new_refGene$name2))
length(unique(JLM_Ref_gene$gene))

#Keep only standard chromosomes in the new reference
standard_chr <- paste0("chr", c(1:22, "X", "Y"))  #keep only genes on standard chromosomes (chr1-chr22, chrX, chrY)

new_refGene_std <- new_refGene %>%
  filter(chrom %in% standard_chr)

#Remove genes that map to more than one standard chromosome
#This avoids ambiguous gene-level collapsing
multi_chr_genes <- new_refGene_std %>%
  group_by(name2) %>%
  summarise(n_chr = n_distinct(chrom), .groups = "drop") %>%
  filter(n_chr > 1) %>%
  pull(name2)

new_refGene_std2 <- new_refGene_std %>%
  filter(!name2 %in% multi_chr_genes)

#Convert transcript-level new reference to gene-level(one row per gene)
new_refGene_genelevel <- new_refGene_std2 %>%
  group_by(name2) %>%
  summarise(
    chrom =  dplyr::first(chrom),                       #takes the chromosome from the first transcript
    txStart_min = min(txStart, na.rm = TRUE),   #gets the earliest transcript start position
    txEnd_max   = max(txEnd, na.rm = TRUE),     #gets the latest transcript end position
    n_transcripts = n(),                        #counts how many transcript rows belong to that gene
    .groups = "drop"
  ) %>%
  rename(gene = name2) %>%
  mutate(
    chr_num = case_when(        #Convert chromosome names to numeric codes for gen.start / gen.end like chr1.
      chrom == "chrX" ~ 23,
      chrom == "chrY" ~ 24,
      TRUE ~ as.numeric(str_remove(chrom, "^chr"))  #For chr1-chr22, remove "chr" prefix and convert to numeric
    ),
    gen.start = chr_num * 1e9 + txStart_min,     #combine chromosome + position into one number for gen.start
    gen.end   = chr_num * 1e9 + txEnd_max        #combine chromosome + position into one number for gen.end
  ) %>%
  select(gene, chrom, txStart_min, txEnd_max, n_transcripts, gen.start, gen.end)

#summaries
cat("Unique genes in new transcript-level reference:",
    length(unique(new_refGene$name2)), "\n")

cat("Unique genes in cleaned new gene-level reference:",
    length(unique(new_refGene_genelevel$gene)), "\n")

cat("Unique genes in McClay old reference:",
    length(unique(JLM_Ref_gene$gene)), "\n")

#compare old and new gene symbols directly
genes_new <- unique(new_refGene_genelevel$gene)
genes_old <- unique(JLM_Ref_gene$gene)

genes_new_only <- setdiff(genes_new, genes_old)     #genes present only in the new reference
genes_old_only <- setdiff(genes_old, genes_new)     #genes present only in the old reference
genes_overlap  <- intersect(genes_new, genes_old)

#First pass symbol conversion with limma (convert outdated names into current approved symbols)
alias_map_limma <- tibble(old_gene = genes_old_only) %>%
  mutate(mapped_symbol_limma = alias2Symbol(old_gene, species = "Hs"))

#the above code is not working because alias2Symbol() does not always return one output for each input gene.
#we need to fix it with a new loop to do it one by one (gene by gene)
mapped_symbol_limma <- rep(NA_character_, length(genes_old_only))

for (i in seq_along(genes_old_only)) {
  out <- alias2Symbol(genes_old_only[i], species = "Hs")
  if (length(out) > 0) {
    mapped_symbol_limma[i] <- out[1]
  }
}

alias_map_limma <- tibble(
  old_gene = genes_old_only,
  mapped_symbol_limma = mapped_symbol_limma)

#sanity check
table(is.na(alias_map_limma$mapped_symbol_limma))
head(alias_map_limma)

#reveal the 838 unresolved genes (I checked with HGNC: all of them still unknown)
still_unmapped_after_limma <- alias_map_limma %>%
  filter(is.na(mapped_symbol_limma) | mapped_symbol_limma == "") %>%
  pull(old_gene)

unresolved_genes <- alias_map_limma %>%
  filter(is.na(mapped_symbol_limma)) %>%
  distinct(old_gene)

nrow(unresolved_genes)
head(unresolved_genes)

#Second pass symbol conversion with org.Hs.eg.db (another lookup using the Bioconductor human annotation database)
alias_map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = still_unmapped_after_limma,
  columns = c("SYMBOL", "ALIAS"),
  keytype = "ALIAS"
)

#### it gaves me an error:Among the genes you gave me, zero are valid ALIAS keys in this database.
#Build a clean conversion table from successful mappings only
conversion_table <- alias_map_limma %>%
  filter(!is.na(mapped_symbol_limma) & mapped_symbol_limma != "") %>%
  rename(mapped_symbol = mapped_symbol_limma) %>%
  distinct(old_gene, .keep_all = TRUE)

nrow(conversion_table)
head(conversion_table)

#Check whether the mapped current symbols are already in the new reference
conversion_table <- alias_map_limma %>%
  filter(!is.na(mapped_symbol_limma) & mapped_symbol_limma != "") %>%  #Keep only successful mappings
  rename(mapped_symbol = mapped_symbol_limma) %>%  #Rename the column for clarity
  distinct(old_gene, .keep_all = TRUE)             #Remove any duplicate old_gene entries

nrow(conversion_table)
head(conversion_table)

#Check whether the mapped current symbols are already in the new reference
conversion_table <- conversion_table %>%
  mutate(
    mapped_symbol_in_new_ref = mapped_symbol %in% new_refGene_genelevel$gene
  )

table(conversion_table$mapped_symbol_in_new_ref)
####result: FALSE  TRUE
#           120    1096
#1096 genes (TRUE):These mapped symbols already exist in our new 2026 reference
#120 genes (FALSE):These mapped symbols do NOT exist in our new reference (potentially missing genes)
#Keep a table of renamed genes for documentation
renamed_genes_table <- conversion_table %>%
  filter(old_gene != mapped_symbol) %>%
  arrange(old_gene)

head(renamed_genes_table)
nrow(renamed_genes_table)
View(renamed_genes_table)

#Make a summary table
conversion_summary <- conversion_table %>%
  summarise(
    n_successfully_mapped = n(),
    n_already_in_new_ref = sum(mapped_symbol_in_new_ref),
    n_not_in_new_ref = sum(!mapped_symbol_in_new_ref),
    n_renamed = sum(old_gene != mapped_symbol),
    n_same_symbol = sum(old_gene == mapped_symbol)
  )

conversion_summary

#Pull unresolved genes from the old reference
old_unresolved_rows <- JLM_Ref_gene %>%
  filter(gene %in% unresolved_genes$old_gene) %>%
  distinct(gene, .keep_all = TRUE)

nrow(old_unresolved_rows)
head(old_unresolved_rows)

#Decode old coordinates into chromosome and position (for old unresolved genes)
old_unresolved_ranges <- old_unresolved_rows %>%
  mutate(
    chr_num_old = floor(gen.start / 1e9),
    pos_start_old = gen.start - chr_num_old * 1e9,
    pos_end_old   = gen.end - chr_num_old * 1e9,
    chrom = case_when(
      chr_num_old == 23 ~ "chrX",
      chr_num_old == 24 ~ "chrY",
      TRUE ~ paste0("chr", chr_num_old)
    )
  ) %>%
  select(gene, chrom, pos_start_old, pos_end_old, gen.start, gen.end)

head(old_unresolved_ranges)


#Compare unresolved old genes to the new reference by coordinate overlap
#we are ask if an unresolved old gene overlap a gene in the new reference on the same chromosome.
coordinate_matches <- old_unresolved_ranges %>%
  inner_join(
    new_refGene_genelevel %>%
      select(gene_new = gene, chrom, txStart_min, txEnd_max),
    by = "chrom"
  ) %>%
  filter(
    pos_start_old <= txEnd_max,
    pos_end_old   >= txStart_min
  ) %>%
  mutate(
    overlap_start = pmax(pos_start_old, txStart_min),
    overlap_end   = pmin(pos_end_old, txEnd_max),
    overlap_bp    = overlap_end - overlap_start + 1,
    old_width     = pos_end_old - pos_start_old + 1,
    new_width     = txEnd_max - txStart_min + 1,
    prop_old_overlap = overlap_bp / old_width,
    prop_new_overlap = overlap_bp / new_width
  ) %>%
  filter(overlap_bp > 0) %>%
  arrange(gene, desc(overlap_bp))

head(coordinate_matches)

#The new genes in the old reference that do not have any coordinate matches in the new reference are
#likely truly missing genes
#we don't need to add to our final reference gene list
#Because almost all of  these 120 genes are antisense, microRNA, small RNAs,pseudogenes, and etc.
conversion_table %>%
  filter(mapped_symbol_in_new_ref == FALSE) %>%
  pull(old_gene)

#Before making a clean gene list, I want to count the exclusion genes
filter_summary <- new_refGene_genelevel %>%
  summarise(
    total_genes = n(),
    n_MIR = sum(str_detect(gene, "^MIR[0-9]")),
    n_SNORA = sum(str_detect(gene, "^SNORA")),
    n_SNORD = sum(str_detect(gene, "^SNORD")),
    n_SNOR  = sum(str_detect(gene, "^SNOR")),
    n_RNU   = sum(str_detect(gene, "^RNU")),
    n_RNA5S = sum(str_detect(gene, "^RNA5S")),
    n_RNA5  = sum(str_detect(gene, "^RNA5-")),
    n_TRNA  = sum(str_detect(gene, "^TRNA")),
    n_antisense = sum(str_detect(gene, "-AS[0-9]*$")),
    n_divergent = sum(str_detect(gene, "-DT$")),
    n_readthrough = sum(str_detect(str_to_lower(gene), "readthrough")),
    n_pseudogene_number = sum(str_detect(gene, "P[0-9]+$"))
  )

View(filter_summary)

final_ref_gene <- new_refGene_genelevel

if (FALSE) {
#make final reference gene list
final_ref_gene <- new_refGene_genelevel %>%
  filter(
    !str_detect(gene, "^MIR[0-9]"),
    !str_detect(gene, "^MIRLET"),
    # Exclude small RNA classes
    !str_detect(gene, "^SNORA"),
    !str_detect(gene, "^SNORD"),
    !str_detect(gene, "^SNOR"),
    !str_detect(gene, "^RNU"),
    !str_detect(gene, "^RNA5S"),
    !str_detect(gene, "^RNA5-"),
    !str_detect(gene, "^TRNA"),
    # Exclude antisense transcripts
    !str_detect(gene, "-AS[0-9]*$"),
    # Exclude divergent transcripts
    !str_detect(gene, "-DT$"),
    # Exclude readthrough transcripts
    !str_detect(str_to_lower(gene), "readthrough"),
    # Exclude likely pseudogenes
    !str_detect(gene, "P[0-9]+$")
    ) %>%
  distinct(gene, .keep_all = TRUE)

}


#Save the filtered final reference gene list
write.csv(final_ref_gene,
          "final_ref_gene_2026_filtered.csv",
          row.names = FALSE)


#Extract the unique gene symbols from the cleaned final reference gene list to create a
#comprehensive list of all genes in the genome for our analysis.
allgenes <- final_ref_gene$gene
allgenes <- unique(allgenes[!is.na(allgenes) & allgenes != ""])
length(allgenes)

#Check again. (remove NAs + empty values + duplicates)
gene_all <- unique(gene_all[!is.na(gene_all) & gene_all != ""])
gene_motif <- unique(gene_motif[!is.na(gene_motif) & gene_motif != ""])
promoter_ebox_genes <- unique(promoter_ebox_genes[!is.na(promoter_ebox_genes) & promoter_ebox_genes != ""])

length(gene_all)
length(gene_motif)
length(promoter_ebox_genes)

#Create logical vectors for our TCF4 gene sets
#For each gene in the genome, check if it’s in our TCF4 list.
is_gene_all <- allgenes %in% gene_all
is_gene_motif <- allgenes %in% gene_motif  #Which genes have TCF4 + E-box
is_promoter_ebox <- allgenes %in% promoter_ebox_genes  #Which genes have TCF4 + E-box + promoter

sum(is_gene_all)
sum(is_gene_motif)
sum(is_promoter_ebox)

#Load external datasets
########################
# Blake TCF4 knockdown #
########################

#Load genes differentially expressed in Blake TCF4 siRNA paper in SH-SY5Y
blake_all <- readLines("Data integration/TCF4 KD - 1205 DEG.txt")
blake_up <- readLines("Data integration/TCF4 blake up n470.txt")
blake_down <- readLines("Data integration/TCF4 blake down n665.txt")

#Clean data
blake_all <- unique(blake_all[blake_all != ""])
blake_up <- unique(blake_up[blake_up != ""])
blake_down <- unique(blake_down[blake_down != ""])

#Create logic vector and count how many genes matched
is_blake_all <- allgenes %in% blake_all
is_blake_up <- allgenes %in% blake_up
is_blake_down <- allgenes %in% blake_down

sum(is_blake_all)
sum(is_blake_up)
sum(is_blake_down)

#Fisher tests
#ALL peak genes
# All peaks vs Blake
fisher_all_blake <- fisher.test(table(is_gene_all, is_blake_all), alternative = "greater")
fisher_all_blake

#Extract the overlapped genes (all ChIP genes and Blake)
overlap_all_blake <- allgenes[is_gene_all & is_blake_all]
overlap_all_blake_up <-  allgenes[is_gene_all & is_blake_up]
overlap_all_blake_down <-  allgenes[is_gene_all & is_blake_down]

# All peaks vs Blake up
fisher_all_blake_up <- fisher.test(table(is_gene_all, is_blake_up), alternative = "greater")
fisher_all_blake_up

# All peaks vs Blake down
fisher_all_blake_down <- fisher.test(table(is_gene_all, is_blake_down), alternative = "greater")
fisher_all_blake_down

# EBOX vs Blake
fisher_motif_blake <- fisher.test(table(is_gene_motif, is_blake_all), alternative = "greater")
fisher_motif_blake

# EBOX vs Blake up
fisher_motif_blake_up <- fisher.test(table(is_gene_motif, is_blake_up), alternative = "greater")
fisher_motif_blake_up

# EBOX vs Blake down
fisher_motif_blake_down <- fisher.test(table(is_gene_motif, is_blake_down), alternative = "greater")
fisher_motif_blake_down

# Extract the overlapped genes (Ebox and Blake)
overlap_motif_blake <- allgenes[is_gene_motif & is_blake_all]
overlap_motif_blake_up <-  allgenes[is_gene_motif & is_blake_up]
overlap_motif_blake_down <-  allgenes[is_gene_motif & is_blake_down]

# Promoter EBOX vs Blake
fisher_promoter_blake <- fisher.test(table(is_promoter_ebox,is_blake_all), alternative = "greater")
fisher_promoter_blake

# Promoter EBOX vs Blake up
fisher_promoter_blake_up <- fisher.test(table(is_promoter_ebox,is_blake_up), alternative = "greater")
fisher_promoter_blake_up

# Promoter EBOX vs Blake down
fisher_promoter_blake_down <- fisher.test(table(is_promoter_ebox,is_blake_down), alternative = "greater")
fisher_promoter_blake_down

# Extract the overlapped genes (Promoter Ebox and Blake)
overlap_PromoterEbox_blake <- allgenes[is_promoter_ebox & is_blake_all]
overlap_PromoterEbox_blake_up <-  allgenes[is_promoter_ebox & is_blake_up]
overlap_PromoterEbox_blake_down <-  allgenes[is_promoter_ebox & is_blake_down]

overlap_PromoterEbox_blake
overlap_PromoterEbox_blake
overlap_PromoterEbox_blake_down
########################
# Fromer Schizophrenia #
########################
# Load genes differentially expressed in Fromer et al gene expression study of SCZ
fromgenes = readLines("Data integration/Fromer_genes.txt")
is_fromer_all = (allgenes %in% fromgenes);
sum(is_fromer_all)

fromgenesup = readLines("Data integration/Fromer-up.txt")
is_fromer_up = (allgenes %in% fromgenesup);
sum(is_fromer_up)

fromgenesdown = readLines("Data integration/Fromer-down.txt")
is_fromer_down = (allgenes %in% fromgenesdown);
sum(is_fromer_down)

#All peak genes
# All peaks vs Fromer
fisher_all_fromer <- fisher.test(table(is_gene_all, is_fromer_all), alternative = "greater")
fisher_all_fromer
table(is_gene_all, is_fromer_all)

# All peaks vs Fromer up
fisher_all_fromer_up <- fisher.test(table(is_gene_all, is_fromer_up), alternative = "greater")
fisher_all_fromer_up
table(is_gene_all, is_fromer_up)

# All peaks vs Fromer down
fisher_all_fromer_down <- fisher.test(table(is_gene_all, is_fromer_down), alternative = "greater")
fisher_all_fromer_down
table(is_gene_all, is_fromer_down)

# Extract the overlapped genes
overlap_all_fromer <- allgenes[is_gene_all & is_fromer_all]
overlap_all_fromer_up <-  allgenes[is_gene_all & is_fromer_up]
overlap_all_fromer_down <-  allgenes[is_gene_all & is_fromer_down]
overlap_all_fromer
overlap_all_fromer_up
overlap_all_fromer_down

# EBOX vs Fromer
fisher_motif_fromer <- fisher.test(table(is_gene_motif, is_fromer_all), alternative = "greater")
fisher_motif_fromer
table(is_gene_motif, is_fromer_all)

# EBOX vs Fromer up
fisher_motif_fromer_up <- fisher.test(table(is_gene_motif, is_fromer_up), alternative = "greater")
fisher_motif_fromer_up
table(is_gene_motif, is_fromer_up)

# EBOX vs Fromer down
fisher_motif_fromer_down <- fisher.test(table(is_gene_motif, is_fromer_down), alternative = "greater")
fisher_motif_fromer_down
table(is_gene_motif, is_fromer_down)

#Extract the overlapped genes
overlap_motif_fromer <- allgenes[is_gene_motif & is_fromer_all]
overlap_motif_fromer_up <-  allgenes[is_gene_motif & is_fromer_up]
overlap_motif_fromer_down <-  allgenes[is_gene_motif & is_fromer_down]
overlap_motif_fromer
overlap_motif_fromer_up
overlap_motif_fromer_down

# Promoter EBOX vs Fromer
fisher_promoter_fromer <- fisher.test(table(is_promoter_ebox, is_fromer_all), alternative = "greater")
fisher_promoter_fromer
table(is_promoter_ebox, is_fromer_all)

# Promoter EBOX vs Fromer up
fisher_promoter_fromer_up <- fisher.test(table(is_promoter_ebox, is_fromer_up), alternative = "greater")
fisher_promoter_fromer_up
table(is_promoter_ebox, is_fromer_up)

# Promoter EBOX vs Fromer down
fisher_promoter_fromer_down <- fisher.test(table(is_promoter_ebox, is_fromer_down), alternative = "greater")
fisher_promoter_fromer_down
table(is_promoter_ebox, is_fromer_down)

#Extract the overlapped genes
overlap_PromoterEbox_fromer <- allgenes[is_promoter_ebox & is_fromer_all]
overlap_PromoterEbox_fromer_up <-  allgenes[is_promoter_ebox & is_fromer_up]
overlap_PromoterEbox_fromer_down <-  allgenes[is_promoter_ebox & is_fromer_down]
overlap_PromoterEbox_fromer
overlap_PromoterEbox_fromer_up
overlap_PromoterEbox_fromer_down

# Create comprehensive TCF4 vs Fromer overlap table
# All peaks
fromer_all_up_df <- data.frame(
  gene = overlap_all_fromer_up,
  TCF4_Set = "All peaks",
  Fromer_Direction = "Up"
)

fromer_all_down_df <- data.frame(
  gene = overlap_all_fromer_down,
  TCF4_Set = "All peaks",
  Fromer_Direction = "Down"
)

# Motif peaks
fromer_motif_up_df <- data.frame(
  gene = overlap_motif_fromer_up,
  TCF4_Set = "Motif peaks",
  Fromer_Direction = "Up"
)

fromer_motif_down_df <- data.frame(
  gene = overlap_motif_fromer_down,
  TCF4_Set = "Motif peaks",
  Fromer_Direction = "Down"
)

# Promoter + E-box peaks
fromer_promoter_up_df <- data.frame(
  gene = overlap_PromoterEbox_fromer_up,
  TCF4_Set = "Promoter E-box",
  Fromer_Direction = "Up"
)

fromer_promoter_down_df <- data.frame(
  gene = overlap_PromoterEbox_fromer_down,
  TCF4_Set = "Promoter E-box",
  Fromer_Direction = "Down"
)

# Combine all overlap tables
TCF4_Fromer_overlap_long <- bind_rows(
  fromer_all_up_df,
  fromer_all_down_df,

  fromer_motif_up_df,
  fromer_motif_down_df,

  fromer_promoter_up_df,
  fromer_promoter_down_df
) %>%
  distinct(gene, TCF4_Set, Fromer_Direction) %>%
  arrange(TCF4_Set, Fromer_Direction, gene)

# View final table
View(TCF4_Fromer_overlap_long)

# Number of rows
nrow(TCF4_Fromer_overlap_long)

# Summary table
TCF4_Fromer_overlap_summary <- data.frame(

  TCF4_Set = c(
    "All peaks",
    "Motif peaks",
    "Promoter E-box"
  ),

  Upregulated_Overlap = c(
    length(overlap_all_fromer_up),
    length(overlap_motif_fromer_up),
    length(overlap_PromoterEbox_fromer_up)
  ),

  Downregulated_Overlap = c(
    length(overlap_all_fromer_down),
    length(overlap_motif_fromer_down),
    length(overlap_PromoterEbox_fromer_down)
  )
)

TCF4_Fromer_overlap_summary

# Save table
write.csv(
  TCF4_Fromer_overlap_long,
  "TCF4_Fromer_overlap.csv",
  row.names = FALSE
)
###############################
# SCZ Transcriptome  (Ruzicka)#
###############################
#https://pmc.ncbi.nlm.nih.gov/articles/PMC12772489/
#Load dataset
SCZ_Trans_data <- read_excel("Transcriptome_SCZ_SingleCell.xlsx")

#Clean dataset and select the relevant columns for SCZ single cell data
View(SCZ_Trans_data)
colnames(SCZ_Trans_data)

In_SCZ <- SCZ_Trans_data %>%
  select(...1,In_SZCS) %>%
  filter (In_SZCS != 0 ) %>%
  arrange(desc(In_SZCS))

Ex_SZCS <- SCZ_Trans_data %>%
  select(...1,Ex_SZCS) %>%
  filter (Ex_SZCS != 0 ) %>%
  arrange(desc(Ex_SZCS))

Ex_SZTR <- SCZ_Trans_data %>%
  select(...1,Ex_SZTR) %>%
  filter (Ex_SZTR != 0 ) %>%
  arrange(desc(Ex_SZTR))

is_In_SCZ_all = (allgenes %in% In_SCZ$...1);
sum(is_In_SCZ_all)

is_Ex_SZCS_all = (allgenes %in% Ex_SZCS$...1);
sum(is_Ex_SZCS_all)

is_Ex_SZTR_all = (allgenes %in% Ex_SZTR$...1);
sum(is_Ex_SZTR_all)

#All peak genes(From ChIP-seq)
# All peaks vs In_SCZ
fisher_all_In_SCZ <- fisher.test(table(is_gene_all, is_In_SCZ_all), alternative = "greater")
fisher_all_In_SCZ
table(is_gene_all, is_In_SCZ_all)

# All peaks vs Ex_SZCS
fisher_all_Ex_SZCS <- fisher.test(table(is_gene_all, is_Ex_SZCS_all), alternative = "greater")
fisher_all_Ex_SZCS
table(is_gene_all, is_Ex_SZCS_all)

# All peaks vs Ex_SZTR
fisher_all_Ex_SZTR <- fisher.test(table(is_gene_all, is_Ex_SZTR_all), alternative = "greater")
fisher_all_Ex_SZTR
table(is_gene_all, is_Ex_SZTR_all)

#Extract the overlapped genes
overlap_all_In_SCZ <- allgenes[is_gene_all & is_In_SCZ_all]
overlap_all_Ex_SZCS <- allgenes[is_gene_all & is_Ex_SZCS_all]
overlap_all_Ex_SZTR <- allgenes[is_gene_all & is_Ex_SZTR_all]
overlap_all_In_SCZ
overlap_all_Ex_SZCS
overlap_all_Ex_SZTR

#Prepare my gene lists
genes_in  <- unique(overlap_all_In_SCZ)
genes_ex  <- unique(overlap_all_Ex_SZCS)
genes_tr  <- unique(overlap_all_Ex_SZTR)

# Define the background universe
universe_genes <- unique(allgenes) #all genes linked to all ChIP-seq peaks

# Convert gene symbols to Entrez IDs
gene_map_in <- bitr(genes_in,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db
)

gene_map_ex <- bitr(genes_ex,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db
)

gene_map_tr <- bitr(genes_tr,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db
)

universe_map <- bitr(unique(universe_genes),
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db
)

# Keep only mapped Entrez IDs
entrez_in <- unique(gene_map_in$ENTREZID)
entrez_ex <- unique(gene_map_ex$ENTREZID)
entrez_tr <- unique(gene_map_tr$ENTREZID)
entrez_universe <- unique(universe_map$ENTREZID)

# GO Biological Process enrichment
ego_in <- enrichGO(
  gene = entrez_in,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_ex <- enrichGO(
  gene = entrez_ex,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_tr <- enrichGO(
  gene = entrez_tr,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG pathway enrichment
kegg_in <- enrichKEGG(
  gene = entrez_in,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_ex <- enrichKEGG(
  gene = entrez_ex,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_tr <- enrichKEGG(
  gene = entrez_tr,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

#Reactome pathway enrichment
react_in <- enrichPathway(
  gene = entrez_in,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

react_ex <- enrichPathway(
  gene = entrez_ex,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

react_tr <- enrichPathway(
  gene = entrez_tr,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#Save results as tables
go_in_df <- as.data.frame(ego_in)
go_ex_df <- as.data.frame(ego_ex)
go_tr_df <- as.data.frame(ego_tr)

react_in_df <- as.data.frame(react_in)
react_ex_df <- as.data.frame(react_ex)
react_tr_df <- as.data.frame(react_tr)

kegg_in_df <- as.data.frame(kegg_in)
kegg_ex_df <- as.data.frame(kegg_ex)
kegg_tr_df <- as.data.frame(kegg_tr)

# Plot top enriched pathways
#GO
dotplot(ego_in, showCategory = 15) + ggtitle("GO BP enrichment: overlap_all_In_SCZ")
dotplot(ego_ex, showCategory = 15) + ggtitle("GO BP enrichment: overlap_all_Ex_SZCS")
dotplot(ego_tr, showCategory = 15) + ggtitle("GO BP enrichment: overlap_all_Ex_SZTR")

#Reactome
dotplot(react_in, showCategory = 15) + ggtitle("Reactome enrichment: overlap_all_In_SCZ")
dotplot(react_ex, showCategory = 15) + ggtitle("Reactome enrichment: overlap_all_Ex_SZCS")
dotplot(react_tr, showCategory = 15) + ggtitle("Reactome enrichment: overlap_all_Ex_SZTR")

#KEGG
dotplot(kegg_in, showCategory = 15) + ggtitle("KEGG enrichment: overlap_all_In_SCZ")
dotplot(kegg_ex, showCategory = 15) + ggtitle("KEGG enrichment: overlap_all_Ex_SZCS")
dotplot(kegg_tr, showCategory = 15) + ggtitle("KEGG enrichment")

#Motif peaks vs In_SCZ
fisher_motif_In_SCZ <- fisher.test(table(is_gene_motif, is_In_SCZ_all), alternative = "greater")
fisher_motif_In_SCZ
table(is_gene_motif, is_In_SCZ_all)

# Motif peaks vs Ex_SZCS
fisher_motif_Ex_SZCS <- fisher.test(table(is_gene_motif, is_Ex_SZCS_all), alternative = "greater")
fisher_motif_Ex_SZCS
table(is_gene_motif, is_Ex_SZCS_all)

# Motif peaks vs Ex_SZTR
fisher_motif_Ex_SZTR <- fisher.test(table(is_gene_motif, is_Ex_SZTR_all), alternative = "greater")
fisher_motif_Ex_SZTR
table(is_gene_motif, is_Ex_SZTR_all)

# Find overlap genes
overlap_motif_In_SCZ <- allgenes[is_gene_motif & is_In_SCZ_all]
overlap_motif_Ex_SZCS <- allgenes[is_gene_motif & is_Ex_SZCS_all]
overlap_motif_Ex_SZTR <- allgenes[is_gene_motif & is_Ex_SZTR_all]
overlap_motif_In_SCZ
overlap_motif_Ex_SZCS
overlap_motif_Ex_SZTR

# Convert motif overlap genes to Entrez IDs
gene_map_motif_in <- bitr(
  overlap_motif_In_SCZ,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_map_motif_ex <- bitr(
  overlap_motif_Ex_SZCS,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_map_motif_tr <- bitr(
  overlap_motif_Ex_SZTR,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Extract Entrez IDs
entrez_motif_in <- unique(gene_map_motif_in$ENTREZID)
entrez_motif_ex <- unique(gene_map_motif_ex$ENTREZID)
entrez_motif_tr <- unique(gene_map_motif_tr$ENTREZID)

# GO enrichment
go_motif_in <- enrichGO(
  gene = entrez_motif_in,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  readable = TRUE
)

go_motif_ex <- enrichGO(
  gene = entrez_motif_ex,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  readable = TRUE
)

go_motif_tr <- enrichGO(
  gene = entrez_motif_tr,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  readable = TRUE
)

# KEGG enrichment
kegg_motif_in <- enrichKEGG(
  gene = entrez_motif_in,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_motif_ex <- enrichKEGG(
  gene = entrez_motif_ex,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_motif_tr <- enrichKEGG(
  gene = entrez_motif_tr,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Reactome enrichment for motif overlap genes
react_motif_in <- enrichPathway(
  gene = entrez_motif_in,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

react_motif_ex <- enrichPathway(
  gene = entrez_motif_ex,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

react_motif_tr <- enrichPathway(
  gene = entrez_motif_tr,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot the results
# Reactome plots
dotplot(react_motif_in, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment: Motif Peaks vs In_SCZ")

dotplot(react_motif_ex, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment: Motif Peaks vs Ex_SZCS")

dotplot(react_motif_tr, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment: Motif Peaks vs Ex_SZTR")


# KEGG plots
dotplot(kegg_motif_in, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment: Motif Peaks vs In_SCZ")

dotplot(kegg_motif_ex, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment: Motif Peaks vs Ex_SZCS")

dotplot(kegg_motif_tr, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment: Motif Peaks vs Ex_SZTR")


# GO plots
dotplot(go_motif_in, showCategory = 15) +
  ggtitle("GO Biological Process: Motif Peaks vs In_SCZ")

dotplot(go_motif_ex, showCategory = 15) +
  ggtitle("GO Biological Process: Motif Peaks vs Ex_SZCS")

dotplot(go_motif_tr, showCategory = 15) +
  ggtitle("GO Biological Process: Motif Peaks vs Ex_SZTR")


# Promoter EBOX vs In_SCZ
fisher_promoterEbox_In_SCZ <- fisher.test(
  table(is_promoter_ebox, is_In_SCZ_all),
  alternative = "greater"
)
fisher_promoterEbox_In_SCZ
table(is_promoter_ebox, is_In_SCZ_all)

# Promoter EBOX vs Ex_SZCS
fisher_promoterEbox_Ex_SZCS <- fisher.test(
  table(is_promoter_ebox, is_Ex_SZCS_all),
  alternative = "greater"
)
fisher_promoterEbox_Ex_SZCS
table(is_promoter_ebox, is_Ex_SZCS_all)

# Promoter EBOX vs Ex_SZTR
fisher_promoterEbox_Ex_SZTR <- fisher.test(
  table(is_promoter_ebox, is_Ex_SZTR_all),
  alternative = "greater"
)
fisher_promoterEbox_Ex_SZTR
table(is_promoter_ebox, is_Ex_SZTR_all)

# Extract the overlap genes
overlap_PromoterEbox_In_SCZ <- unique(allgenes[is_promoter_ebox & is_In_SCZ_all])
overlap_PromoterEbox_Ex_SZCS <- unique(allgenes[is_promoter_ebox & is_Ex_SZCS_all])
overlap_PromoterEbox_Ex_SZTR <- unique(allgenes[is_promoter_ebox & is_Ex_SZTR_all])

overlap_PromoterEbox_In_SCZ
overlap_PromoterEbox_Ex_SZCS
overlap_PromoterEbox_Ex_SZTR

#Pathway Analysis
# Convert overlap genes to Entrez IDs
gene_map_promoter_in <- bitr(
  overlap_PromoterEbox_In_SCZ,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_map_promoter_ex <- bitr(
  overlap_PromoterEbox_Ex_SZCS,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_map_promoter_tr <- bitr(
  overlap_PromoterEbox_Ex_SZTR,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Extract Entrez IDs
entrez_promoter_in <- unique(gene_map_promoter_in$ENTREZID)
entrez_promoter_ex <- unique(gene_map_promoter_ex$ENTREZID)
entrez_promoter_tr <- unique(gene_map_promoter_tr$ENTREZID)

# GO enrichment
go_promoter_in <- enrichGO(
  gene = entrez_promoter_in,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  readable = TRUE
)

go_promoter_ex <- enrichGO(
  gene = entrez_promoter_ex,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  readable = TRUE
)

go_promoter_tr <- enrichGO(
  gene = entrez_promoter_tr,
  universe = entrez_universe,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  readable = TRUE
)

# KEGG enrichment
kegg_promoter_in <- enrichKEGG(
  gene = entrez_promoter_in,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_promoter_ex <- enrichKEGG(
  gene = entrez_promoter_ex,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_promoter_tr <- enrichKEGG(
  gene = entrez_promoter_tr,
  universe = entrez_universe,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Reactome enrichment
react_promoter_in <- enrichPathway(
  gene = entrez_promoter_in,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

react_promoter_ex <- enrichPathway(
  gene = entrez_promoter_ex,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

react_promoter_tr <- enrichPathway(
  gene = entrez_promoter_tr,
  universe = entrez_universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot
dotplot(go_promoter_in, showCategory = 15) +
  ggtitle("GO Biological Process: Promoter EBOX ∩ In_SCZ")

dotplot(go_promoter_ex, showCategory = 15) +
  ggtitle("GO Biological Process: Promoter EBOX ∩ Ex_SZCS")

dotplot(go_promoter_tr, showCategory = 15) +
  ggtitle("GO Biological Process: Promoter EBOX ∩ Ex_SZTR")

dotplot(react_promoter_in, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment: Promoter EBOX ∩ In_SCZ")

dotplot(react_promoter_ex, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment: Promoter EBOX ∩ Ex_SZCS")

dotplot(react_promoter_tr, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment: Promoter EBOX ∩ Ex_SZTR")

dotplot(kegg_promoter_in, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment: Promoter EBOX ∩ In_SCZ")

dotplot(kegg_promoter_ex, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment: Promoter EBOX ∩ Ex_SZCS")

dotplot(kegg_promoter_tr, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment: Promoter EBOX ∩ Ex_SZTR")

#########################################
# SCZ, ASD and BPD Transcriptome(Gandal)#
#########################################
#https://pmc.ncbi.nlm.nih.gov/articles/PMC6443102/#SD8

#Load data
SCZ_AD_BPD_data <- read_excel("Transcriptome_SCZ_AD_BPD.xlsx", sheet = "DGE")
View(SCZ_AD_BPD_data)
colnames(SCZ_AD_BPD_data)

#Slect significant genes for each disorder
SCZ_genes <- SCZ_AD_BPD_data %>%
  filter(SCZ.p.value < 0.05,
         SCZ.fdr < 0.05,) %>%
  select(chr, start, stop, gene_name, gene_type, SCZ.p.value, SCZ.fdr, SCZ.log2FC)
View(SCZ_genes)

ASD_genes <- SCZ_AD_BPD_data %>%
  filter(ASD.p.value < 0.05,
         ASD.fdr < 0.05) %>%
  select(chr, start, stop, gene_name, gene_type, ASD.p.value, ASD.fdr, ASD.log2FC)
View(ASD_genes)

BPD_genes <- SCZ_AD_BPD_data %>%
  filter(BD.p.value < 0.05,
         BD.fdr < 0.05) %>%
  select(chr, start, stop, gene_name, gene_type, BD.p.value, BD.fdr, BD.log2FC)
View(BPD_genes)

## Compare with our TCF4 gene set
# Check overlap of allgenes with each disorder gene list
# I didn't separate the upregulated and downregulated genes in these datasets
is_SCZ_all <- allgenes %in% SCZ_genes$gene_name
sum(is_SCZ_all)

is_ASD_all <- allgenes %in% ASD_genes$gene_name
sum(is_ASD_all)

is_BPD_all <- allgenes %in% BPD_genes$gene_name
sum(is_BPD_all)

# Fisher test for all peak genes
# All peaks vs SCZ
fisher_all_SCZ <- fisher.test(table(is_gene_all, is_SCZ_all), alternative = "greater")
fisher_all_SCZ
table(is_gene_all, is_SCZ_all)

# All peaks vs ASD
fisher_all_ASD <- fisher.test(table(is_gene_all, is_ASD_all), alternative = "greater")
fisher_all_ASD
table(is_gene_all, is_ASD_all)

# All peaks vs BPD
fisher_all_BPD <- fisher.test(table(is_gene_all, is_BPD_all), alternative = "greater")
fisher_all_BPD
table(is_gene_all, is_BPD_all)

# Extract overlapped genes
overlap_all_SCZ <- allgenes[is_gene_all & is_SCZ_all]
overlap_all_ASD <- allgenes[is_gene_all & is_ASD_all]
overlap_all_BPD <- allgenes[is_gene_all & is_BPD_all]

overlap_all_SCZ
overlap_all_ASD
overlap_all_BPD

###########################################
# SCZ GWAS Prioritized Genes (Trubetskoy) #
###########################################
# Load data
Trsky_data <- read_excel("Supplementary Table 12.xlsx", sheet = "ST12 all criteria")
View(Trsky_data)
colnames(Trsky_data)

#Clean gene list
Trsky_genes <- Trsky_data %>%
  select(Symbol.ID) %>%
  rename(gene = Symbol.ID) %>%
  mutate(gene = str_trim(gene)) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct()

View(Trsky_genes)

#Compare with our TCF4 gene set
#Check overlap of allgenes with prioritized SCZ gene list
is_Trsky_all <- allgenes %in% Trsky_genes$gene
sum(is_Trsky_all)

#Fisher test
#ALL peak genes
fisher_all_Trsky <- fisher.test(table(is_gene_all, is_Trsky_all), alternative = "greater")
fisher_all_Trsky
table(is_gene_all, is_Trsky_all)

# EBOX genes
fisher_motif_Trsky <- fisher.test(table(is_gene_motif, is_Trsky_all), alternative = "greater")
fisher_motif_Trsky
table(is_gene_motif, is_Trsky_all)

# Promoter EBOX genes
fisher_promoter_Trsky <- fisher.test(table(is_promoter_ebox, is_Trsky_all), alternative = "greater")
fisher_promoter_Trsky
table(is_promoter_ebox, is_Trsky_all)

# Extract overlapped genes
overlap_all_Trsky <- allgenes[is_gene_all & is_Trsky_all]
overlap_motif_Trsky <- allgenes[is_gene_motif & is_Trsky_all]
overlap_promoterEbox_Trsky <- allgenes[is_promoter_ebox & is_Trsky_all]

overlap_all_Trsky
overlap_motif_Trsky
overlap_promoterEbox_Trsky
