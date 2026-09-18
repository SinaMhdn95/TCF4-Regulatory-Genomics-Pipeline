
## Transcription factor 4 (TCF4) data integration
## This script integrates and tests overlap of three datasets:
#  1. chromatin immunoprecipitation sequencing analysis of TCF4 bindings
#  2. microarray analysis of differential gene expression following siRNA knockdown of TCF4
#  3. RNA-seq analysis of genes differentially expressed in postmortem brain of schizophrenia patients

## There are two goals. First, we want to see if genes with TCF4 binding from ChIP-seq are also
## differentially regulated when TCF4 is knocked down. We expect to see a significant overlap
## so this is really just a check of data validity. Then, we want to see if genes regulated by
## TCF4 (a schizophrenia risk factor) are differentially expressed in patients.


## Configure project-relative data paths. Override TCF4_DATA_DIR when data live elsewhere.
source("scripts/00_setup.R")
analysis_dir <- file.path(data_root, "expression", "siRNA_data_integration")

## Read in all genes in genome from RefGene file
genfil = read.csv(file.path(analysis_dir, "refGene_autoX_geneposns_redux.csv"), stringsAsFactors = FALSE);
head(genfil)
allgenes <- genfil$gene
head(allgenes)

## Read in ChIP-seq data for TCF4 and extract genes
tcf4res <- read.table(file.path(data_root, "peaks", "Overlap_4749_exact_newgenes.dat"), header = F)
head(tcf4res)
colnames(tcf4res) <- c("chr.gene", "start.gene", "end.gene", "gene", "longstart.gene", "longend.gene", "nalign", "longstart.K12",
                       "longend.K12", "chr.K12", "start.K12", "end.K12", "score.K12", "Qval.K12", "overlap", "chr.N16",
                       "start.N16", "end.N16", "best.score.N16", "Qval.N16", "mean.score")
head(tcf4res)
ord.tcf4res <- tcf4res[order(tcf4res$mean.score, decreasing = T),]
head(ord.tcf4res)
tgenall <- unique(as.character(ord.tcf4res$gene))

## Create vector for each gene in "allgenes" where 1 = gene with TCF4 binding site and 0 = not
istcf4 = (allgenes %in% tgenall);
sum(istcf4)

## Load genes differentially expressed in Blake TCF4 siRNA paper in SH-SY5Y
blakegenes = readLines(file.path(analysis_dir, "TCF4 KD - 1205 DEG.txt"))
isblake = (allgenes %in% blakegenes);
sum(isblake)

blakeup = readLines(file.path(analysis_dir, "TCF4 blake up n470.txt"))
isblakeup = (allgenes %in% blakeup);
sum(isblakeup)

blakedown = readLines(file.path(analysis_dir, "TCF4 blake down n665.txt"))
isblakedown = (allgenes %in% blakedown);
sum(isblakedown)


## Load genes differentially expressed in Fromer et al gene expression study of SCZ
fromgenes = readLines(file.path(analysis_dir, "Fromer_genes.txt"))
isfrom = (allgenes %in% fromgenes);
sum(isfrom)

fromgenesup = readLines(file.path(analysis_dir, "Fromer-up.txt"))
isfromup = (allgenes %in% fromgenesup);
sum(isfromup)

fromgenesdown = readLines(file.path(analysis_dir, "Fromer-down.txt"))
isfromdown = (allgenes %in% fromgenesdown);
sum(isfromdown)


# Fisher test of crosstab
fisher.test(table(istcf4, isblake), alternative = "g")
fisher.test(table(istcf4, isblakeup), alternative = "g")
fisher.test(table(istcf4, isblakedown), alternative = "g")
fisher.test(table(istcf4, isfrom), alternative = "g")
fisher.test(table(istcf4, isfromup), alternative = "g")
fisher.test(table(istcf4, isfromdown), alternative = "g")
