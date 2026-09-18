# Analysis pipeline

## 1. Acquire and quality-check inputs

Collect TCF4 ChIP-seq peak data, neural-cell epigenomic annotations, gene-expression studies, and disease-genetics resources from their original repositories or authorized project storage. Raw FASTQ/BAM files and downloaded literature are maintained outside GitHub.

For datasets reprocessed from sequencing reads, alignments and replicate-level peaks are evaluated before downstream integration. MACS2/MACS3 is used for peak calling and IDR is used where replicate information permits reproducibility assessment.

## 2. Harmonize TCF4 peak sets

Peak files from Forrest, McClay, and NPC/iPSC-NPC analyses are standardized to consistent BED conventions and genome builds. Coordinate conversion is performed when required, and consensus or summit-centered intervals are generated for comparable analyses.

Primary scripts:

- `scripts/01_peaks/01_build_mcclay_consensus_bed.R`
- `scripts/01_peaks/02_annotate_npc_peaks.R`
- `scripts/01_peaks/03_annotate_ipsc_npc_peaks.R`
- `scripts/01_peaks/04_extract_peak_genes.R`
- `scripts/01_peaks/05_compare_peaks_and_genes.R`

## 3. Annotate regulatory context and motifs

ChIPseeker and human TxDb/org.Hs.eg.db resources assign peaks to genomic features and nearby genes. E-box sequences are scanned or intersected with peak intervals. MEME-ChIP, FIMO, and CentriMo outputs may be used for motif enrichment and positional analysis.

Regulatory context is expanded with CTCF, histone marks (including H3K27ac, H3K4me3, and H3K27me3), EP300, other signaling-related factors, and chromatin-domain information where available.

Primary scripts:

- `scripts/01_peaks/06_integrate_epigenomic_features.R`
- `scripts/01_peaks/07_scan_ebox_reference.R`
- `scripts/01_peaks/08_intersect_ebox_intervals.R`

## 4. Integrate TCF4 perturbation and expression studies

Gene-expression signatures from TCF4 knockdown or relevant neural differentiation experiments are reanalyzed or compared with peak-associated genes. The current collection includes analyses based on Blake, Hill, Doostparast, Wang, WNT5A/SCZ, and reference datasets.

Scripts are organized in `scripts/02_expression/`. Depending on the source study, analyses use GEOquery, limma, lumi, DESeq2, or published differential-expression tables.

## 5. Connect TCF4 regulation with disease evidence

Peak-associated genes are evaluated against psychiatric and neurodevelopmental resources, including schizophrenia GWAS/locus evidence, cross-ancestry developing-brain analyses, Gandal transcriptomic findings, Ruzicka single-cell results, and additional disease-enrichment analyses.

Scripts are organized in `scripts/03_disease_integration/`.

## 6. Test enrichment and summarize results

Analyses use overlap counts, Fisher's exact tests, direction-of-effect comparisons, and functional enrichment with Gene Ontology, Reactome, and related resources. Selected visualizations are generated under `scripts/04_reporting/`, and an exploratory Shiny prototype is retained in `app/`.

## Current reproducibility boundary

This repository captures the analysis logic and selected outputs. Because the study is ongoing and combines multiple sources with different usage terms, it is not yet a fully automated end-to-end workflow. A future milestone is to replace remaining script-specific paths with the central configuration, add data checksums, lock package versions with `renv`, and add synthetic-data integration tests.
