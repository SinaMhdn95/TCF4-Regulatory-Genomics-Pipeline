# Data sources and provenance

The analyses reference multiple public or collaborator-maintained datasets. Always verify the original publication, accession record, genome build, sample metadata, and current data-use terms before reproducing or sharing an analysis.

| Resource represented in the code | Identifier or role | Typical use in this project |
|---|---|---|
| Forrest TCF4 ChIP-seq | GEO: GSE96915 | TCF4 binding and peak-set comparison |
| McClay TCF4 study | GEO: GSE112704 | Replicate/consensus TCF4 peaks |
| Blake TCF4 knockdown | GEO: GSE48367 | Perturbation-expression signature |
| Hill TCF4 knockdown | GEO: GSE62085 | Neural progenitor perturbation analysis |
| Doostparast neural differentiation | GEO: GSE128333 | Differentiation-stage expression and pathways |
| Gandal brain transcriptomics | Study supplementary tables | DTE/DTU and psychiatric-disorder integration |
| Ruzicka single-cell resource | Study supplementary data | Cell-type and disease-expression enrichment |
| Psychiatric GWAS resources | Ripke/Trubetskoy and cross-ancestry tables | Locus/gene overlap and enrichment |
| Cistrome resources | Search/export results | Regulatory-profile similarity |
| Human genome annotations | hg19 and hg38 | Coordinates, genes, and genomic features |

Additional resources referenced by individual scripts include ENCODE-derived regulatory features, SCHEMA evidence, pathway databases, and study-specific supplementary tables. These inputs are excluded from the repository because of size, licensing, or provenance requirements.

Before any public release, replace this working inventory with a versioned manifest containing exact download URLs, access dates, checksums, genome builds, and primary citations for every input.
