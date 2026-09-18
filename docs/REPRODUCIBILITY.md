# Reproducibility guide

## Software layers

- **R (4.3 or later recommended):** statistical analysis, Bioconductor genomics, enrichment, visualization, and the Shiny prototype.
- **Python (3.10 or later recommended):** available for workflow helpers and command-line tool interoperability.
- **Command-line bioinformatics:** MACS2/MACS3, IDR, BEDTools, SAMtools, and MEME Suite.

Create the command-line environment with:

```bash
conda env create -f environment.yml
conda activate tcf4-regulatory-genomics
```

R dependencies are declared in `DESCRIPTION`. Bioconductor packages should be installed with `BiocManager`; CRAN packages can be installed with `install.packages`. Because the code was developed iteratively, install packages as required by the selected script and record the resulting `sessionInfo()` with each formal result.

## Configuration

Copy `config/example_config.yml` to `config/config.yml`, then update local paths. The real configuration is ignored by Git so local paths and storage details are not published. Scripts being modernized should source `scripts/00_setup.R` and derive paths from `project_root` or `TCF4_DATA_DIR`.

## Recommended execution order

1. Prepare authorized source datasets outside version control.
2. Generate or validate peak sets with MACS2/MACS3 and IDR.
3. Run `scripts/01_peaks/` in numeric order where the required inputs are available.
4. Run the relevant dataset-specific scripts in `scripts/02_expression/`.
5. Run disease-integration analyses in `scripts/03_disease_integration/`.
6. Generate summary figures with `scripts/04_reporting/`.

Not every numbered script depends on every preceding script; several are parallel, dataset-specific analyses.

## Quality-control checklist

- Record genome build for every genomic interval file.
- Confirm chromosome naming conventions before overlap operations.
- Keep raw, interim, and processed data distinct.
- Record source accession, download date, checksum, and license/data-use terms.
- Set random seeds for stochastic procedures.
- Save `sessionInfo()` and command-line tool versions with formal outputs.
- Review all preliminary results before publication or external sharing.
