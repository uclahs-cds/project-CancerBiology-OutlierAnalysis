## Getting Started

Most of the shared analysis code has been abstracted into an `outlierAnalysisSupport` R package, which includes all of the appropriate dependencies. This project uses [renv](https://rstudio.github.io/renv/index.html) for package isolation and reproducibility.

Bootstrapping a new development environment requires the following steps:

1. Set the environment variables `OUTLIER_DATA_DIR` and `OUTLIER_DATA_FILENAME` to reference the original outlier datafile. These can be set up in an `.Renviron` file as well.
2. Within R, call `renv::restore()` to install all of the snapshotted packages.
3. Within R, call `devtools::install('../outlierAnalysisSupport')` to install the helper package.

### Analysis

The scripts in this repository are run in two phases: performing the common analysis, then plotting individual figures. The common analysis scripts use the helper function `cache.multiple.computed.variables` to save individual computed variables into `./output/variable-cache/*.rda`. Those variables are reloaded by later scripts using `load.multiple.computed.variables`.

The analysis scripts must be run in order:
1. `1.preparation.R`
2. `2.copy.number.analysis.R`
3. `3.dna.methylation.analysis.R`
4. `4.cell.line.analysis.R`
5. `5.crispr.rnai.analysis.R`

The individual plotting scripts (`Figure*/Figure*.R`) have no further dependencies and may be run in any order. Each `Figure*/Figure*.R` script produces one or more figures in the `./output/` directory, along with a `./output/Figure*.txt` session information file.


