# Figures for "An Outlier Detection Method in Breast Cancer Analysis"

## Overview

This repository contains the figure-generating scripts for our academic paper titled:

"OutSeekR: A Novel Outlier Detection Method for Comprehensive Gene Expression Analysis in Breast Cancer"

This code is supplementary to the main analysis package, `package-OutSeekR` ([GitHub](https://github.com/uclahs-cds/package-OutSeekR)).

## Repository Contents

* `Figure/`: R project workspace with individual plotting scripts and `renv.lock`
* `outlierAnalysisSupport/`: Common R package used by plotting scripts

## Getting Started

### Prerequisites

* R (version >= 4.3)

### Bootstrapping

Bootstrapping a new development environment requires the following steps:

1. Clone this repository: `git clone https://github.com/uclahs-cds/project-CancerBiology-OutlierAnalysis.git`.
1. Navigate to the `Figure/` subdirectory.
1. Set the environment variables `OUTLIER_DATA_DIR` and `OUTLIER_DATA_FILENAME` to reference the original outlier datafile. These can be set in an `.Renviron` file.
2. Within R, call `renv::restore()` to install the snapshotted packages.

### Usage

The scripts in this repository are run in two phases: first performing the common analysis, then plotting individual figures. The common analysis scripts use the helper function `cache.multiple.computed.variables` to save individual computed variables into `./output/variable-cache/*.rda`. Those variables are reloaded by later scripts using `load.multiple.computed.variables`.

The analysis scripts must be run in order:

1. `1.preparation.R`
2. `2.copy.number.analysis.R`
3. `3.dna.methylation.analysis.R`
4. `4.cell.line.analysis.R`
5. `5.crispr.rnai.analysis.R`

The individual plotting scripts (`Figure/Figure*/Figure*.R`) have no further dependencies and may be run in any order. Each `Figure/Figure*/Figure*.R` script produces one or more figures in the `Figure/output/` directory, along with a `Figure/output/Figure*.txt` session information file.

## License

This project is licensed under the GNU General Public License version 2. See the file LICENSE.md for the terms of the GNU GPL license.

Copyright (C) 2024 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
