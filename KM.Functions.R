### Functions for Survival Analysis on TCGA Data ###################################################

### Preamble #######################################################################################
# Install dependencies (devtools, deldir).

# library(devtools);

# install.packages('deldir');
# Install Survival Analysis Packages
# install_github('uclahs-cds/public-R-BoutrosLab-statistics-survival');
# install_github('uclahs-cds/public-R-BoutrosLab-prognosticsignature-general');
# install_github('uclahs-cds/public-R-BoutrosLab-plotting-survival');

# Load the BL General Packages.
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.plotting.general);

# Load the BL Survival Analysis Packages.
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.statistics.survival);
library(BoutrosLab.plotting.survival);

# Load the RDA file
load("/Users/amaanjsattar/Desktop/2023-07-07_TCGA_BRCA_Outlier.rda");

### Data Cleaning ##################################################################################

### FUNCTION: subtype.cleaning ######################################
# Description:
# Function description
# Input variables:
# Input variable names, descriptions
# Output variables:
# Output variable names, descriptions

subtype.cleaning <- function(brca.clinic.data) {
    subtype.labels <- c('BRCA_Basal', 'BRCA_Her2', 'BRCA_LumA', 'BRCA_LumB', 'BRCA_Normal')
    brca.clinic.subtypes <- brca.clinic.data[brca.clinic.data$Subtype %in% subtype.labels, ]
    return(brca.clinic.subtypes)
    };
subtype.cleaning(brca.clinic)

subtype.cleaning <- function(brca.clinic.data) {
    subtype.labels <- c('BRCA_Basal', 'BRCA_Her2', 'BRCA_LumA', 'BRCA_LumB', 'BRCA_Normal')
    subtypes.only <- brca.clinic.data[brca.clinic.data$Subtype %in% subtype.labels, ]
    return(subtypes.only)
    };

brca.clinic.subtypes <- subtype.cleaning(brca.clinic)

