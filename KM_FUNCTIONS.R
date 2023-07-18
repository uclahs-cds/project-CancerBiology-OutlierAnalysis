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

subtype.labels <- c('BRCA_Basal', 'BRCA_Her2', 'BRCA_LumA', 'BRCA_LumB', 'BRCA_Normal');

subtype.cleaning <- function(brca.clinic.data) {
    subtypes.only <- brca.clinic.data[brca.clinic.data$Subtype %in% subtype.labels, ]
    return(subtypes.only)
    };

brca.clinic.subtypes <- subtype.cleaning(brca.clinic);

### Creating a Binary Variable for Overall Survival Status
binary.survival <- function(brca.clinic.data) {
    brca.clinic.data$Survival.Binary <- ifelse(
        brca.clinic.data$Overall.Survival.Status == "0:LIVING", 
        0, 1
    )
    return(brca.clinic.data)
    };

brca.clinic.subtypes <- binary.survival(brca.clinic.subtypes);

### Creating a Survival Object and Plotting Function Combination
create.surv <- function(brca.clinic.data) {
    surv.obj <- Surv(brca.clinic.data$Overall.Survival..Months.,
                     brca.clinic.data$Survival.Binary)
    return(surv.obj)
    };

tcga.surv <- create.surv(brca.clinic.subtypes);

tcga.subtypes <- factor(brca.clinic.subtypes$Subtype,
                        levels = c('BRCA_Basal', 'BRCA_Her2', 
                                   'BRCA_LumA', 'BRCA_LumB', 'BRCA_Normal')
                        );

### Create plot
Subtype.KM.Grouped <- BoutrosLab.plotting.survival::create.km.plot(
    survival.object = tcga.surv,
    patient.groups = tcga.subtypes,
    main = 'TCGA Breast Cancer Subtype-Specific Patient Survival',
    main.cex = 1.5,
    height = 12,
    width = 12,
    # xat = c(0, 50, 100, 150, 200, 250, 300, 350, 400),
    xaxis.cex = 1,
    yaxis.cex = 1,
    ylab.cex = 1.7,
    xlab.cex = 1.7,
    key.groups.cex = 1,
    key.groups.title = 'Subtype',
    key.groups.title.cex = 1,
    resolution = 300,
    key.stats.cex = 1,
    key.stats.corner = c(1, -15), 
    statistical.method = 'logrank',
    filename = '/Users/amaanjsattar/Desktop/tcga.km.grouped.tiff'
);

