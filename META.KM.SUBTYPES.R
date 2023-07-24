### METABRIC Dataset EDA ###########################################################################

# The following script was utilized to create a Kaplan-Meier plot,
# with survival curves for each of the breast cancer subtypes in the METABRIC dataset.

### Preamble #######################################################################################

# Read in the data 
# load("/Users/amaanjsattar/Desktop/2023-07-07_Metabric_Outlier.rda");

# Install dependencies (devtools, deldir).
library(devtools);


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

### Data Cleaning ##################################################################################

# Remove any NA Values.
meta.clinic.patient.rm <- na.omit(meta.clinic.patient);

# Remove any empty string values in the subtype column.
meta.clinic.patient.rm <- meta.clinic.patient.rm[
    meta.clinic.patient.rm$Pam50...Claudin.low.subtype != '',
    ];

# Create a separate dataframe with the metadata.
meta.clinic.patient.metadata <- meta.clinic.patient[1:4,];

# Remove the metadata from meta.clinic.patient.rm.
meta.clinic.patient.rm <- meta.clinic.patient.rm[-c(1:4), ];

### Subtype Identification and Grouping ############################################################

# Create a frequency table to identify the different subtypes in the dataset.
unique.subtypes <- table(meta.clinic.patient.rm$Pam50...Claudin.low.subtype);

# Create a vector of subtype labels.
subtype.labels <- c('LumA', 'LumB', 'claudin-low', 'Her2', 'Basal', 'Normal');

# Filter the data to include the rows where subtype is in subtype.labels.
meta.clinic.subtypes <- meta.clinic.patient.rm[
    meta.clinic.patient.rm$Pam50...Claudin.low.subtype %in% subtype.labels, ];

# Create a new column to represent survival status as a binary variable.
meta.clinic.subtypes$Survival.Status.Binary <- ifelse(
    meta.clinic.subtypes$Overall.Survival.Status == "0:LIVING", 
    0, 1
    );

### Setup for Kaplan-Meier Plot Inputs #############################################################

# Subset the dataframe, extracting columns corresponding to:
# 1) SUBTYPE LABEL,
# 2) SURVIVAL TIME,
# 3) SURVIVAL STATUS (BINARY)

kmplot.cols <- subset(
    meta.clinic.subtypes,
    select = c('Pam50...Claudin.low.subtype', 
               'Overall.Survival..Months.', 
               'Survival.Status.Binary')
    );

# Make Overall.Survival..Months. numeric-type
kmplot.cols$Overall.Survival..Months. <- as.numeric(kmplot.cols$Overall.Survival..Months.);

# Create Survival Object (to be used as input for kmplot)
# Survival Object Inputs include:
# 1) SURVIVAL TIME
# 2) SURVIVAL STATUS
meta.surv.obj <- Surv(
    kmplot.cols$Overall.Survival..Months.,
    kmplot.cols$Survival.Status.Binary
                      );

# Make subtype a factor variable with levels corresponding to each subtype label
subtype.groups <- factor(kmplot.cols$Pam50...Claudin.low.subtype,
                         levels = unique(kmplot.cols$Pam50...Claudin.low.subtype)
                        );

### Plot Creation ##################################################################################
# Create Survival Plot, grouping patients by the subtype factor variable
Subtype.KM.Grouped <- BoutrosLab.plotting.survival::create.km.plot(
    survival.object = Surv(
        kmplot.cols$Overall.Survival..Months.,
        kmplot.cols$Survival.Status.Binary),
    patient.groups = subtype.groups,
    main = 'METABRIC Breast Cancer Subtype-Specific Patient Survival',
    main.cex = 1.5,
    height = 12,
    width = 12,
    xat = c(0, 50, 100, 150, 200, 250, 300, 350, 400),
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
    filename = '/Users/amaanjsattar/Desktop/simple.survival.plot.six.groups.tiff'
);


