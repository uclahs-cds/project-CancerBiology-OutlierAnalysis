### METABRIC Dataset EDA ###########################################################################

# The following script was utilized to create a Kaplan-Meier plot,
# with survival curves for each of the breast cancer subtypes in the METABRIC dataset.

### Preamble #######################################################################################

# Read in the data 
load("/Users/amaanjsattar/Desktop/2023-07-07_Metabric_Outlier.rda");
source('/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/TCGA.FUNCTIONS.R');
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

# Change subtype labels to official names
meta.clinic.subtypes$Pam50...Claudin.low.subtype <- gsub(
    'claudin-low',
    'Claudin-low', 
    meta.clinic.subtypes$Pam50...Claudin.low.subtype
    );

meta.clinic.subtypes$Pam50...Claudin.low.subtype <- gsub(
    'LumA',
    'Luminal A', 
    meta.clinic.subtypes$Pam50...Claudin.low.subtype
    );

meta.clinic.subtypes$Pam50...Claudin.low.subtype <- gsub(
    'LumB',
    'Luminal B', 
    meta.clinic.subtypes$Pam50...Claudin.low.subtype
    );
# Create a new column to represent survival status as a binary variable.
meta.clinic.subtypes$Survival.Status.Binary <- ifelse(
    meta.clinic.subtypes$Overall.Survival.Status == "0:LIVING", 
    0, 1
    );

####################################################################################################
### KM Plot Setup Part 1: Overall Survival #########################################################
####################################################################################################

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

# Order the rows based on the alphabetical order of the Subtype column
kmplot.cols <- kmplot.cols[order(kmplot.cols$Pam50...Claudin.low.subtype), ]


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
                         levels = c('Basal', 'Claudin-low', 'Her2', 'Luminal A', 
                                    'Luminal B', 'Normal')
                        );

### Plot Creation ##################################################################################
# Create Survival Plot, grouping patients by the subtype factor variable
Subtype.KM.Grouped <- BoutrosLab.plotting.survival::create.km.plot(
    survival.object = Surv(
        kmplot.cols$Overall.Survival..Months.,
        kmplot.cols$Survival.Status.Binary),
    patient.groups = subtype.groups,
    main = 'Overall Survival by PAM50 Subtype: METABRIC Patients',
    main.cex = 2,
    height = 12,
    width = 12,
    xat = c(0, 50, 100, 150, 200, 250, 300, 350, 400),
    xaxis.cex = 1,
    yaxis.cex = 1,
    ylab.cex = 1.5,
    xlab.cex = 1.5,
    xlab.label = 'Overall Survival Time (Months)',
    ylab.label = 'Overall Survival Probability',
    key.groups.cex = 1,
    key.groups.title = 'PAM50 Subtype',
    key.groups.corner = c(-0.2, 0),
    key.groups.title.cex = 1,
    resolution = 300,
    top.padding = 5,
    left.padding = 15,
    key.stats.cex = 1,
    key.stats.corner = c(1, -35),
    statistical.method = 'logrank',
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/META.KM.OS.SUBTYPES.pdf'
);

####################################################################################################
### KM Plot Setup Part 2: Relapse-Free Survival ####################################################
####################################################################################################

# Variables: Disease Free Status, Disease Free Months, Subtypes
    # 
# Variables: Progression Free Status, Progression Free Months, Subtypes
# Variables: Disease Specific Status, Disease Specific Months, Subtypes

# Subtype Data: meta.clinic.subtypes

meta.clinic.subtypes <- binary.survival(
    meta.clinic.subtypes, 
    'Relapse.Free.Status', 
    'Relapse.Free.Binary'
    );

relapse.surv <- create.surv(
    meta.clinic.subtypes, 
    'Relapse.Free.Status..Months.',
    'Relapse.Free.Binary'
    );

# subtype.groups already made
BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/META.KM.RELAPSE.SUBTYPES.tiff',
    survival.object = relapse.surv,
    patient.groups = subtype.groups,
    statistical.method = 'logrank',
    main = 'Relapse-Free Survival by PAM50 Subtype: METABRIC Patients',
    ylab.label = 'Relapse-Free Survival Probability',
    xlab.label = 'Relapse-Free Survival Time (Months)',
    main.cex = 2,
    xaxis.cex = 1,
    xlab.cex = 1.5,
    yaxis.cex = 1,
    ylab.cex = 1.5,
    key.groups.title = 'PAM50 Subtype',
    key.groups.title.cex = 1.1,
    key.groups.cex = 1.4,
    key.groups.corner = c(-0.2, -0.1),
    ylab.axis.padding = 2,
    left.padding = 17,
    top.padding = 1,
    bottom.padding = 10,
    resolution = 300,
    key.stats.cex = 1.5,
    key.stats.corner = c(1, -25)
)


# Suggested TCGA Application (replace last argument, filename):

subtype.km.grouped(tcga.surv,
                   tcga.subtypes,
                   'Subtype-Specific Kaplan-Meier Survival Estimates: TCGA Patients',
                   '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/TCGA.KM.SUBTYPES.pdf'
);

####################################################################################################
### KM Plot Setup Part 3: Disease-Specific Survival ####################################################
####################################################################################################

# Data Cleaning: Identify unique values of Patient Vital Status
# table(meta.clinic.subtypes$Patient.s.Vital.Status)
# Noticing there is 1 empty string value, find the row in which it is contained
# rows_with_empty_strings <- which(meta.clinic.subtypes$Patient.s.Vital.Status == "")
# meta.clinic.subtypes[rows_with_empty_strings, ]
# We see that its Survival.Status.Binary = 1, which means this patient died. 
# As we are not sure how this patient died, we will censor this value.
# Populate the column with a 0 when we convert vital status to binary

# Change the empty string to 'living'
meta.clinic.subtypes$Patient.s.Vital.Status <- ifelse(
    meta.clinic.subtypes$Patient.s.Vital.Status == "",
    "Living", 
    meta.clinic.subtypes$Patient.s.Vital.Status
    );

# Convert to binary and make new column
meta.clinic.subtypes$DSS.Binary <- ifelse(
    meta.clinic.subtypes$Patient.s.Vital.Status == "Died of Disease",
    1, 0
    );

# Create DSS Survival Object
dss.surv <- create.surv(
    meta.clinic.subtypes,
    'Overall.Survival..Months.',
    'DSS.Binary'
    );

# Create KM Plot
BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/META.KM.DSS.SUBTYPES.pdf',
    survival.object = dss.surv,
    patient.groups = subtype.groups,
    statistical.method = 'logrank',
    main = 'Disease-Specific Survival by PAM50 Subtype: METABRIC Patients',
    ylab.label = 'Disease-Specific Survival Probability',
    xlab.label = 'Disease-Specific Survival Time (Months)',
    main.cex = 2,
    xaxis.cex = 1.2,
    xlab.cex = 1.5,
    yaxis.cex = 1.2,
    ylab.cex = 1.5,
    key.groups.title = 'PAM50 Subtype',
    key.groups.title.cex = 1.1,
    key.groups.cex = 1.4,
    key.groups.corner = c(-0.2, -0.1),
    ylab.axis.padding = 2,
    left.padding = 17,
    top.padding = 1,
    bottom.padding = 10,
    resolution = 300,
    key.stats.cex = 1.5,
    key.stats.corner = c(1, -25)
)

