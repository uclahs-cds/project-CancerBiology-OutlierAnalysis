# Load the BL General Packages.
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.plotting.general);

# Load the BL Survival Analysis Packages.
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.statistics.survival);
library(BoutrosLab.plotting.survival);

# Load necessary functions
# source('~/Desktop/BIGSUMMER.PROJ/KM_FUNCTIONS.R')

# Load the RDA file
# load("/Users/amaanjsattar/Desktop/2023-07-07_Metabric_Outlier.rda");

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

# sum up all of the columns to get per-patient outlier counts, and 
# # create a new dataframe with indices as patient ID and values as outlier counts
outlier.counts <- as.data.frame(
    colSums(outlier.patient.tag.01)
    );

# create a new dataframe with 
    # indices corresponding to patient ID and values as per-patient outlier counts
outlier.counts.df <- as.data.frame(outlier.counts);

# name the column with outlier counts 'Outliers'
colnames(outlier.counts.df) <- 'Outliers';

# replace values >= 3 with '3+'
outlier.counts.df$Outliers[outlier.counts.df$Outliers >= 3] <- "3+";

# Convert row names in 'outlier.counts.df' to use dashes instead of periods
rownames(outlier.counts.df) <- gsub(
    "\\.",
    "-",
    rownames(outlier.counts.df)
    );

# Get the common patient IDs between the two dataframes
merged.data <- merge(
    meta.clinic.patient.rm,
    outlier.counts.df,
    by = 0
    );


# Set the patient IDs as row names of the merged dataframe
rownames(merged.data) <- merged.data$Row.names;

# Remove the redundant column containing patient IDs
merged.data$Row.names <- NULL;

# Create a column with binary Survival Status in merged_data
merged.data <- binary.survival(
    merged.data, 
    'Overall.Survival.Status',
    'Overall.Survival.Binary'
    );

# Subset the dataframe to only contain necessary columns
merged.data <- merged.data[, c('Overall.Survival..Months.',
                                    'Overall.Survival.Binary',
                                    'Outliers')
                                ];

# Create Survival Object
merged.data.surv <- create.surv(
    merged.data, 
    'Overall.Survival..Months.', 
    'Overall.Survival.Binary'
    );

# Create a factor for the outlier groupings, making sure the levels correspond correctly
outlier.groups <- factor(merged.data$Outliers,
                         levels = unique(merged.data$Outliers));

km.outliers.meta <- subtype.km.grouped(merged.data.surv,
                                       outlier.groups,
                                       'Kaplan-Meier Plot Grouped by Subtype',
                                       '/Users/amaanjsattar/Desktop/testplot.tiff')
