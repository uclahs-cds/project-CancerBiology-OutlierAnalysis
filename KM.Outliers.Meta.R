# Load the BL General Packages.
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.plotting.general);

# Load the BL Survival Analysis Packages.
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.statistics.survival);
library(BoutrosLab.plotting.survival);

# Load necessary functions
source('~/Desktop/BIGSUMMER.PROJ/KM_FUNCTIONS.R')

# Load the RDA file
load("/Users/amaanjsattar/Desktop/2023-07-07_Metabric_Outlier.rda");

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


# 1991 patients, 4782 outlier genes
# sum up all of the columns to get per-patient outlier counts
outlier.counts <- colSums(outlier.patient.tag.01);

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
merged_data <- merge(
    meta.clinic.patient.rm,
    outlier.counts.df,
    by = 0
    );


# Set the patient IDs as row names of the merged dataframe
rownames(merged_data) <- merged_data$Row.names;

# Remove the redundant column containing patient IDs
merged_data$Row.names <- NULL;

merged.data.surv <- merged_data[, c('Overall.Survival..Months.',
                                    'Overall.Survival.Status',
                                    'Pam50...Claudin.low.subtype',
                                    'Outliers')
                                ];
colnames(merged.data.surv)[3] <- 'Subtype'



merged.data.surv <- subtype.cleaning(merged.data.surv, c('Basal',
                                                         'claudin-low',
                                                        'Her2',
                                                        'LumA',
                                                        'LumB',
                                                        'Normal')
                                     );
                                  
