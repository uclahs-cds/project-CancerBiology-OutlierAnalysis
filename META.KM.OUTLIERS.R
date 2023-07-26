# Load the BL General Packages.
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.plotting.general);

# Load the BL Survival Analysis Packages.
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.statistics.survival);
library(BoutrosLab.plotting.survival);

# Load necessary functions
source('/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/TCGA.FUNCTIONS.R')

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



order_levels <- c("3+", "2", "1", "0")

# Reorder the dataframe based on the "Outliers" column
merged.data$Outliers <- factor(merged.data$Outliers, levels = order_levels)
merged.data <- merged.data[order(merged.data$Outliers), ]

####################################################################################################
### KM Plot Setup Part 1: Overall Survival #########################################################
####################################################################################################
# Create a column with binary OverallSurvival Status in merged_data
merged.data <- binary.survival(
    merged.data, 
    'Overall.Survival.Status',
    'Overall.Survival.Binary'
    );


                               

# Create Survival Object
merged.data.surv <- create.surv(
    merged.data, 
    'Overall.Survival..Months.', 
    'Overall.Survival.Binary'
    );

# Create a factor for the outlier groupings, making sure the levels correspond correctly
outlier.groups <- factor(merged.data$Outliers,
                         levels = unique(merged.data$Outliers)
                         );

km.outliers.meta <- subtype.km.grouped(merged.data.surv,
                                       outlier.groups,
                                       'Kaplan-Meier Plot Grouped by Subtype',
                                       '/Users/amaanjsattar/Desktop/META.KM.OUTLIERS.tiff'
                                       );


BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/Plots/META.KM.OS.OUTLIERS.pdf',
    survival.object = merged.data.surv,
    patient.groups = outlier.groups,
    statistical.method = 'logrank',
    main = 'Overall Survival by Outlier Gene Count: METABRIC Patients',
    ylab.label = 'Overall Survival Probability',
    xlab.label = 'Overall Survival Time: Months',
    main.cex = 2,
    xaxis.cex = 1.2,
    xlab.cex = 1.7,
    yaxis.cex = 1.2,
    ylab.cex = 1.7,
    key.groups.title = '# Outlier Genes',
    key.groups.title.cex = 1,
    key.groups.cex = 1.4,
    key.groups.corner = c(-0.4, -0.2),
    ylab.axis.padding = 0,
    left.padding = 10,
    top.padding = 5,
    risk.label.pos = -40,
    resolution = 400,
    key.stats.cex = 1.2,
    key.stats.corner = c(1, -45)
);

####################################################################################################
### KM OutlierPlot Setup Part 2: Relapse-Free Survival #############################################
####################################################################################################
merged.data <- binary.survival(
    merged.data, 
    'Relapse.Free.Status', 
    'Relapse.Free.Binary'
    );

outlier.relapse.surv <- create.surv(
    merged.data, 
    'Relapse.Free.Status..Months.',
    'Relapse.Free.Binary'
    );

BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/META.KM.RFS.OUTLIERS.pdf',
    survival.object = outlier.relapse.surv,
    patient.groups = outlier.groups,
    statistical.method = 'logrank',
    main = 'Relapse-Free Survival by Outlier Gene Count: METABRIC Patients',
    ylab.label = 'Relapse-Free Survival Probability',
    xlab.label = 'Relapse-Free Survival Time (Months)',
    main.cex = 2,
    xaxis.cex = 1.2,
    xlab.cex = 1.7,
    yaxis.cex = 1.2,
    ylab.cex = 1.7,
    key.groups.title = '# Outlier Genes',
    key.groups.title.cex = 1,
    key.groups.cex = 1.4,
    key.groups.corner = c(-0.4, -0.2),
    ylab.axis.padding = 0,
    left.padding = 10,
    top.padding = 5,
    risk.label.pos = -40,
    resolution = 400,
    key.stats.cex = 1.2,
    key.stats.corner = c(1, -45)
    );

####################################################################################################
### KM OutlierPlot Setup Part 2: Disease-Specific Survival #########################################
####################################################################################################

merged.data$DSS.Binary <- ifelse(
    merged.data$Patient.s.Vital.Status == "Died of Disease",
    1, 0
    );

dss.outlier.surv <- create.surv(merged.data, 
                                'Overall.Survival..Months.',
                                'DSS.Binary'
                                );
# Create KM Plot
BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/Plots/META.KM.DSS.OUTLIERS.pdf',
    survival.object = dss.outlier.surv,
    patient.groups = outlier.groups,
    statistical.method = 'logrank',
    main = 'Disease-Specific Survival by Outlier Gene Count: METABRIC Patients',
    ylab.label = 'Disease-Specific Survival Probability',
    xlab.label = 'Disease-Specific Survival Time (Months)',
    main.cex = 2,
    xaxis.cex = 1.2,
    xlab.cex = 1.5,
    yaxis.cex = 1.2,
    ylab.cex = 1.5,
    key.groups.title = '# Outlier Genes',
    key.groups.title.cex = 1.1,
    key.groups.cex = 1.4,
    key.groups.corner = c(-0.5, -0.1),
    ylab.axis.padding = 2,
    left.padding = 17,
    top.padding = 1,
    bottom.padding = 10,
    resolution = 300,
    key.stats.cex = 1.5,
    key.stats.corner = c(1, -38)
    );
