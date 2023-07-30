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


### FUNCTION: subtype.cleaning #####################################################################

# Description: 
# A function that takes in a dataframe and returns rows corresponding to 
    # subtypes in subtype.labels, c('BRCA_Basal', 'BRCA_Her2', 'BRCA_LumA', 
    # 'BRCA_LumB', 'BRCA_Normal')

# Input variables: 
# (1) brca.data <- [dataframe] a breast cancer dataframe with column 'Subtype' 
# (2) subtype.labels <- [vector of strings] a vector containing subtype labels in brca.data
    # *** In this case, 
    # subtype.labels <- c('BRCA_Basal', 'BRCA_Her2', 'BRCA_LumA', 'BRCA_LumB', 'BRCA_Normal');

# Output variables:
# (1) subtypes.only <- [dataframe] a subset of the input dataframe with rows corresponding to 
    # subtypes in subtype.labels

subtype.cleaning <- function(brca.data, subtype.labels) {
    subtypes.only <- brca.data[brca.data$Subtype %in% subtype.labels, ]
    # Remove rows with NA values in "column_name" using na.omit()
    subtypes.only <- subtypes.only[!is.na(subtypes.only[, 'Overall.Survival..Months.']), ]
    
    return(subtypes.only)
    };

# Suggested TCGA Application: 
brca.clinic.subtypes <- subtype.cleaning(brca.clinic, c('BRCA_Basal',
                                                        'BRCA_Her2',
                                                        'BRCA_LumA',
                                                        'BRCA_LumB',
                                                        'BRCA_Normal')
                                         );

### RENAME SUBTYPES 
# Data Cleaning: Rename Subtypes
brca.clinic.subtypes$Subtype <- gsub(
    'BRCA_Basal',
    'Basal', 
    brca.clinic.subtypes$Subtype
);

# Data Cleaning: Rename Subtypes
brca.clinic.subtypes$Subtype <- gsub(
    'BRCA_Her2',
    'Her2', 
    brca.clinic.subtypes$Subtype
);

# Data Cleaning: Rename Subtypes
brca.clinic.subtypes$Subtype <- gsub(
    'BRCA_LumA',
    'Luminal A', 
    brca.clinic.subtypes$Subtype
);

# Data Cleaning: Rename Subtypes
brca.clinic.subtypes$Subtype <- gsub(
    'BRCA_LumB',
    'Luminal B', 
    brca.clinic.subtypes$Subtype
);

# Data Cleaning: Rename Subtypes
brca.clinic.subtypes$Subtype <- gsub(
    'BRCA_Normal',
    'Normal', 
    brca.clinic.subtypes$Subtype
);

### Creating a Binary Variable for Overall Survival Status #########################################

### FUNCTION: binary.survival ######################################################################
# Description: 
# A function that takes in a dataframe and survival status column,
#   and adds a column to the dataframe containing survival status as a binary variable

# Input variables: 
# (1) brca.data <- [dataframe] a breast cancer dataframe with n columns containing a survival column 
    # with str values '0:LIVING', '1:DECEASED'
# (2) survival.str.col <- [string] the current survival status column in brca.data
# (3) binary.survival.col.name <- [string] the desired name of a new, binary survival status column

# Output variables:
# (1) brca.data <- [dataframe] returns the augmented dataframe, which should have n + 1 columns, 


# Function to convert the column with values in the format '0:somestring' and '1:somestring' to a binary column
binary.survival <- function(brca.data, old.col.name, new.col.name) {
    brca.data[[new.col.name]] <- ifelse(
        grepl('^1:', brca.data[[old.col.name]]), 1, 0)
    return(brca.data)
    };


# OLD IMPLEMENTATION:
# binary.survival <- function(brca.data, str.survival.col.name, binary.survival.col.name) {
    # brca.data[binary.survival.col.name] <- ifelse(
        # brca.data[[str.survival.col.name]] == '0:LIVING', 
        # 0, 1
    # )
    # return(brca.data)
    # };

# Suggested TCGA Application: 
brca.clinic.subtypes <- binary.survival(brca.clinic.subtypes, 
                                        'Overall.Survival.Status',
                                        'Overall.Survival.Binary'
                                        );

### Creating a Survival Object and Plotting Function Combination ###################################

### FUNCTION: create.surv ##########################################################################
# Description: 
# A function that takes in a dataframe and specified columns corresponding to 
    # survival time and survival status. Returns a survival object

# Input variables: 
# (1) brca.data <- [dataframe] a breast cancer dataframe containing columns corresponding to 
    # patient survival time and survival status
# (2) survival.time <- [string] the name of a column in brca.data containing patient survival times
# (3) survival.status <- [string] the name of a column in brca.data 
    # containing **binary** patient survival status

# Output variables:
# (1) surv.obj <-[survival object] a survival object created based on the input columns 

create.surv <- function(brca.data, survival.time, survival.status) {
    surv.obj <- Surv(as.numeric(brca.data[[survival.time]]),
                     brca.data[[survival.status]])
    return(surv.obj)
    };

# Suggested TCGA Application: 
tcga.surv <- create.surv(brca.clinic.subtypes,
                         'Overall.Survival..Months.',
                         'Overall.Survival.Binary');


### FUNCTION: create.subtype.groups ################################################################

# Description: 
# A function that takes in a dataframe, a 'group column' containing group values for each item, 
    # and a vector of levels with values that match the values in the 'group column'

# Input variables: 
# (1) brca.data <- [dataframe] a breast cancer dataframe containing a column 
    # identifying rows with groups
# (2) subtype.col <- [string] a column name pertaining to the column with each row's group value
# (3) subtype.levels <- a vector of the unique group names in subtype.col

# Output variables:
# (1) tcga.subtypes <- [factor] a factor with length corresponding to 
    # the number of rows in brca.data and levels corresponding to subtype.levels
    # [to be used as value of the patient.groups variable in the Kaplan-Meier Plot]

create.subtype.groups <- function(brca.data, subtype.col, subtype.levels) {
    tcga.subtypes <- factor(brca.data[[subtype.col]],
                            levels = subtype.levels)
    return(tcga.subtypes)
    };

# Suggested TCGA Application: 
tcga.subtypes <- create.subtype.groups(brca.clinic.subtypes, 
                                       'Subtype', 
                                       unique(brca.clinic.subtypes$Subtype)
                                       );


### FUNCTION: subtype.km.grouped ###################################################################

# Description: Produces a grouped Kaplan-Meier plot from a survival object, specified groups, and 
    # a title. The plot is saved to a file of the user's choice. 

# Input variables:
# (1) surv.obj <-[survival object] the survival object upon which the plot is based
# (2) subtype.groups <- [factor] a factor with levels corresponding to each specified group 
    # [corresponds to 'patient.groups' in the BL create.km.plot function]
# (3) title <- [string] title of the plot [corresponds to 'main' in the BL create.km.plot function]
# (4) file.path <- [string] a file path and name for the plot 
    # [corresponds to 'filename' in the BL create.km.plot function]

# Output variables:
# (1) grouped.plot <- [saved as filename to specified path] the resulting Kaplan-Meier plot, with 
    # parameters specified as inputs to the wrapper function and 
    # as defaults specified in the function

# ***NOTE*** Formatting (e.g. label sizing, resolution, dimensions, etc.) included in function body.
# To modify formatting, add arguments to the function or modify existing attached argument values


subtype.km.grouped <- function(
        surv.obj,
        subtype.groups,
        title,
        file.path) {
    grouped.plot <- BoutrosLab.plotting.survival::create.km.plot(
        height = 12,
        width = 12,
        filename = file.path,
        survival.object = surv.obj,
        patient.groups = subtype.groups,
        statistical.method = 'logrank',
        main = title,
        main.cex = 1.5,
        xaxis.cex = 1,
        xlab.cex = 1.7,
        yaxis.cex = 1,
        ylab.cex = 1.7,
        key.groups.title = 'PAM50 Subtype',
        key.groups.title.cex = 1,
        key.groups.cex = 1,
        resolution = 300,
        key.stats.cex = 1,
        key.stats.corner = c(1, -15)
    )
    
    return(grouped.plot)
        };

# Suggested TCGA Application (replace last argument, filename):

# subtype.km.grouped(tcga.surv,
                   # tcga.subtypes,
                   # 'Subtype-Specific Kaplan-Meier Survival Estimates: TCGA Patients',
                  #  '/Users/amaanjsattar/Desktop/TCGA.KM.SUBTYPES.tiff'
                   # );

### TCGA Overall Survival by Subtype ###############################################################

BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/TCGA.KM.OS.SUBTYPES.pdf',
    survival.object = tcga.surv,
    patient.groups = tcga.subtypes,
    statistical.method = 'logrank',
    main = 'Overall Survival by PAM50 Subtype: TCGA Patients',
    main.cex = 2,
    xaxis.cex = 1,
    xlab.cex = 1.7,
    yaxis.cex = 1,
    ylab.cex = 1.7,
    xlab.label = 'Overall Survival Time (Months)',
    ylab.label = 'Overall Survival Probability',
    key.groups.title = 'PAM50 Subtype',
    key.groups.title.cex = 1,
    key.groups.cex = 1,
    top.padding = 1,
    left.padding = 20,
    resolution = 300,
    key.stats.cex = 1,
    key.stats.corner = c(1, -50),
    key.groups.corner = c(-0.2, 0)
)

### TCGA Progression-Free Survival by Subtype ######################################################



# Create a new column to represent progression-free survival status as a binary variable.
brca.clinic.subtypes <- binary.survival(
    brca.clinic.subtypes,
    'Progression.Free.Status',
    'PFS.Binary'
    );

tcga.pfs.surv <- create.surv(
    brca.clinic.subtypes,
    'Progress.Free.Survival..Months.',
    'PFS.Binary'
    );


BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/TCGA.KM.PFS.SUBTYPES.pdf',
    survival.object = tcga.pfs.surv,
    patient.groups = tcga.subtypes,
    statistical.method = 'logrank',
    main = 'Progression-Free Survival by PAM50 Subtype: TCGA Patients',
    main.cex = 2,
    xaxis.cex = 1,
    xlab.cex = 1.7,
    yaxis.cex = 1,
    ylab.cex = 1.7,
    xlab.label = 'Progression-Free Survival Time (Months)',
    ylab.label = 'Progression-Free Survival Probability',
    key.groups.title = 'PAM50 Subtype',
    key.groups.title.cex = 1,
    key.groups.cex = 1,
    top.padding = 1,
    left.padding = 20,
    resolution = 300,
    key.stats.cex = 1,
    key.stats.corner = c(1, -50),
    key.groups.corner = c(-0.2, 0)
    );

### TCGA Disease-Specific Survival by Subtype ######################################################
brca.clinic.subtypes <- binary.survival(
    brca.clinic.subtypes,
    'Disease.specific.Survival.status',
    'DSS.Binary'
    );

tcga.dss.surv <- create.surv(
    brca.clinic.subtypes,
    'Months.of.disease.specific.survival',
    'DSS.Binary'
    );

BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/TCGA.KM.DSS.SUBTYPES.pdf',
    survival.object = tcga.dss.surv,
    patient.groups = tcga.subtypes,
    statistical.method = 'logrank',
    main = 'Disease-Specific Survival by PAM50 Subtype: TCGA Patients',
    main.cex = 2,
    xaxis.cex = 1,
    xlab.cex = 1.7,
    yaxis.cex = 1,
    ylab.cex = 1.7,
    xlab.label = 'Disease-Specific Survival Time (Months)',
    ylab.label = 'Disease-Specific Survival Probability',
    key.groups.title = 'PAM50 Subtype',
    key.groups.title.cex = 1,
    key.groups.cex = 1,
    top.padding = 1,
    left.padding = 20,
    resolution = 300,
    key.stats.cex = 1,
    key.stats.corner = c(1, -50),
    key.groups.corner = c(-0.2, 0)
);


### TCGA Disease-Free Survival by Subtype ######################################################
brca.clinic.subtypes <- binary.survival(
    brca.clinic.subtypes,
    'Disease.Free.Status',
    'DFS.Binary'
)

tcga.dfs.surv <- create.surv(
    brca.clinic.subtypes,
    'Disease.Free..Months.',
    'DFS.Binary'
)

BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/TCGA.KM.DFS.SUBTYPES.pdf',
    survival.object = tcga.dfs.surv,
    patient.groups = tcga.subtypes,
    statistical.method = 'logrank',
    main = 'Disease-Free Survival by PAM50 Subtype: TCGA Patients',
    main.cex = 2,
    xaxis.cex = 1,
    xlab.cex = 1.7,
    yaxis.cex = 1,
    ylab.cex = 1.7,
    xlab.label = 'Disease-Free Survival Time (Months)',
    ylab.label = 'Disease-Free Survival Probability',
    key.groups.title = 'PAM50 Subtype',
    key.groups.title.cex = 1,
    key.groups.cex = 1,
    top.padding = 1,
    left.padding = 20,
    resolution = 300,
    key.stats.cex = 1,
    key.stats.corner = c(1, -50),
    key.groups.corner = c(-0.2, 0)
);


### KM PLOTS BY OUTLIER GENE COUNT #################################################################


### DATA CLEANING ##################################################################################
tcga.outlier.counts <- as.data.frame(
    colSums(outlier.patient.tag.01.t.p.order)
);

colnames(tcga.outlier.counts) <- 'Outliers';

tcga.outlier.counts$Outliers[tcga.outlier.counts$Outliers >= 3] <- "3+";

tcga.merged <- merge(
    brca.clinic.subtypes,
    tcga.outlier.counts,
    by = 0
);

# Convert row names in 'outlier.counts.df' to use dashes instead of periods
rownames(tcga.outlier.counts) <- gsub(
    "\\.",
    "-",
    rownames(tcga.outlier.counts)
);

tcga.merged <- merge(
    transform(
        tcga.outlier.counts, 
        ModifiedSampleID = substr(rownames(tcga.outlier.counts), 
                                  1, 
                                  nchar(rownames(tcga.outlier.counts)) - 1)),
                     transform(
                         brca.clinic.subtypes, 
                         ModifiedSampleID = substr(Sample.ID, 
                                                   1, 
                                                   nchar(Sample.ID) - 1)),
                     by = "ModifiedSampleID",
                     all.x = TRUE
    );

tcga.outlier.counts$Sample.ID <- rownames(tcga.outlier.counts)

# Function to extract the modified 'Sample.ID' (excluding the last character)
getModifiedSampleID <- function(x) substr(x, 1, nchar(x) - 1)

# Extract modified 'Sample.ID' tcga.outlier.counts (Match with brca.clinic.subtypes)
tcga.outlier.counts$ModifiedSampleID <- getModifiedSampleID(tcga.outlier.counts$Sample.ID)

# Verify matches
matches <- tcga.outlier.counts$ModifiedSampleID %in% brca.clinic.subtypes$Sample.ID

tcga.merged <- merge(tcga.outlier.counts[matches, ], brca.clinic.subtypes,
                     by.x = "ModifiedSampleID", by.y = "Sample.ID", all.x = TRUE)


rownames(tcga.merged) <- tcga.merged$Sample.ID;

tcga.merged$Sample.ID <- NULL;

order.levels <- c('3+', '2', '1', '0')


# Reorder the dataframe based on the "Outliers" column
tcga.merged$Outliers <- factor(tcga.merged$Outliers, levels = order.levels);

tcga.merged <- tcga.merged[order(tcga.merged$Outliers), ];

####################################################################################################
### KM Plot Setup Part 1: Overall Survival #########################################################
####################################################################################################

tcga.merged <- binary.survival(
    tcga.merged,
    'Overall.Survival.Status',
    'OSS.Binary'
);

tcga.outlier.os.surv <- create.surv(
    tcga.merged, 
    'Overall.Survival..Months.', 
    'OSS.Binary'
);

outlier.groups <- factor(tcga.merged$Outliers,
                         levels = order.levels)

BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/Plots/TCGA.KM.OS.OUTLIERS.pdf',
    survival.object = tcga.outlier.os.surv,
    patient.groups = outlier.groups,
    statistical.method = 'logrank',
    main = 'Overall Survival by Outlier Gene Count: TCGA Patients',
    ylab.label = 'Overall Survival Probability',
    xlab.label = 'Overall Survival Time (Months)',
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
### KM Plot Setup Part 2: Disease-Specific Survival #########################################################
####################################################################################################

tcga.merged <- binary.survival(tcga.merged,
                               'Disease.specific.Survival.status',
                                 'DSS.Binary')

tcga.outlier.dss.surv <- create.surv(tcga.merged,
                                     'Months.of.disease.specific.survival',
                                     'DSS.Binary')

BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/Plots/TCGA.KM.DSS.OUTLIERS.pdf',
    survival.object = tcga.outlier.dss.surv,
    patient.groups = outlier.groups,
    statistical.method = 'logrank',
    main = 'Disease-Specific Survival by Outlier Gene Count: TCGA Patients',
    ylab.label = 'Disease-Specific Survival Probability',
    xlab.label = 'Disease-Specific Survival Time (Months)',
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
### KM Plot Setup Part 3: Progression-Free Survival #########################################################
####################################################################################################

tcga.merged <- binary.survival(
    tcga.merged,
    'Progression.Free.Status',
    'PFS.Binary'
);

tcga.outlier.pfs.surv <- create.surv(
    tcga.merged,
    'Progress.Free.Survival..Months.',
    'PFS.Binary'
);

BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/Plots/TCGA.KM.PFS.OUTLIERS.pdf',
    survival.object = tcga.outlier.pfs.surv,
    patient.groups = outlier.groups,
    statistical.method = 'logrank',
    main = 'Progression-Free Survival by Outlier Gene Count: TCGA Patients',
    ylab.label = 'Progression-Free Survival Probability',
    xlab.label = 'Progression-Free Survival Time (Months)',
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
### KM Plot Setup Part 3: Disease-Free Survival #########################################################
####################################################################################################

tcga.merged <- binary.survival(
    tcga.merged,
    'Disease.Free.Status',
    'DFS.Binary'
    );

tcga.outlier.dfs.surv <- create.surv(
    tcga.merged,
    'Disease.Free..Months.',
    'DFS.Binary'
    );


BoutrosLab.plotting.survival::create.km.plot(
    height = 12,
    width = 12,
    filename = '/Users/amaanjsattar/Desktop/BIGSUMMER.PROJ/Plots/TCGA.KM.DFS.OUTLIERS.pdf',
    survival.object = tcga.outlier.dfs.surv,
    patient.groups = outlier.groups,
    statistical.method = 'logrank',
    main = 'Disease-Free Survival by Outlier Gene Count: TCGA Patients',
    ylab.label = 'Disease-Free Survival Probability',
    xlab.label = 'Disease-Free Survival Time (Months)',
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
