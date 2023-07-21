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
    return(subtypes.only)
    };

# Suggested TCGA Application: 
brca.clinic.subtypes <- subtype.cleaning(brca.clinic, c('BRCA_Basal',
                                                        'BRCA_Her2',
                                                        'BRCA_LumA',
                                                        'BRCA_LumB',
                                                        'BRCA_Normal')
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

binary.survival <- function(brca.data, str.survival.col.name, binary.survival.col.name) {
    brca.data[binary.survival.col.name] <- ifelse(
        brca.data[[str.survival.col.name]] == '0:LIVING', 
        0, 1
    )
    return(brca.data)
    };

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
                                       c('BRCA_Basal',
                                         'BRCA_Her2',
                                         'BRCA_LumA',
                                         'BRCA_LumB',
                                         'BRCA_Normal')
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

subtype.km.grouped(tcga.surv,
                   tcga.subtypes,
                   'TCGA Breast Cancer Subtype-Specific Patient Survival',
                   '/Users/amaanjsattar/Desktop/tcga.function.plot.tiff'
                   );
