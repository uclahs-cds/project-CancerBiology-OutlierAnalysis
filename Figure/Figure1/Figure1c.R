### HISTORY #####################################################################
# This script generates a scatter plot with a smoothed spline curve, which visualizes
# the number of patients required to observe an outlier gene. The x-axis represents
# the required sample size on a log scale, and the y-axis represents the percentage
# of patients. An extrapolation line is added to the plot.
# Date: 2024-08-12

### DESCRIPTION #################################################################
# The script processes patient data to calculate the fraction of patients per required
# sample size, applies a logarithmic transformation, and then fits a smoothing spline
# to the data. The resulting plot shows the relationship between the required sample
# size and the percentage of patients needed to observe an outlier gene.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-09-10_Figure1.rda'));

### DATA PREPARATION ############################################################

five.data.outlier.symbol.na <- na.omit(five.data.outlier.symbol);
# Remove NA values from 'five.data.outlier.symbol'

### 1. MATADOR
# Find the position of "_" in Metador data row names
outlier.patient.tag.01.metador.pos <- which(
    strsplit(rownames(outlier.patient.tag.01.metador), '')[[1]] == '_'
    );
# Extract the substring after "_"
outlier.patient.tag.01.metador.symbol <- substring(
    rownames(outlier.patient.tag.01.metador),
    outlier.patient.tag.01.metador.pos + 1
    );
# Match and extract rows from Metador data based on 'five.data.outlier.symbol.na'
outlier.patient.tag.01.metador.match.five <- outlier.patient.tag.01.metador[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.metador.symbol),
    ];

### 2. TCGA-BRCA
outlier.patient.tag.01.brca.symbol <- fpkm.tumor.symbol.filter.brca[
    rownames(outlier.patient.tag.01.brca),
    ]$Symbol;
outlier.patient.tag.01.brca.match.five <- outlier.patient.tag.01.brca[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.brca.symbol),
    ];


### 3. METABRIC
outlier.patient.tag.01.meta.symbol <- fpkm.tumor.symbol.filter.meta.symbol[
    rownames(outlier.patient.tag.01.meta),
    ]$Symbol;
outlier.patient.tag.01.meta.match.five <- outlier.patient.tag.01.meta[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.meta.symbol),
    ];


### 4. ISPY-2
outlier.patient.tag.01.ispy.symbol <- rownames(outlier.patient.tag.01.ispy);
outlier.patient.tag.01.ispy.match.five <- outlier.patient.tag.01.ispy[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.ispy.symbol),
    ];


### 5. ICGC BRCA-EU
outlier.patient.tag.01.icgc.symbol <- fpkm.tumor.symbol.filter.symbol.icgc[
    rownames(outlier.patient.tag.01.icgc),
    ]$Symbol;
# Match and extract rows from ICGC data based on 'five.data.outlier.symbol.na'
outlier.patient.tag.01.icgc.match.five <- outlier.patient.tag.01.icgc[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.icgc.symbol),
    ];



# Combine matched data from all sources into a single data frame
outlier.patient.all.five.01 <- data.frame(
    cbind(
        outlier.patient.tag.01.brca.match.five,
        outlier.patient.tag.01.meta.match.five,
        outlier.patient.tag.01.ispy.match.five,
        outlier.patient.tag.01.metador.match.five,
        outlier.patient.tag.01.icgc.match.five
        )
    );

# Set row names of the combined data frame
rownames(outlier.patient.all.five.01) <- five.data.outlier.symbol.na;

# Calculate the sum of non-NA values for each row
outlier.patient.all.five.01.sum <- apply(
    outlier.patient.all.five.01,
    1,
    function(x) {
        sum(na.omit(x))
        }
    );

# Calculate the fraction of the sum relative to the total number of columns
outlier.patient.all.five.01.sum.fraction <- as.numeric(outlier.patient.all.five.01.sum) / ncol(outlier.patient.all.five.01) * 100;

# Calculate the number of patients required to observe outlier gene
outlier.patient.all.five.01.sum.fraction.number.patient <- 100 / outlier.patient.all.five.01.sum.fraction;

# Create a table of the log-transformed number of patients
outlier.patient.all.five.01.sum.fraction.number.patient.table <- data.frame(
    table(
        log10(outlier.patient.all.five.01.sum.fraction.number.patient)
        )
    );

# Calculate the percentage for each frequency
outlier.patient.all.five.01.sum.fraction.number.patient.table$percent <-
    outlier.patient.all.five.01.sum.fraction.number.patient.table$Freq /
        sum(outlier.patient.all.five.01.sum.fraction.number.patient.table$Freq);

# Convert the first column to numeric
outlier.patient.all.five.01.sum.fraction.number.patient.table$Var1 <- as.numeric(
    as.vector(outlier.patient.all.five.01.sum.fraction.number.patient.table$Var1)
    );

# Smooth the line using a spline
smooth.line <- smooth.spline(
    outlier.patient.all.five.01.sum.fraction.number.patient.table$Var1,
    outlier.patient.all.five.01.sum.fraction.number.patient.table$percent * 100,
    spar = 0.7
    );

### PREDICTION ##################################################################

# Generate predicted values over a sequence of Var1
predicted.values <- data.frame(
    Var1 = seq(
        min(outlier.patient.all.five.01.sum.fraction.number.patient.table$Var1),
        max(outlier.patient.all.five.01.sum.fraction.number.patient.table$Var1),
        length.out = 100
        )
    );

# Predict y-values using the smooth line
predicted.values$predicted <- predict(smooth.line, x = predicted.values$Var1)$y;

### SMOOTHING FUNCTION ##########################################################

# Convert the smoothed spline into a function
smooth.func <- list(
    function(x) {
        predict(smooth.line, x)$y
        }
    );


### PLOTTING ####################################################################

# Create a scatter plot with the smoothed spline curve
scatter.smooth.line <- create.scatterplot(
    percent * 100 ~ Var1,
    data = outlier.patient.all.five.01.sum.fraction.number.patient.table,
    type = c('p'),
    main = expression('Number of patients needed to observe outlier gene'),
    main.cex = 1.3,
    ylab.label = expression('Percent'),
    xlab.lab = expression('Outlier Frequency'),
    xlimits = c(-0.02, 4.09),
    ylimits = c(-2.5, 59),
    xat = c(log10(1), log10(10), log10(100), log10(1000), log10(10000)),
    xaxis.lab = c(0, 10, 100, 1000, 10000),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    cex = 1.4,
    pch = 21,
    col = 'orange',
    col.border = 'black',
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    width = 10,
    height = 10,
    lwd = 3,
    abline.h = 0,
    abline.col = 'grey30',
    abline.lwd = 1.5,
    add.curves = TRUE,
    curves.exprs = smooth.func,
    curves.from = min(outlier.patient.all.five.01.sum.fraction.number.patient.table$Var1),
    curves.to = max(outlier.patient.all.five.01.sum.fraction.number.patient.table$Var1),
    curves.col = 'grey15',
    curves.lwd = 3.5,
    curves.lty = 1,
    add.grid = TRUE,
    ygrid.at = seq(0, 60, 10),
    xgrid.at = c(log10(10), log10(100), log10(1000), log10(10000)),
    grid.colour = 'grey90'
    );

# Display the plot
scatter.smooth.line;


### OUTPUT ######################################################################

save.outlier.figure(
    scatter.smooth.line,
    c('Figure1c', '5_patient_per_outlier_ratio_needed_patient_percent_smoothline', 'scatter'),
    width = 5.5,
    height = 5
    );

save.session.profile(file.path('output', 'Figure1c.txt'));
