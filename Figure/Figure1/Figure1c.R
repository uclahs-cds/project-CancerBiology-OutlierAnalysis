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
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

# Save these variables for later scripts
load.multiple.computed.variables(c(
    'outlier.symbol',
    'outlier.gene.fdr.01',
    'outlier.patient.all.five.01'
    ));

# Calculate the sum of non-NA values for each row
outlier.patient.all.five.01.sum <- apply(
    outlier.patient.all.five.01,
    1,
    function(x) {
        sum(na.omit(x))
        }
    );

# Calculate the fraction of the sum relative to the total number of columns
outlier.patient.all.five.01.sum.fraction <- (as.numeric(outlier.patient.all.five.01.sum) / ncol(outlier.patient.all.five.01)) * 100;

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
