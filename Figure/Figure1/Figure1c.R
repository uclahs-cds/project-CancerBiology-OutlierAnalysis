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

### DATA PREPARATION ############################################################

# Create a summary table of the overlap data
outlier.patient.overlap.summary <- data.frame(
    table(outlier.patient.all.five.01.overlap.sum)
    );

# Rename the columns
colnames(outlier.patient.overlap.summary) <- c('patient', 'number');

# Add extrapolated values and their logarithmic transformation
outlier.patient.overlap.summary$ext <- 100000 / as.numeric(outlier.patient.overlap.summary$patient);
outlier.patient.overlap.summary$ext.log <- log10(outlier.patient.overlap.summary$ext);
outlier.patient.overlap.summary$fraction <- outlier.patient.overlap.summary$number / sum(outlier.patient.overlap.summary$number);

# Fit a smoothing spline to the data
smooth.line.ext <- smooth.spline(
    outlier.patient.overlap.summary$ext.log, 
    outlier.patient.overlap.summary$fraction * 100, 
    spar = 0.7
    );

### PREDICTION ##################################################################

# Generate predicted values for the extrapolated line
predicted.values.ext <- data.frame(
    Var1 = seq(
        min(outlier.patient.overlap.summary$ext.log), 
        max(outlier.patient.overlap.summary$ext.log), 
        length.out = 100
        )
    );
predicted.values.ext$predicted <- predict(smooth.line.ext, x = predicted.values.ext$Var1)$y;

### SMOOTHING FUNCTION ##########################################################

# Convert the smoothed spline into a function
smooth.func.ext <- list(function(x) {
    predict(smooth.line.ext, x)$y;
    });

### PLOTTING ####################################################################

# Create a scatter plot with the smoothed spline curve
scatter.smooth.line.ext <- create.scatterplot(
    formula = fraction * 100 ~ ext.log,
    data = outlier.patient.overlap.summary,
    type = c("p"),
    main = expression('Number of patients needed to observe outlier gene'),
    main.cex = 1.3,
    ylab.label = expression('Percent'),
    xlab.lab = expression('Required sample size'),
    xlimits = c(-0.02, 5.4),
    ylimits = c(-2.5, 65),
    xat = c(log10(1), log10(10), log10(100), log10(1000), log10(10000), log10(100000)),
    xaxis.lab = c(0, 10, 100, 1000, 10000, 100000),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    cex = 1.2,
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
    curves.exprs = smooth.func.ext,
    curves.from = min(outlier.patient.overlap.summary$ext.log),
    curves.to = max(outlier.patient.overlap.summary$ext.log),
    curves.col = 'grey15',
    curves.lwd = 3.5,
    curves.lty = 1,
    add.grid = TRUE,
    ygrid.at = seq(0, 60, 10),
    xgrid.at = c(log10(10), log10(100), log10(1000), log10(10000)),
    grid.colour = 'grey90'
    );

# Display the plot
scatter.smooth.line.ext;

### OUTPUT ######################################################################

# Save the plot as a PDF
pdf(
    file = generate.filename(
        '5_patient_per_outlier_ratio_needed_patient_percent_smoothline_fromzero', 
        'scatter', 
        'pdf'
        ), 
    width = 5.5, 
    height = 5
    );
scatter.smooth.line.ext;
dev.off();

# Save the plot as a PNG
png(
    file = generate.filename(
        '5_patient_per_outlier_ratio_needed_patient_percent_smoothline_fromzero', 
        'scatter', 
        'png'
        ), 
    width = 5.5, 
    height = 5, 
    unit = 'in', 
    res = 1200
    );
scatter.smooth.line.ext;
dev.off();
