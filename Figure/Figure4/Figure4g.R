### HISTORY #####################################################################
# This script analyzes RNAi gene dependency scores for outlier genes identified
# in CCLE.
# Date: 2024-08-16

### DESCRIPTION #################################################################
# The script processes RNAi gene dependency scores from CCLE data. It filters genes
# with FDR < 0.05, matches the gene names with those in the RNAi dataset, and calculates
# the mean RNAi scores for outliers and non-outliers. It then computes the difference
# between the RNAi scores of outliers and non-outliers. The results are visualized
# with a scatter plot, highlighting interesting genes based on specific thresholds.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'sample.outlier.05.overlap.na',
    'gene.rnai.diff.matrix.05.overlap'
    ));

# Set colors for the scatter plot

interesting.points <- gene.rnai.diff.matrix.05.overlap$diff < -0.5;
text.x <- na.omit(gene.rnai.diff.matrix.05.overlap$non[interesting.points]);
text.y <- na.omit(gene.rnai.diff.matrix.05.overlap$out[interesting.points]);
text.labels <- na.omit(gene.rnai.diff.matrix.05.overlap$symbol[interesting.points]);
dot.colours <- rep('grey30', nrow(gene.rnai.diff.matrix.05.overlap));
dot.colours[gene.rnai.diff.matrix.05.overlap$diff < -0.5] <- 'dodgerblue3';
dot.colours[gene.rnai.diff.matrix.05.overlap$diff > 0.5] <- 'red2';

# Create the scatter plot
gene.scatter.05.minus.overlap.label <- create.scatterplot(
    formula = diff ~ non,
    data = gene.rnai.diff.matrix.05.overlap,
    col = dot.colours,
    alpha = .6,
    ylimits = c(-1.3, 1.05),
    xlimits = c(-1.4, 0.5),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = FALSE,
    cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1.2,
    main.cex = 1.5,
    left.padding = 0,
    main = expression('Gene effect score of the overlap outliers'),
    xlab.label = expression(paste('Mean of gene rnai score of non-outlier samples')),
    ylab.label = expression('Difference of gene effect score'),
    add.text = TRUE,
    text.x = text.x,
    text.y = text.y,
    text.labels = text.labels,
    text.fontface = 1,
    abline.h = c(0, -0.5),
    abline.col = c('grey30', 'dodgerblue3'),
    abline.lwd = c(1.5, 2),
    abline.lty = c(1, 3)
    );


save.outlier.figure(
    gene.scatter.05.minus.overlap.label,
    c('Figure4g', 'gene', 'dependency', 'diff', 'rnai', 'scatter'),
    width = 6,
    height = 5
    );

save.session.profile(file.path('output', 'Figure4g.txt'));
