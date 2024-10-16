### HISTORY #####################################################################
# This script analyzes RNAi gene dependency scores for outlier genes identified
# in CCLE.
# Date: 2024-08-16

### DESCRIPTION #################################################################
# The script processes RNAi gene dependency scores for outliers from CCLE data.
# It filters out genes with a difference in RNAi score < -0.4, orders them based
# on the score, and prepares the data for plotting. The resulting plot shows gene
# effect scores for outliers, with highlighted genes that meet specific criteria.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'rnai.05.box'
    ));

dot.colours <- vector(length = nrow(rnai.05.box));
dot.colours <- rep('grey70', nrow(rnai.05.box));
dot.colours[rnai.05.box$status == 1] <- 'dodgerblue2';

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure4h')));

rnai.05.box.plot <- BoutrosLab.plotting.general::create.boxplot(
    formula = score ~ gene,
    data = rnai.05.box,
    main = expression('Gene effect score of outlier genes'),
    outlier = TRUE,
    add.stripplot = TRUE,
    add.rectangle = TRUE,
    xleft.rectangle = seq(0.5, 14.5, 2),
    xright.rectangle = seq(1.5, 15.5, 2),
    ybottom.rectangle = -3,
    ytop.rectangle = 5,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    main.cex = 1.5,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Gene effect score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.rot = 90,
    sample.order = 'decreasing',
    points.pch = 16,
    points.cex = 1,
    points.col = dot.colours,
    lwd = 1.2,
    col = c('gold2'),
    alpha = 0.25
    );

save.outlier.figure(
    rnai.05.box.plot,
    c('Figure4h', 'gene', 'effect', 'example', 'rnai', 'box'),
    width = 6,
    height = 6
    );

save.session.profile(file.path('output', 'Figure4h.txt'));
