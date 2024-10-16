### HISTORY ######################################################################
# This script processes and merges DNA methylation data.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script analyzes DNA methylation data for breast cancer. It merges data
# from multiple sources, processes outlier statuses, and creates visualizations
# of DNA methylation differences between outlier and non-outlier samples. The script
# focuses on promoter regions (TSS ~ +500bp) and includes statistical analysis.

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'mean.minus.ma.merge.two.500',
    'p.me'
    ));

# Delta beta > 0.2
dot.colours <- vector(
    length = length(mean.minus.ma.merge.two.500$minus.beta)
    );

dot.colours <- rep(
    'grey60',
    length(mean.minus.ma.merge.two.500$minus.beta)
    );

dot.colours[mean.minus.ma.merge.two.500$minus.beta < -0.2] <- 'dodgerblue3';
dot.colours[mean.minus.ma.merge.two.500$minus.beta > 0.2] <- 'red2';


text.pvalue <- display.statistical.result(
    x = p.me$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

key.minus.beta <- list(
    text = list(
        lab = text.pvalue,
        cex = 1
        ),
    x = 0.25,
    y = 0.95
    );

me.merge.scatter <- create.scatterplot(
    formula = minus.beta ~ mean.beta,
    data = mean.minus.ma.merge.two.500,
    col = dot.colours,
    alpha = .80,
    ylimits = c(-1.01, 1.01),
    xlimits = c(-0.01, 1.01),
    yat = seq(-1.0, 1.0, 0.5),
    xat = seq(0.0, 1.0, 0.2),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = TRUE,
    grid.colour = 'grey80',
    cex = 0.9,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1.1,
    main.cex = 1.5,
    left.padding = 0,
    main = expression('DNA methylation on promoter of outlier patients'),
    xlab.label = expression(paste('Mean of ', beta, ' value')),
    ylab.label = expression(paste(Delta, ' Mean of ', beta, ' value')),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = key.minus.beta
                ),
            x = 0.71,
            y = 0.95,
            corner = c(0, 1)
            )
        ),
    abline.h = 0,
    abline.col = 'black',
    abline.lwd = 2,
    abline.lty = 3
    );

save.outlier.figure(
    me.merge.scatter,
    c('Figure2h', 'methylation', 'diff', 'merge', 'scatter'),
    width = 6,
    height = 4.8
    );

save.session.profile(file.path('output', 'Figure2h.txt'));
