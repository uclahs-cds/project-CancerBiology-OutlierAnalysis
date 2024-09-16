### HISTORY #####################################################################
# This script processes methylation data to analyze outlier and non-outlier
# genes in  outier and non-outlier patients.
# Date: 2024-08-14

# Load necessary libraries
library(stats);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
args <- commandArgs();
source(file.path(
    dirname(dirname(normalizePath(sub('^--file=', '', args[grep('^--file=', args)])))),
    'common_functions.R'
    ));
# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-08-26_Figure2h-l_input.rda'));

# Create histograms and calculate percentages for each group
breaks <- seq(0, 1, length.out = 31);

create.histogram.df <- function(data) {
    hist.data <- hist(data, breaks = breaks, plot = FALSE)
    percentages <- hist.data$counts / sum(hist.data$counts) * 100
    data.frame(
        bin_start = hist.data$breaks[-length(hist.data$breaks)],
        bin_end = hist.data$breaks[-1],
        percentage = percentages
        );
    }

percent.beta.out.out.500.percentages.df <- create.histogram.df(percent.beta.out.out.500);
percent.beta.non.out.500.percentages.df <- create.histogram.df(percent.beta.non.out.500);
percent.beta.out.non.500.percentages.df <- create.histogram.df(percent.beta.out.non.500);
percent.beta.non.non.500.percentages.df <- create.histogram.df(percent.beta.non.non.500);

# Combine percentages for all groups
percent.merge.two.four.group.500 <- cbind(
    group1 = percent.beta.out.out.500.percentages.df$percentage,
    group2 = percent.beta.non.out.500.percentages.df$percentage,
    group3 = percent.beta.out.non.500.percentages.df$percentage,
    group4 = percent.beta.non.non.500.percentages.df$percentage
    );

# Prepare data for heatmap
heat.df <- data.frame(percent.merge.two.four.group.500);
heat.df.rev <- heat.df[rev(seq(nrow(heat.df))), ];

# Set up legend
legend.col <- list(
    legend = list(
        colours = c('#2166ac', 'white', '#b2182b'),
        title = expression(underline('Percentage')),
        labels = c(0, 10),
        size = 3,
        label.cex = 1,
        continuous = TRUE,
        height = 3
        )
    );

# Create heatmap
heat.out <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(heat.df.rev),
    clustering.method = 'none',
    colour.scheme = c('#2166ac', 'white', '#b2182b'),
    col.colour = 'white',
    grid.row = FALSE,
    grid.col = TRUE,
    yaxis.tck = 0,
    xaxis.tck = 0,
    xaxis.lab = NULL,
    xlab.label = expression('Outlier Genes'),
    yaxis.cex = 1.2,
    xaxis.cex = 1.2,
    yaxis.rot = 0,
    xaxis.rot = 90,
    xlab.cex = 1.2,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    colour.centering.value = 5,
    at = seq(0, 10, 0.01),
    covariate.legend = legend.col,
    legend.cex = 1,
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );


save.outlier.figure(
    heat.out,
    c('Figure2j', 'merge', 'me', 'quantile', 'heatmap'),
    width = 6,
    height = 4.5
    );

save.session.profile(file.path('output', 'Figure2j.txt'));
