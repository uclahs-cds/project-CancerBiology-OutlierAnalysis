### HISTORY #####################################################################
# This script processes methylation data to analyze outlier and non-outlier
# genes in outlier and non-outlier patients.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script analyzes methylation data for outlier and non-outlier genes in
# different patient groups. It creates histograms, calculates percentages,
# and generates a heatmap visualization of the DNA methylation patterns across
# different gene and patient categories.

### PREAMBLE #####################################################################
# Load necessary libraries
library(stats);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-10-03_Figure1_2_3_4_min_input.rda'));

load.multiple.computed.variables(c(
    'brca.outlier.non.promoter.symbol.sample.match.merge.500',
    'me.out.symbol.two.500',
    'meta.outlier.non.promoter.symbol.sample.match.merge.500',
    'non.outlier.sample.me.two.500',
    'outlier.sample.me.two.500',
    'two.outlier.patient.status.merge.filter.500',
    'two.outlier.promoter.symbol.sample.match.merge.filter.500'
    ));


unequal.quan <- rev(seq(0, 0.9, 0.1));

# 1. Analyze outlier probes in outlier patients
percent.beta.out.out.500 <- NULL;
for (i in 1:nrow(two.outlier.promoter.symbol.sample.match.merge.filter.500)) {
    value.vector <- as.numeric(two.outlier.promoter.symbol.sample.match.merge.filter.500[i,]);
    percentile <- ecdf(value.vector);
    percent.value <- percentile(outlier.sample.me.two.500[[i]]);
    percent.value.mean <- na.omit(percent.value);
    percent.beta.out.out.500 <- c(percent.beta.out.out.500, percent.value.mean);
    }

# 2. Analyze outlier probes in non-outlier patients
percent.beta.non.out.500 <- NULL;
for (i in 1:nrow(two.outlier.promoter.symbol.sample.match.merge.filter.500)) {
    value.vector <- as.numeric(two.outlier.promoter.symbol.sample.match.merge.filter.500[i,]);
    percentile <- ecdf(value.vector);
    percent.value <- percentile(non.outlier.sample.me.two.500[[i]]);
    percent.value.mean <- na.omit(percent.value);
    percent.beta.non.out.500 <- c(percent.beta.non.out.500, percent.value.mean);
    }

# 3. & 4. Analyze non-outlier probes in outlier and non-outlier patients
two.outlier.patient.status.merge.filter.sum.500 <- apply(
    two.outlier.patient.status.merge.filter.500, 
    2, 
    function(x) { sum(na.omit(as.numeric(x))); }
    );

me.non.out.symbol.two.500 <- unique(c(
    rownames(brca.outlier.non.promoter.symbol.sample.match.merge.500), 
    rownames(meta.outlier.non.promoter.symbol.sample.match.merge.500)
    ));

me.non.out.symbol.two.500 <- me.non.out.symbol.two.500[
    !(me.non.out.symbol.two.500 %in% me.out.symbol.two.500)
    ];

brca.outlier.non.promoter.symbol.sample.match.merge.500 <- data.frame(
    brca.outlier.non.promoter.symbol.sample.match.merge.500
    );

meta.outlier.non.promoter.symbol.sample.match.merge.500 <- data.frame(
    meta.outlier.non.promoter.symbol.sample.match.merge.500
    );


percent.beta.out.non.500 <- list();
percent.beta.non.non.500 <- list();

for (i in 1:length(me.non.out.symbol.two.500)) {
    i.symbol <- me.non.out.symbol.two.500[i];
    value.vector.brca <- as.numeric(
        brca.outlier.non.promoter.symbol.sample.match.merge.500[rownames(brca.outlier.non.promoter.symbol.sample.match.merge.500) %in% i.symbol,]
        );
    value.vector.meta <- as.numeric(
        meta.outlier.non.promoter.symbol.sample.match.merge.500[rownames(meta.outlier.non.promoter.symbol.sample.match.merge.500) %in% i.symbol,]
        );
    
    value.vector <- na.omit(c(value.vector.brca, value.vector.meta));
    if (length(value.vector) == 0) {
        percent.beta.out.non.500[[i]] <- NA;
        percent.beta.non.non.500[[i]] <- NA;
        }
    else {
        percentile <- ecdf(value.vector);
        
        percent.value.out.non <- percentile(
            value.vector[two.outlier.patient.status.merge.filter.sum.500 > 0]
            );
        percent.beta.out.non.500[[i]] <- na.omit(percent.value.out.non);
        
        percent.value.non.non <- percentile(
            value.vector[two.outlier.patient.status.merge.filter.sum.500 == 0]
            );
        percent.beta.non.non.500[[i]] <- na.omit(percent.value.non.non);
        }
    }

percent.beta.out.non.500 <- unlist(percent.beta.out.non.500);
percent.beta.non.non.500 <- unlist(percent.beta.non.non.500);


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
