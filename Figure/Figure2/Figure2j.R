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
library(stats)
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

# Load required variables
load.multiple.computed.variables(c(
    'brca.outlier.non.promoter.symbol.sample.match.merge.500',
    'me.out.symbol.two.500',
    'meta.outlier.non.promoter.symbol.sample.match.merge.500',
    'non.outlier.sample.me.two.500',
    'outlier.sample.me.two.500',
    'two.outlier.patient.status.merge.filter.500',
    'two.outlier.promoter.symbol.sample.match.merge.filter.500'
    ))

# Define constants
unequal.quan <- rev(seq(0, 0.9, 0.1))

# Function to analyze beta values
analyze_beta_values <- function(sample_list, data_matrix) {
    lapply(seq_len(nrow(data_matrix)), function(i) {
        value.vector <- as.numeric(data_matrix[i, ])
        percentile <- ecdf(value.vector)
        percent.value <- percentile(sample_list[[i]])
        na.omit(percent.value)
        }) |> unlist()
    }

# Compute the sum of the patient status filter
two.outlier.patient.status.merge.filter.sum.500 <- apply(
    two.outlier.patient.status.merge.filter.500,
    2,
    function(x) {
        sum(na.omit(as.numeric(x)))
        }
    )

# Analyze outlier and non-outlier probes in outlier and non-outlier patients
percent.beta.out.out.500 <- analyze_beta_values(outlier.sample.me.two.500, two.outlier.promoter.symbol.sample.match.merge.filter.500)
percent.beta.non.out.500 <- analyze_beta_values(non.outlier.sample.me.two.500, two.outlier.promoter.symbol.sample.match.merge.filter.500)

# Analyze non-outlier probes
me.non.out.symbol.two.500 <- setdiff(
    unique(c(rownames(brca.outlier.non.promoter.symbol.sample.match.merge.500), rownames(meta.outlier.non.promoter.symbol.sample.match.merge.500))),
    me.out.symbol.two.500
    )

# Analyze non-outlier probes in outlier and non-outlier patients
analyze_non_outlier_beta <- function(gene_symbol, brca_data, meta_data, patient_status_sum) {
    value.vector.brca <- as.numeric(brca_data[rownames(brca_data) == gene_symbol, ])
    value.vector.meta <- as.numeric(meta_data[rownames(meta_data) == gene_symbol, ])
    value.vector <- na.omit(c(value.vector.brca, value.vector.meta))

    if (length(value.vector) == 0) {
        return(list(out.non = NA, non.non = NA))
        }

    percentile <- ecdf(value.vector)
    list(
        out.non = na.omit(percentile(value.vector[patient_status_sum > 0])),
        non.non = na.omit(percentile(value.vector[patient_status_sum == 0]))
        )
    }

# Compute beta values for non-outlier probes
non.outlier.results <- lapply(me.non.out.symbol.two.500, function(gene) {
    analyze_non_outlier_beta(gene, brca.outlier.non.promoter.symbol.sample.match.merge.500, meta.outlier.non.promoter.symbol.sample.match.merge.500, two.outlier.patient.status.merge.filter.sum.500)
    })

# Unlist results
percent.beta.out.non.500 <- unlist(lapply(non.outlier.results, `[[`, 'out.non'))
percent.beta.non.non.500 <- unlist(lapply(non.outlier.results, `[[`, 'non.non'))

# Function to create histograms and calculate percentages
create_histogram_df <- function(data, breaks = seq(0, 1, length.out = 31)) {
    hist.data <- hist(data, breaks = breaks, plot = FALSE)
    percentages <- hist.data$counts / sum(hist.data$counts) * 100
    data.frame(
        bin_start = hist.data$breaks[-length(hist.data$breaks)],
        bin_end = hist.data$breaks[-1],
        percentage = percentages
        )
    }

# Calculate percentages for each group
percent.beta.out.out.500.percentages.df <- create_histogram_df(percent.beta.out.out.500)
percent.beta.non.out.500.percentages.df <- create_histogram_df(percent.beta.non.out.500)
percent.beta.out.non.500.percentages.df <- create_histogram_df(percent.beta.out.non.500)
percent.beta.non.non.500.percentages.df <- create_histogram_df(percent.beta.non.non.500)

# Combine percentages for all groups
percent.merge.two.four.group.500 <- cbind(
    group1 = percent.beta.out.out.500.percentages.df$percentage,
    group2 = percent.beta.non.out.500.percentages.df$percentage,
    group3 = percent.beta.out.non.500.percentages.df$percentage,
    group4 = percent.beta.non.non.500.percentages.df$percentage
    )

# Prepare data for heatmap
heat.df <- data.frame(percent.merge.two.four.group.500)
heat.df.rev <- heat.df[rev(seq(nrow(heat.df))), ]

# Set up heatmap legend
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
    )

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
    )

# Save the heatmap
save.outlier.figure(heat.out, c('Figure2j', 'merge', 'me', 'quantile', 'heatmap'), width = 6, height = 4.5)
save.session.profile(file.path('output', 'Figure2j.txt'))
