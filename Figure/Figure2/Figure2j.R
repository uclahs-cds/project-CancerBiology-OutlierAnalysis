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


# Ensure that the data matrix is numeric
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- apply(two.outlier.promoter.symbol.sample.match.merge.filter.500, 2, as.numeric)

# Function to analyze beta values using lapply (with proper indexing)
analyze_beta_values <- function(sample_list, data_matrix) {
    lapply(seq_len(nrow(data_matrix)), function(i) {
        value.vector <- data_matrix[i, ] # Retrieve the ith row
        percentile <- ecdf(value.vector) # Compute the ecdf for the row
        percent_value <- percentile(sample_list[[i]]) # Apply ecdf to corresponding sample
        na.omit(percent_value) # Remove any NA values
        }) |> unlist()
    }

# Compute the sum of the patient status filter
two.outlier.patient.status.merge.filter.sum.500 <- colSums(
    apply(two.outlier.patient.status.merge.filter.500, 2, as.numeric),
    na.rm = TRUE
    )

# Analyze outlier and non-outlier probes in outlier and non-outlier patients
percent.beta.out.out.500 <- analyze_beta_values(outlier.sample.me.two.500, two.outlier.promoter.symbol.sample.match.merge.filter.500)
percent.beta.non.out.500 <- analyze_beta_values(non.outlier.sample.me.two.500, two.outlier.promoter.symbol.sample.match.merge.filter.500)

# Analyze non-outlier probes
me.non.out.symbol.two.500 <- setdiff(
    unique(c(rownames(brca.outlier.non.promoter.symbol.sample.match.merge.500), rownames(meta.outlier.non.promoter.symbol.sample.match.merge.500))),
    me.out.symbol.two.500
    )

# Precompute row indices for brca and meta data
brca_gene_indices <- match(me.non.out.symbol.two.500, rownames(brca.outlier.non.promoter.symbol.sample.match.merge.500))
meta_gene_indices <- match(me.non.out.symbol.two.500, rownames(meta.outlier.non.promoter.symbol.sample.match.merge.500))

analyze_non_outlier_beta <- function(gene_index_brca, gene_index_meta, brca_data, meta_data, patient_status_sum) {
    # Pre-extracted rows based on gene indices (instead of searching by rownames)
    value.vector.brca <- as.numeric(brca_data[gene_index_brca, , drop = FALSE])
    value.vector.meta <- as.numeric(meta_data[gene_index_meta, , drop = FALSE])
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

# Compute beta values for non-outlier probes using precomputed row indices
non.outlier.results <- mapply(
    analyze_non_outlier_beta,
    brca_gene_indices,
    meta_gene_indices,
    MoreArgs = list(
        brca_data = brca.outlier.non.promoter.symbol.sample.match.merge.500,
        meta_data = meta.outlier.non.promoter.symbol.sample.match.merge.500,
        patient_status_sum = two.outlier.patient.status.merge.filter.sum.500
        ),
    SIMPLIFY = FALSE
    )

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
