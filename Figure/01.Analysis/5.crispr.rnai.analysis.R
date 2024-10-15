# Load necessary utility functions and helper libraries
library(BoutrosLab.utilities)
library(outlierAnalysisSupport)

# Attach the directory containing the outlier data path
attach(get.outlier.data.path())

# Load precomputed variables that will be needed for the analysis
load.multiple.computed.variables(c(
    'outlier.symbol', # List of outlier gene symbols
    'ccle.sample.outlier.status.overlap.na', # Outlier status for CCLE samples
    'gene.rnai.diff.matrix.05.overlap' # Differential RNAi matrix for overlap genes
    ))

### GENE EFFECT ANALYSIS USING CAS-CRISPR DATASET ###############################

# Extract the gene effect scores from the Cas-CRISPR dataset for the selected outlier genes
# Match the rownames to 'ccle.outlier.rank.fdr.05' to only include outlier genes with FDR < 0.05
cas.effect.breast.05 <- cas.effect.breast[rownames(ccle.outlier.rank.fdr.05), ]

# Remove any rows with missing values
cas.effect.breast.05.na <- na.omit(cas.effect.breast.05)

# Initialize lists to store quantiles, outlier gene effect scores, and non-outlier gene effect scores
effect.quantile.05 <- list()
outlier.gene.effect.score.05 <- list()
nonoutlier.gene.effect.score.05 <- list()

# Loop through each sample in the dataset
for (i in 1:nrow(ccle.sample.outlier.status.overlap.na)) {
    # Extract gene effect scores for outlier samples
    outlier.gene <- cas.effect.breast.05.na[i, which(ccle.sample.outlier.status.overlap.na[i, ] == 1), drop = FALSE]

    # Extract gene effect scores for non-outlier samples
    non.outlier.gene <- cas.effect.breast.05.na[i, -(which(ccle.sample.outlier.status.overlap.na[i, ] == 1))]

    # Compute the empirical CDF (quantile function) based on non-outlier gene effect scores
    ecdf.obj <- ecdf(as.numeric(non.outlier.gene))
    quantile.value <- t(data.frame(ecdf.obj(outlier.gene)))

    # Set column and row names for the quantile values
    colnames(quantile.value) <- colnames(outlier.gene)
    rownames(quantile.value) <- rownames(outlier.gene)

    # Store quantile values and gene effect scores in the corresponding lists
    effect.quantile.05[[i]] <- quantile.value
    outlier.gene.effect.score.05[[i]] <- outlier.gene
    nonoutlier.gene.effect.score.05[[i]] <- non.outlier.gene
    }

### MEAN GENE EFFECT SCORE CALCULATION ##########################################

# Calculate the mean gene effect score for outlier samples
outlier.gene.effect.score.05.mean <- sapply(outlier.gene.effect.score.05, function(x) {
    mean(unlist(x)) # Convert each list to a numeric vector and calculate mean
    })
outlier.gene.effect.score.05.mean <- data.frame(effect = outlier.gene.effect.score.05.mean)
rownames(outlier.gene.effect.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na)

# Calculate the mean gene effect score for non-outlier samples
nonoutlier.gene.effect.score.05.mean <- sapply(nonoutlier.gene.effect.score.05, function(x) {
    mean(unlist(x)) # Convert each list to a numeric vector and calculate mean
    })
nonoutlier.gene.effect.score.05.mean <- data.frame(effect = nonoutlier.gene.effect.score.05.mean)
rownames(nonoutlier.gene.effect.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na)

### DIFFERENCE IN GENE EFFECT SCORES ############################################

# Create a difference matrix to store the effect scores for outliers and non-outliers
# along with the calculated differences
gene.effect.diff.matrix.05 <- data.frame(
    out = outlier.gene.effect.score.05.mean$effect, # Outlier gene effect scores
    non = nonoutlier.gene.effect.score.05.mean$effect, # Non-outlier gene effect scores
    diff = outlier.gene.effect.score.05.mean$effect - nonoutlier.gene.effect.score.05.mean$effect, # Difference
    symbol = sub('\\..*', '', rownames(outlier.gene.effect.score.05.mean)) # Extract symbol without version
    )
rownames(gene.effect.diff.matrix.05) <- rownames(outlier.gene.effect.score.05.mean)

### FILTER FOR OVERLAPPING GENES ################################################

# Filter the difference matrix to only include the genes that overlap with outlier symbols
gene.effect.diff.matrix.05.overlap <- gene.effect.diff.matrix.05[
    gene.effect.diff.matrix.05$symbol %in% outlier.symbol$unique,
    ]

### SAVE RESULTS ################################################################

# Cache computed variables for further analysis in other scripts
cache.multiple.computed.variables(c(
    'cas.effect.breast.05.na', # Cleaned Cas-CRISPR gene effect scores
    'gene.effect.diff.matrix.05.overlap' # Filtered matrix for overlapping genes
    ))

# Save the session profile for reproducibility
save.session.profile(file.path('output', '5.crispr.rnai.analysis.txt'))
