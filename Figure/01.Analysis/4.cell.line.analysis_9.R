library(BoutrosLab.utilities);

# Source the helper library for outlier analysis
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
# Attach the outlier data path
attach(get.outlier.data.path());

# Load precomputed variables, including outlier symbols
load.multiple.computed.variables(c(
    'outlier.symbol'
    ));

# Transpose gene dependency data for breast cancer and convert to numeric
gene.dependency.breast.t <- t(gene.dependency.breast)
gene.dependency.breast.t.num.match <- as.data.frame(apply(gene.dependency.breast.t, 2, as.numeric))

# Set row and column names for the numeric data frame
rownames(gene.dependency.breast.t.num.match) <- rownames(gene.dependency.breast.t)
colnames(gene.dependency.breast.t.num.match) <- colnames(gene.dependency.breast.t)

# Filter the gene dependency data to match the sample names in the CCLE dataset
gene.dependency.breast.t.num.match <- gene.dependency.breast.t.num.match[, colnames(fpkm.tumor.symbol.filter.ccle)]

# Match the outlier rank FDR < 0.05 genes with the gene dependency data
gene.dependency.breast.t.num.match.05 <- gene.dependency.breast.t.num.match[rownames(ccle.outlier.rank.fdr.05), ]
gene.dependency.breast.t.num.match.05.na <- na.omit(gene.dependency.breast.t.num.match.05)

# Find the overlapping outlier status data
ccle.sample.outlier.status.overlap <- ccle.sample.outlier.status[rownames(ccle.outlier.rank.fdr.05), ]
ccle.sample.outlier.status.overlap.na <- ccle.sample.outlier.status.overlap[rownames(gene.dependency.breast.t.num.match.05.na), ]

ccle.sample.outlier.status.fdr.05 <- ccle.sample.outlier.status[rownames(ccle.outlier.rank.fdr.05), ];
ccle.sample.outlier.status.fdr.05.five <- ccle.sample.outlier.status.fdr.05[sub('\\..*', '', rownames(ccle.sample.outlier.status.fdr.05)) %in% outlier.symbol$unique, ];
ccle.sample.outlier.status.fdr.05.five.symbol <- sub('\\..*', '', rownames(ccle.sample.outlier.status.fdr.05.five));

# Match gene names with RNAi data and filter for FDR < 0.05
gene.rnai.breast.t.num.match.05 <- rnai.effect.breast[
    rownames(rnai.effect.breast) %in% gsub('\\..*$', '', rownames(ccle.outlier.rank.fdr.05)),
    ]

# Ensure no missing data in the RNAi dataset
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05

# Filter overlapping outlier statuses
sample.outlier.05.overlap.na <- ccle.sample.outlier.status.overlap[match(rownames(gene.rnai.breast.t.num.match.05.na), gsub('\\..*$', '', rownames(ccle.sample.outlier.status.overlap))), ]
sample.outlier.05.overlap.na <- sample.outlier.05.overlap.na[, colnames(gene.rnai.breast.t.num.match.05.na)]

# Rename rownames of the RNAi dataset to match the sample overlap
rownames(gene.rnai.breast.t.num.match.05.na) <- rownames(sample.outlier.05.overlap.na)

# Calculate the sum of outliers for each sample and filter based on those with outliers
sample.outlier.05.overlap.na.sum <- apply(sample.outlier.05.overlap.na, 1, sum)
sample.outlier.05.overlap.na <- sample.outlier.05.overlap.na[sample.outlier.05.overlap.na.sum > 0, ]
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05.na[sample.outlier.05.overlap.na.sum > 0, ]

### RNAi Quantile Calculations ###################################################

# Initialize lists to store quantiles, RNAi scores for outliers, and non-outliers
rnai.quantile.05 <- list()
outlier.gene.rnai.score.05 <- list()
nonoutlier.gene.rnai.score.05 <- list()

# Loop over the RNAi data to calculate the ECDF (quantile) for outlier genes
for (i in 1:nrow(sample.outlier.05.overlap.na)) {
    outlier.gene <- gene.rnai.breast.t.num.match.05.na[i, which(sample.outlier.05.overlap.na[i, ] == 1), drop = FALSE]
    non.outlier.gene <- gene.rnai.breast.t.num.match.05.na[i, -(which(sample.outlier.05.overlap.na[i, ] == 1))]

    ecdf.obj <- ecdf(as.numeric(non.outlier.gene)) # Compute empirical CDF for non-outliers
    quantile.value <- t(data.frame(ecdf.obj(outlier.gene))) # Get quantile values for outliers

    # Set column and row names for quantile values
    colnames(quantile.value) <- colnames(outlier.gene)
    rownames(quantile.value) <- rownames(outlier.gene)

    # Store calculated values in the lists
    rnai.quantile.05[[i]] <- quantile.value
    outlier.gene.rnai.score.05[[i]] <- outlier.gene
    nonoutlier.gene.rnai.score.05[[i]] <- non.outlier.gene
    }

### RNAi Mean and Differences Calculation ########################################

# Calculate mean RNAi scores for outliers and non-outliers
outlier.gene.rnai.score.05.mean <- sapply(outlier.gene.rnai.score.05, function(x) mean(na.omit(unlist(x))))
outlier.gene.rnai.score.05.mean <- data.frame(rnai = outlier.gene.rnai.score.05.mean)
rownames(outlier.gene.rnai.score.05.mean) <- rownames(sample.outlier.05.overlap.na)

nonoutlier.gene.rnai.score.05.mean <- sapply(nonoutlier.gene.rnai.score.05, function(x) mean(na.omit(unlist(x))))
nonoutlier.gene.rnai.score.05.mean <- data.frame(rnai = nonoutlier.gene.rnai.score.05.mean)
rownames(nonoutlier.gene.rnai.score.05.mean) <- rownames(sample.outlier.05.overlap.na)

# Create a difference matrix to store outlier and non-outlier RNAi score differences
gene.rnai.diff.matrix.05 <- data.frame(
    out = outlier.gene.rnai.score.05.mean$rnai,
    non = nonoutlier.gene.rnai.score.05.mean$rnai,
    diff = outlier.gene.rnai.score.05.mean$rnai - nonoutlier.gene.rnai.score.05.mean$rnai,
    symbol = sub('\\..*', '', rownames(outlier.gene.rnai.score.05.mean))
    )
rownames(gene.rnai.diff.matrix.05) <- rownames(outlier.gene.rnai.score.05.mean)

# Filter genes overlapping with the outlier symbols
gene.rnai.diff.matrix.05.overlap <- gene.rnai.diff.matrix.05[gene.rnai.diff.matrix.05$symbol %in% outlier.symbol$unique, ]

### RNAi Score Overlap and Difference Calculation ################################

# Filter for genes with a significant difference (e.g., difference < -0.4)
gene.rnai.diff.matrix.05.overlap.minus.05 <- gene.rnai.diff.matrix.05.overlap[gene.rnai.diff.matrix.05.overlap$diff < -0.4, ]
gene.rnai.diff.matrix.05.overlap.minus.05 <- na.omit(gene.rnai.diff.matrix.05.overlap.minus.05)

# Order the genes by their difference values
gene.rnai.diff.matrix.05.overlap.minus.05.order <- gene.rnai.diff.matrix.05.overlap.minus.05[order(gene.rnai.diff.matrix.05.overlap.minus.05$diff), ]

# Extract the corresponding RNAi scores and outlier statuses
rnai.score.05.overlap.minus.05 <- gene.rnai.breast.t.num.match.05.na[rownames(gene.rnai.diff.matrix.05.overlap.minus.05.order), ]
outlier.status.05.overlap.minus.05 <- ccle.sample.outlier.status.overlap.na[rownames(gene.rnai.diff.matrix.05.overlap.minus.05.order), colnames(rnai.score.05.overlap.minus.05)]

# Prepare data for boxplot creation by replicating gene names for each column
gene.name.box.05 <- unlist(lapply(1:nrow(outlier.status.05.overlap.minus.05), function(i) {
    rep(sub('\\..*', '', rownames(rnai.score.05.overlap.minus.05))[i], ncol(rnai.score.05.overlap.minus.05))
    }))

# Combine RNAi score and status into a data frame for plotting
rnai.05.box <- data.frame(
    score = as.numeric(unlist(t(rnai.score.05.overlap.minus.05))),
    gene = gene.name.box.05,
    status = as.numeric(unlist(t(outlier.status.05.overlap.minus.05)))
    )
rnai.05.box$score <- as.numeric(rnai.05.box$score)
rnai.05.box$status <- as.numeric(rnai.05.box$status)

### Protein Data Matching ########################################################

# Process protein data and match gene symbols to the protein dataset
protein.info.breast.num.symbol <- sapply(strsplit(rownames(protein.info.breast.num), '\\|'), function(x) x[3])
protein.info.breast.num.symbol <- sub('_HUMAN', '', protein.info.breast.num.symbol)

# Match protein data to CCLE sample columns
protein.info.breast.num.match <- protein.info.breast.num[, colnames(protein.info.breast.num) %in% colnames(fpkm.tumor.symbol.filter.ccle)]
rownames(protein.info.breast.num.match) <- protein.info.breast.num.symbol

# Filter protein data based on outlier status
protein.info.breast.num.match.05 <- protein.info.breast.num.match[rownames(protein.info.breast.num.match) %in% ccle.sample.outlier.status.fdr.05.five.symbol, ]

### Caching and Saving Session Data ##############################################

# Cache computed variables for future use
cache.multiple.computed.variables(c(
    'ccle.sample.outlier.status.fdr.05.five',
    'ccle.sample.outlier.status.fdr.05.five.symbol',
    'ccle.sample.outlier.status.overlap',
    'ccle.sample.outlier.status.overlap.na',
    'gene.dependency.breast.t.num.match.05.na',
    'gene.rnai.breast.t.num.match.05.na',
    'gene.rnai.diff.matrix.05.overlap',
    'protein.info.breast.num.match',
    'protein.info.breast.num.match.05',
    'protein.info.breast.num.symbol',
    'rnai.05.box',
    'rnai.score.05.overlap.minus.05',
    'sample.outlier.05.overlap.na'
    ))

# Save the session profile
save.session.profile(file.path('output', '4.cell.line.analysis.txt'))
