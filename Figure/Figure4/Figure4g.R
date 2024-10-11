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
    'five.data.outlier.symbol',
    'ccle.sample.outlier.status.overlap'
    ));


# Filter for FDR < 0.05 and match gene names
gene.rnai.breast.t.num.match.05 <- rnai.effect.breast[
    rownames(rnai.effect.breast) %in% gsub('\\..*$', '', rownames(ccle.outlier.rank.fdr.05)),
    ];
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05;
sample.outlier.05.overlap <- ccle.sample.outlier.status.overlap
sample.outlier.05.overlap.na <- sample.outlier.05.overlap[match(rownames(gene.rnai.breast.t.num.match.05.na), gsub('\\..*$', '', rownames(sample.outlier.05.overlap))), ];
sample.outlier.05.overlap.na <- sample.outlier.05.overlap.na[, colnames(gene.rnai.breast.t.num.match.05.na)]; # should chnage the name
sample.outlier.05.overlap.na.rnai <- sample.outlier.05.overlap.na[, colnames(gene.rnai.breast.t.num.match.05.na)];
rownames(gene.rnai.breast.t.num.match.05.na) <- rownames(sample.outlier.05.overlap.na);
sample.outlier.05.overlap.na.sum <- apply(sample.outlier.05.overlap.na, 1, sum);

sample.outlier.05.overlap.na <- sample.outlier.05.overlap.na[sample.outlier.05.overlap.na.sum > 0, ];
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05.na[sample.outlier.05.overlap.na.sum > 0, ];


rnai.quantile.05 <- list();
outlier.gene.rnai.score.05 <- list();
nonoutlier.gene.rnai.score.05 <- list();
for (i in 1:nrow(sample.outlier.05.overlap.na)) {
    outlier.gene <- gene.rnai.breast.t.num.match.05.na[i, which(sample.outlier.05.overlap.na[i, ] == 1), drop = FALSE];
    non.outlier.gene <- gene.rnai.breast.t.num.match.05.na[i, -(which(sample.outlier.05.overlap.na[i, ] == 1))];

    ecdf.obj <- ecdf(as.numeric(non.outlier.gene));
    quantile.value <- ecdf.obj(outlier.gene);
    quantile.value <- t(data.frame(quantile.value));
    colnames(quantile.value) <- colnames(outlier.gene);
    rownames(quantile.value) <- rownames(outlier.gene);
    rnai.quantile.05[[i]] <- quantile.value;

    outlier.gene.rnai.score.05[[i]] <- outlier.gene;
    nonoutlier.gene.rnai.score.05[[i]] <- non.outlier.gene;
    }



# Calculate mean RNAi scores and differences
outlier.gene.rnai.score.05.mean <- sapply(outlier.gene.rnai.score.05, function(x) {
    mean(na.omit(unlist(x)));
    });
outlier.gene.rnai.score.05.mean <- data.frame(rnai = outlier.gene.rnai.score.05.mean);
rownames(outlier.gene.rnai.score.05.mean) <- rownames(sample.outlier.05.overlap.na);

nonoutlier.gene.rnai.score.05.mean <- sapply(nonoutlier.gene.rnai.score.05, function(x) {
    mean(na.omit(unlist(x)));
    });
nonoutlier.gene.rnai.score.05.mean <- data.frame(rnai = nonoutlier.gene.rnai.score.05.mean);
rownames(nonoutlier.gene.rnai.score.05.mean) <- rownames(sample.outlier.05.overlap.na);

# Create difference matrix for RNAi scores
gene.rnai.diff.matrix.05 <- data.frame(
    out = outlier.gene.rnai.score.05.mean$rnai,
    non = nonoutlier.gene.rnai.score.05.mean$rnai,
    diff = outlier.gene.rnai.score.05.mean$rnai - nonoutlier.gene.rnai.score.05.mean$rnai,
    symbol = sub('\\..*', '', rownames(outlier.gene.rnai.score.05.mean))
    );
rownames(gene.rnai.diff.matrix.05) <- rownames(outlier.gene.rnai.score.05.mean);
gene.rnai.diff.matrix.05.overlap <- gene.rnai.diff.matrix.05[gene.rnai.diff.matrix.05$symbol %in% five.data.outlier.symbol, ];

cache.multiple.computed.variables(c(
    'gene.rnai.breast.t.num.match.05.na',
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
