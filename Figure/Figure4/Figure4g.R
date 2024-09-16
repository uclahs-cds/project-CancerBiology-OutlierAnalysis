### HISTORY #####################################################################
# This script analyzes RNAi gene dependency scores for outlier genes identified
# in CCLE.
# Date: 2024-08-16

library(BoutrosLab.plotting.general);

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));


# Filter for FDR < 0.05 and match gene names
gene.rnai.breast.t.num.match.05 <- rnai.effect.breast[
    rownames(rnai.effect.breast) %in% gsub('\\..*$', '', rownames(ccle.outlier.rank.fdr.05)),
    ];
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05;
ccle.sample.outlier.status.overlap <- ccle.sample.outlier.status[rownames(ccle.outlier.rank.fdr.05), ];
ccle.sample.outlier.status.overlap.na <- ccle.sample.outlier.status.overlap[
    match(rownames(gene.rnai.breast.t.num.match.05.na), gsub('\\..*$', '', rownames(ccle.sample.outlier.status.overlap))),
    ];
ccle.sample.outlier.status.overlap.na <- ccle.sample.outlier.status.overlap.na[, colnames(gene.rnai.breast.t.num.match.05.na)];
rownames(gene.rnai.breast.t.num.match.05.na) <- rownames(ccle.sample.outlier.status.overlap.na);

# Sum and filter out genes with no outliers
ccle.sample.outlier.status.overlap.na.sum <- apply(ccle.sample.outlier.status.overlap.na, 1, sum);
ccle.sample.outlier.status.overlap.na <- ccle.sample.outlier.status.overlap.na[
    ccle.sample.outlier.status.overlap.na.sum > 0,
    ];
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05.na[
    ccle.sample.outlier.status.overlap.na.sum > 0,
    ];


rnai.quantile.05 <- list();
outlier.gene.rnai.score.05 <- list();
nonoutlier.gene.rnai.score.05 <- list();
for (i in 1:nrow(ccle.sample.outlier.status.overlap.na)) {
    outlier.gene <- gene.rnai.breast.t.num.match.05.na[i, which(ccle.sample.outlier.status.overlap.na[i, ] == 1), drop = FALSE];
    non.outlier.gene <- gene.rnai.breast.t.num.match.05.na[i, -(which(ccle.sample.outlier.status.overlap.na[i, ] == 1))];

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
rownames(outlier.gene.rnai.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na);

nonoutlier.gene.rnai.score.05.mean <- sapply(nonoutlier.gene.rnai.score.05, function(x) {
    mean(na.omit(unlist(x)));
    });
nonoutlier.gene.rnai.score.05.mean <- data.frame(rnai = nonoutlier.gene.rnai.score.05.mean);
rownames(nonoutlier.gene.rnai.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na);

# Create difference matrix for RNAi scores
gene.rnai.diff.matrix.05 <- data.frame(
    out = outlier.gene.rnai.score.05.mean$rnai,
    non = nonoutlier.gene.rnai.score.05.mean$rnai,
    diff = outlier.gene.rnai.score.05.mean$rnai - nonoutlier.gene.rnai.score.05.mean$rnai,
    symbol = sub('\\..*', '', rownames(outlier.gene.rnai.score.05.mean))
    );
rownames(gene.rnai.diff.matrix.05) <- rownames(outlier.gene.rnai.score.05.mean);

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
    c('gene', 'dependency', 'diff', 'rnai', 'scatter'),
    width = 6,
    height = 5
    );
