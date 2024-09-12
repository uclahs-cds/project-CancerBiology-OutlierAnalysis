### HISTORY #####################################################################
# This script analyzes gene dependency scores for outlier genes identified in 
# CCLE across multiple cancer datasets. It compares the gene dependency scores 
# between outlier and non-outlier samples.
# Date: 2024-08-16

library(BoutrosLab.plotting.general);

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));

gene.dependency.breast.t <- t(gene.dependency.breast);
gene.dependency.breast.t.num.match <- as.data.frame(apply(gene.dependency.breast.t, 2, as.numeric));
rownames(gene.dependency.breast.t.num.match) <- rownames(gene.dependency.breast.t);
colnames(gene.dependency.breast.t.num.match) <- colnames(gene.dependency.breast.t);
gene.dependency.breast.t.num.match <- gene.dependency.breast.t.num.match[
    , colnames(fpkm.tumor.symbol.filter.ccle)
    ];

# Filter for FDR < 0.05
ccle.sample.outlier.status.overlap <- ccle.sample.outlier.status[
    rownames(ccle.outlier.rank.fdr.05), 
    ];
gene.dependency.breast.t.num.match.05 <- gene.dependency.breast.t.num.match[
    rownames(ccle.outlier.rank.fdr.05), 
    ];
gene.dependency.breast.t.num.match.05.na <- na.omit(gene.dependency.breast.t.num.match.05);
ccle.sample.outlier.status.overlap.na <- ccle.sample.outlier.status.overlap[
    rownames(gene.dependency.breast.t.num.match.05.na), 
    ];

ccle.sample.outlier.status.na <- ccle.sample.outlier.status[rownames(gene.dependency.breast.t.num.match.05.na),];


dependency.quantile.05 <- list();
outlier.gene.dependency.score.05 <- list();
nonoutlier.gene.dependency.score.05 <- list();
for (i in 1:nrow(ccle.sample.outlier.status.na)) {
    outlier.gene <- gene.dependency.breast.t.num.match.05.na[i, which(ccle.sample.outlier.status.na[i,] == 1), drop = FALSE];
    non.outlier.gene <- gene.dependency.breast.t.num.match.05.na[i, -(which(ccle.sample.outlier.status.na[i,] == 1))];
    ecdf.obj <- ecdf(as.numeric(non.outlier.gene));
    quantile.value <- ecdf.obj(outlier.gene);
    ecdf.obj <- ecdf(as.numeric(non.outlier.gene));
    quantile.value <- ecdf.obj(outlier.gene);
    quantile.value <- t(data.frame(quantile.value));
    colnames(quantile.value) <- colnames(outlier.gene);
    rownames(quantile.value) <- rownames(outlier.gene);
    dependency.quantile.05[[i]] <- quantile.value;
    
    outlier.gene.dependency.score.05[[i]] <- outlier.gene;
    nonoutlier.gene.dependency.score.05[[i]] <- non.outlier.gene;
    }


# Calculate gene dependency score differences
calculte.mean <- function(scores) {
    sapply(scores, function(x) mean(unlist(x)));
    }

outlier.gene.dependency.score.05.mean <- data.frame(
    dependency = calculte.mean(outlier.gene.dependency.score.05)
    );
nonoutlier.gene.dependency.score.05.mean <- data.frame(
    dependency = calculte.mean(nonoutlier.gene.dependency.score.05)
    );

rownames(outlier.gene.dependency.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na);
rownames(nonoutlier.gene.dependency.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na);

p.dependency.05 <- wilcox.test(
    outlier.gene.dependency.score.05.mean$dependency, 
    nonoutlier.gene.dependency.score.05.mean$dependency, 
    alternative = "two.sided", conf.int = TRUE
    );

gene.dependency.diff.matrix.05 <- data.frame(
    out = outlier.gene.dependency.score.05.mean$dependency,
    non = nonoutlier.gene.dependency.score.05.mean$dependency,
    diff = outlier.gene.dependency.score.05.mean$dependency - nonoutlier.gene.dependency.score.05.mean$dependency,
    symbol = sub("\\..*", "", rownames(outlier.gene.dependency.score.05.mean))
    );

# Filter overlapping genes
gene.dependency.diff.matrix.05.overlap <- gene.dependency.diff.matrix.05[
    gene.dependency.diff.matrix.05$symbol %in% five.data.outlier.symbol, 
    ];

# Set colors
dot.colours <- ifelse(
    gene.dependency.diff.matrix.05.overlap$diff < -0.5, 'dodgerblue3',
    ifelse(gene.dependency.diff.matrix.05.overlap$diff > 0.5, 'red2', 'grey30')
    );

# Label interesting points
interesting.points <- gene.dependency.diff.matrix.05.overlap$diff > 0.75;
text.x <- na.omit(gene.dependency.diff.matrix.05.overlap$non[interesting.points]);
text.y <- na.omit(gene.dependency.diff.matrix.05.overlap$out[interesting.points]);
text.labels <- na.omit(gene.dependency.diff.matrix.05.overlap$symbol[interesting.points]);

# Create scatter plot
gene.scatter.05.minus.overlap.label <- create.scatterplot(
    formula = diff ~ non,
    data = gene.dependency.diff.matrix.05.overlap,
    col = dot.colours,
    alpha = 0.6,
    ylimits = c(-1.05, 1.2),
    xlimits = c(-0.07, 1.07),
    yat = seq(-1.0, 1.0, 0.5),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = FALSE,
    grid.colour = 'grey80',
    cex = 1,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1,
    main.cex = 1.5,
    left.padding = 0,
    main = expression('Gene dependency score of overlap outliers'),
    xlab.label = expression(paste('Mean of gene dependency score of non-outlier samples')),
    ylab.label = expression('Difference of gene dependency score'),
    add.text = TRUE,
    text.x = text.x,
    text.y = text.y,
    text.labels = text.labels,
    text.fontface = 1,
    abline.h = c(0, 0.5),
    abline.col = c('grey30', 'red2'),
    abline.lwd = c(1.5, 2),
    abline.lty = c(1, 3)
    );


save.outlier.figure(
    gene.scatter.05.minus.overlap.label,
    c('gene', 'dependency', 'diff', 'scatter'),
    width = 6,
    height = 5
    );
