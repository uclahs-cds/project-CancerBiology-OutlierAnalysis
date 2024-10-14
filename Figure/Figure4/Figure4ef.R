### HISTORY #####################################################################
# This script analyzes gene dependency scores for outlier genes identified in
# CCLE across multiple cancer datasets. It compares the gene dependency scores
# between outlier and non-outlier samples.
# Date: 2024-08-16

### DESCRIPTION #################################################################
# This script processes data from the CCLE to analyze
# gene dependency scores of identified outlier genes. It compares the scores between
# outlier and non-outlier samples and visualizes the results using a scatter plot.

### PREAMBLE ####################################################################
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

gene.dependency.breast.t <- t(gene.dependency.breast);
gene.dependency.breast.t.num.match <- as.data.frame(apply(gene.dependency.breast.t, 2, as.numeric));
rownames(gene.dependency.breast.t.num.match) <- rownames(gene.dependency.breast.t);
colnames(gene.dependency.breast.t.num.match) <- colnames(gene.dependency.breast.t);
gene.dependency.breast.t.num.match <- gene.dependency.breast.t.num.match[
    , colnames(fpkm.tumor.symbol.filter.ccle)
    ];

gene.dependency.breast.t.num.match.05 <- gene.dependency.breast.t.num.match[rownames(ccle.outlier.rank.fdr.05), ];
gene.dependency.breast.t.num.match.05.na <- na.omit(gene.dependency.breast.t.num.match.05);
ccle.sample.outlier.status.overlap <- ccle.sample.outlier.status[rownames(ccle.outlier.rank.fdr.05), ];
ccle.sample.outlier.status.overlap.na <- ccle.sample.outlier.status.overlap[rownames(gene.dependency.breast.t.num.match.05.na), ];



dependency.quantile.05 <- list();
outlier.gene.dependency.score.05 <- list();
nonoutlier.gene.dependency.score.05 <- list();
for (i in 1:nrow(ccle.sample.outlier.status.overlap.na)) {
    outlier.gene <- gene.dependency.breast.t.num.match.05.na[i, which(ccle.sample.outlier.status.overlap.na[i, ] == 1), drop = FALSE];
    non.outlier.gene <- gene.dependency.breast.t.num.match.05.na[i, -(which(ccle.sample.outlier.status.overlap.na[i, ] == 1))];
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


gene.dependency.diff.matrix.05 <- data.frame(
    out = outlier.gene.dependency.score.05.mean$dependency,
    non = nonoutlier.gene.dependency.score.05.mean$dependency,
    diff = outlier.gene.dependency.score.05.mean$dependency - nonoutlier.gene.dependency.score.05.mean$dependency,
    symbol = sub('\\..*', '', rownames(outlier.gene.dependency.score.05.mean))
    );

load.multiple.computed.variables(c(
    'outlier.symbol'
    ));

gene.dependency.diff.matrix.05.overlap <- gene.dependency.diff.matrix.05[
    gene.dependency.diff.matrix.05$symbol %in% outlier.symbol$unique,
    ];

cache.multiple.computed.variables(c(
    'ccle.sample.outlier.status.overlap'
    ));

# Filter overlapping genes

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
    c('Figure4ef', 'gene', 'dependency', 'diff', 'scatter'),
    width = 6,
    height = 5
    );

# Filter for genes with dependency score difference > 0.5
gene.dependency.diff.matrix.05.overlap.minus.05 <- gene.dependency.diff.matrix.05.overlap[
    gene.dependency.diff.matrix.05.overlap$diff > 0.5,
    ];
gene.dependency.diff.matrix.05.overlap.minus.05.order <- gene.dependency.diff.matrix.05.overlap.minus.05[
    order(gene.dependency.diff.matrix.05.overlap.minus.05$diff),
    ];

# Prepare data for box plot
dependency.score.05.overlap.minus.05 <- gene.dependency.breast.t.num.match.05.na[
    as.numeric(rownames(gene.dependency.diff.matrix.05.overlap.minus.05.order)),
    ];
outlier.status.05.overlap.minus.05 <- ccle.sample.outlier.status.overlap.na[
    as.numeric(rownames(gene.dependency.diff.matrix.05.overlap.minus.05.order)),
    ];

# Create gene name box
gene.name.box.05 <- NULL;
for (i in 1:nrow(outlier.status.05.overlap.minus.05)) {
    gene.name <- rep(sub('\\..*', '', rownames(dependency.score.05.overlap.minus.05))[i], ncol(dependency.score.05.overlap.minus.05));
    gene.name.box.05 <- c(gene.name.box.05, gene.name);
    }

# Combine dependency score and outlier status into a data frame
dependency.05.box <- data.frame(
    cbind(
        score = as.numeric(unlist(t(dependency.score.05.overlap.minus.05))),
        gene = gene.name.box.05,
        status = as.numeric(unlist(t(outlier.status.05.overlap.minus.05)))
        )
    );
dependency.05.box$score <- as.numeric(dependency.05.box$score);
dependency.05.box$status <- as.numeric(dependency.05.box$status);

# Filter out specific genes from the data
dependency.05.box.part <- dependency.05.box[!(
    dependency.05.box$gene %in% c('CASC3', 'TSEN54', 'MRPL21')
    ), ];

# Add specific genes to the box plot data
gene.dependency.diff.matrix.05.overlap.plus.05 <- gene.dependency.diff.matrix.05.overlap[
    gene.dependency.diff.matrix.05.overlap$symbol %in% c('TACC3', 'CCT2'),
    ];
dependency.score.05.overlap.plus.05 <- gene.dependency.breast.t.num.match.05.na[
    rownames(gene.dependency.diff.matrix.05.overlap.plus.05),
    ];
outlier.status.05.overlap.plus.05 <- ccle.sample.outlier.status.overlap.na[
    rownames(gene.dependency.diff.matrix.05.overlap.plus.05),
    ];
dependency.05.box.plus <- data.frame(
    cbind(
        score = as.numeric(unlist(t(dependency.score.05.overlap.plus.05))),
        gene = c(
            rep(gene.dependency.diff.matrix.05.overlap.plus.05$symbol[1], ncol(dependency.score.05.overlap.minus.05)),
            rep(gene.dependency.diff.matrix.05.overlap.plus.05$symbol[2], ncol(dependency.score.05.overlap.minus.05))
            ),
        status = as.numeric(unlist(t(outlier.status.05.overlap.plus.05)))
        )
    );

# Combine the additional data into the main data frame
dependency.05.box.part <- rbind(dependency.05.box.part, dependency.05.box.plus);

# Convert score and status to numeric
dependency.05.box.part$score <- as.numeric(dependency.05.box.part$score);
dependency.05.box.part$status <- as.numeric(dependency.05.box.part$status);

# Further filter out specific genes from the data
dependency.05.box.part.4 <- dependency.05.box.part[!(
    dependency.05.box.part$gene %in% c('CCT2', 'TACC3', 'MSL1', 'RTN4IP1')
    ), ];

# Set colors for the box plot
dot.colours <- rep('grey70', nrow(dependency.05.box.part.4));
dot.colours[dependency.05.box.part.4$status == 1] <- 'red2';

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure4f')));

# Create the box plot
dependency.05.box.plot <- BoutrosLab.plotting.general::create.boxplot(
    formula = score ~ gene,
    data = dependency.05.box.part.4,
    main = expression('Gene dependency score of outlier genes'),
    outlier = TRUE,
    add.stripplot = TRUE,
    add.rectangle = TRUE,
    xleft.rectangle = seq(0.5, 14.5, 2),
    xright.rectangle = seq(1.5, 15.5, 2),
    ybottom.rectangle = -3,
    ytop.rectangle = 5,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    main.cex = 1.5,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Gene dependency score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.rot = 90,
    sample.order = 'increasing',
    points.pch = 16,
    points.cex = 1,
    points.col = dot.colours,
    lwd = 1.2,
    col = c('gold2'),
    alpha = 0.25
    );


save.outlier.figure(
    dependency.05.box.plot,
    c('Figure4ef', 'gene', 'dependency', 'example', 'box'),
    width = 6,
    height = 5
    );

save.session.profile(file.path('output', 'Figure4ef.txt'));
