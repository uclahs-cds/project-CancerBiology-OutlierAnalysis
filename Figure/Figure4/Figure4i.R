### HISTORY #####################################################################  
# This script analyzes the correlation between Cas-CRISPR and RNAi gene  
# dependency scores.  
# Date: 2024-08-16  

### DESCRIPTION #################################################################  
# This script processes Cas-CRISPR and RNAi gene dependency scores, calculating  
# the differences in gene effect scores for outliers and non-outliers. The data  
# is filtered and matched across both datasets. The correlation between Cas-CRISPR  
# and RNAi gene dependency score differences is then visualized in a scatter plot,  
# highlighting genes of interest based on specific thresholds.  

### PREAMBLE ####################################################################  
# Load necessary libraries  
library(BoutrosLab.plotting.general);  
library(BoutrosLab.utilities);  

# Source the helper library
args <- commandArgs();
source(file.path(
    dirname(dirname(normalizePath(sub('^--file=', '', args[grep('^--file=', args)])))),
    'common_functions.R'
    ));
# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-09-10_Figure4.rda'));


# Get gene effect score from Cas-CRISPR dataset

effect.quantile.05 <- list();
outlier.gene.effect.score.05 <- list();
nonoutlier.gene.effect.score.05 <- list();
for (i in 1:nrow(ccle.sample.outlier.status.overlap.na)) {
    outlier.gene <- cas.effect.breast.05.na[i, which(ccle.sample.outlier.status.overlap.na[i, ] == 1), drop = FALSE];
    non.outlier.gene <- cas.effect.breast.05.na[i, -(which(ccle.sample.outlier.status.overlap.na[i, ] == 1))];
    ecdf.obj <- ecdf(as.numeric(non.outlier.gene));
    quantile.value <- ecdf.obj(outlier.gene);
    ecdf.obj <- ecdf(as.numeric(non.outlier.gene));
    quantile.value <- ecdf.obj(outlier.gene);
    quantile.value <- t(data.frame(quantile.value));
    colnames(quantile.value) <- colnames(outlier.gene);
    rownames(quantile.value) <- rownames(outlier.gene);
    effect.quantile.05[[i]] <- quantile.value;

    outlier.gene.effect.score.05[[i]] <- outlier.gene;
    nonoutlier.gene.effect.score.05[[i]] <- non.outlier.gene;
    }



#   - Use absolute difference of gene effect
outlier.gene.effect.score.05.mean <- sapply(outlier.gene.effect.score.05, function(x) {
    mean(unlist(x)) # Convert inner list to a numeric vector and calculate mean
    });
outlier.gene.effect.score.05.mean <- data.frame(effect = outlier.gene.effect.score.05.mean);
rownames(outlier.gene.effect.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na);

nonoutlier.gene.effect.score.05.mean <- sapply(nonoutlier.gene.effect.score.05, function(x) {
    mean(unlist(x)) # Convert inner list to a numeric vector and calculate mean
    });
nonoutlier.gene.effect.score.05.mean <- data.frame(effect = nonoutlier.gene.effect.score.05.mean);
rownames(nonoutlier.gene.effect.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na);

gene.effect.diff.matrix.05 <- data.frame(
    out = outlier.gene.effect.score.05.mean$effect,
    non = nonoutlier.gene.effect.score.05.mean$effect,
    diff = outlier.gene.effect.score.05.mean$effect - nonoutlier.gene.effect.score.05.mean$effect,
    symbol = sub('\\..*', '', rownames(outlier.gene.effect.score.05.mean))
    );
rownames(gene.effect.diff.matrix.05) <- rownames(outlier.gene.effect.score.05.mean);

# Prepare the datasets by matching and removing NA values
gene.rnai.diff.matrix.05.overlap.na <- na.omit(gene.rnai.diff.matrix.05.overlap);
gene.effect.diff.matrix.05.overlap.match.rnai <- gene.effect.diff.matrix.05.overlap[rownames(gene.rnai.diff.matrix.05.overlap.na), ];
gene.effect.diff.matrix.05.overlap.match.rnai.na <- na.omit(gene.effect.diff.matrix.05.overlap.match.rnai);
gene.rnai.diff.matrix.05.overlap.match.cas <- gene.rnai.diff.matrix.05.overlap[rownames(gene.effect.diff.matrix.05.overlap.match.rnai.na), ];

# Create a data frame with RNAi and Cas-CRISPR differences and gene symbols
cor.diff.cas.rnai.effect <- data.frame(cbind(
    rnai = gene.rnai.diff.matrix.05.overlap.match.cas$diff,
    cas = gene.effect.diff.matrix.05.overlap.match.rnai.na$diff,
    symbol = gene.effect.diff.matrix.05.overlap.match.rnai.na$symbol
    ));
cor.diff.cas.rnai.effect$rnai <- as.numeric(cor.diff.cas.rnai.effect$rnai);
cor.diff.cas.rnai.effect$cas <- as.numeric(cor.diff.cas.rnai.effect$cas);

# Set dot colors
dot.colours <- vector(length = nrow(cor.diff.cas.rnai.effect));
dot.colours <- rep('grey60', nrow(cor.diff.cas.rnai.effect));
dot.colours[cor.diff.cas.rnai.effect$cas < -0.5 & cor.diff.cas.rnai.effect$rnai < -0.4] <- 'dodgerblue2';


interesting.points <- cor.diff.cas.rnai.effect$cas < -0.5 & cor.diff.cas.rnai.effect$rnai < -0.4;
text.x <- na.omit(cor.diff.cas.rnai.effect$cas[interesting.points]);
text.y <- na.omit(cor.diff.cas.rnai.effect$rnai[interesting.points]);
text.labels <- na.omit(cor.diff.cas.rnai.effect$symbol[interesting.points]);

# Create the scatter plot
cor.scatter.effect.cas.rnai <- BoutrosLab.plotting.general::create.scatterplot(
    formula = rnai ~ cas,
    data = cor.diff.cas.rnai.effect,
    main = expression('Correlation between Cas-CRISPR and RNAi'),
    ylab.label = expression('Difference of gene effect score (RNAi)'),
    xlab.label = expression('Difference of gene effect score (Cas-CRISPR)'),
    col = dot.colours,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    alpha = 0.8,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.5,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    ylimits = c(-1.2, 1.05),
    xlimits = c(-2.1, 0.9),
    type = c('p', 'r', 'g'),
    lwd = 1.5,
    add.grid = TRUE,
    grid.colour = 'grey80',
    cex = 1,
    add.text = TRUE,
    text.x = text.x,
    text.y = text.y,
    text.labels = text.labels,
    text.fontface = 1,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = cor.diff.cas.rnai.effect$rnai,
                    y = cor.diff.cas.rnai.effect$cas,
                    label.items = c('spearman', 'spearman.p'),
                    alpha.background = 0,
                    key.cex = 1.1
                    )
                ),
            x = 0.03,
            y = 0.95,
            corner = c(0, 1)
            )
        ),
    resolution = 300
    );


save.outlier.figure(
    cor.scatter.effect.cas.rnai,
    c('Figure4i', 'cas', 'rnai', 'cor', 'scatter'),
    width = 5,
    height = 5
    );

save.session.profile(file.path('output', 'Figure4i.txt'));
