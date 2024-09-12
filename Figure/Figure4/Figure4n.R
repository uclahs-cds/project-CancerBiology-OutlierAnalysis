### HISTORY #####################################################################
# This script generates a boxplot comparing gene effect scores between CRISPR 
# and RNAi for the gene SYK.
# Date: 2024-08-16

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));

i <-'SYK';

# Extract RNAi and CRISPR data
rnai.score.i <- gene.rnai.breast.t.num.match.05.na[sub("\\..*","", rownames(gene.rnai.breast.t.num.match.05.na)) %in% i,];
rnai.status.i <- ccle.sample.outlier.status.overlap.na.rnai[sub("\\..*","", rownames(ccle.sample.outlier.status.overlap.na.rnai)) %in% i,];
effect.score.i <- cas.effect.breast.05.na[sub("\\..*","", rownames(cas.effect.breast.05.na)) %in% i,];
effect.status.i <- ccle.sample.outlier.status.overlap.na[sub("\\..*","", rownames(ccle.sample.outlier.status.overlap.na)) %in% i,];

# Combine RNAi and CRISPR data into a single data frame
rnai.cas.effect.score.drug.gene.1 <- data.frame(
    score = c(
        unlist(rnai.score.i), 
        unlist(effect.score.i)
        ),
    status = c(
        unlist(rnai.status.i), 
        unlist(effect.status.i)
        ),
    group = c(
        rep('b',length(rnai.score.i)),
        rep('a',length(effect.score.i))
        )
    );

# Set colors for the points based on outlier status
dot.colours <- vector(length= nrow(rnai.cas.effect.score.drug.gene.1));
dot.colours <- rep('grey70', nrow(rnai.cas.effect.score.drug.gene.1));
dot.colours[rnai.cas.effect.score.drug.gene.1$status == 1] <- 'dodgerblue3';

# Create the boxplot
rnai.cas.effect.box.plot <- BoutrosLab.plotting.general::create.boxplot(
    formula = score ~ group,
    data = rnai.cas.effect.score.drug.gene.1,
    main = expression('Gene effect score of outlier genes'),
    outlier = TRUE,
    add.stripplot = TRUE,
    add.rectangle = TRUE,
    xleft.rectangle = 1.5,
    xright.rectangle = 3,
    ybottom.rectangle = -3,
    ytop.rectangle = 5,
    col.rectangle = "grey",
    alpha.rectangle = 0.25,
    main.cex = 1.5,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Gene effect score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.lab = c('CRISPR','RNAi'),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    xaxis.rot = 90,
    yat = seq(-2.5, 1.0, 0.5),
    sample.order = "none",
    points.pch = 16,
    points.cex = 1,
    points.col = dot.colours,
    lwd = 1.2,
    col = c('red3','dodgerblue3'),
    alpha = 0.25
    );

save.outlier.figure(
    rnai.cas.effect.box.plot,
    c(i, 'box', 'cas', 'rnai'),
    width = 3.8,
    height = 6
    );
