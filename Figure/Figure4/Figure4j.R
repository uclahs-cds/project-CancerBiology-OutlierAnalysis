### HISTORY #####################################################################
# This script compares gene effect scores of outlier genes from both RNAi and
# Cas-CRISPR datasets, focusing on four specific genes: MECOM, FGFR2, FOXP4, WIPF2.
# Date: 2024-08-16

### DESCRIPTION #################################################################
# This script processes RNAi and Cas-CRISPR gene dependency data, specifically for four
# key outlier genes (MECOM, FGFR2, FOXP4, WIPF2). It calculates the gene effect scores
# for outliers and visualizes the comparison between the RNAi and Cas-CRISPR datasets
# using a box plot, where RNAi and Cas-CRISPR scores are grouped for each gene.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'p.me',
    'ccle.sample.outlier.status.overlap.na',
    'rnai.05.box',
    'cas.effect.breast.05.na',
    'gene.effect.diff.matrix.05.overlap'
    ));


gene.effect.diff.matrix.05.overlap.minus.05 <- gene.effect.diff.matrix.05.overlap[gene.effect.diff.matrix.05.overlap$diff < -0.5, ];
gene.effect.diff.matrix.05.overlap.minus.05.order <- gene.effect.diff.matrix.05.overlap.minus.05[order(gene.effect.diff.matrix.05.overlap.minus.05$diff), ];


effect.score.05.overlap.minus.05 <- cas.effect.breast.05.na[rownames(gene.effect.diff.matrix.05.overlap.minus.05.order), ];
outlier.status.05.overlap.minus.05 <- ccle.sample.outlier.status.overlap.na[rownames(gene.effect.diff.matrix.05.overlap.minus.05.order), ];

gene.name.box.05 <- NULL;
for (i in 1:nrow(outlier.status.05.overlap.minus.05)) {
    gene.name <- rep(sub('\\..*', '', rownames(effect.score.05.overlap.minus.05))[i], ncol(effect.score.05.overlap.minus.05));
    gene.name.box.05 <- c(gene.name.box.05, gene.name);
    }



effect.05.box <- data.frame(cbind(
    score = as.numeric(unlist(t(effect.score.05.overlap.minus.05))),
    gene = gene.name.box.05,
    status = as.numeric(unlist(t(outlier.status.05.overlap.minus.05)))
    ));
effect.05.box$score <- as.numeric(effect.05.box$score);
effect.05.box$status <- as.numeric(effect.05.box$status);


# 
# interesting.points <- cor.diff.cas.rnai.effect$cas < -0.5 & cor.diff.cas.rnai.effect$rnai < -0.4;
# gene.five.cas.rnai <-  na.omit(cor.diff.cas.rnai.effect$symbol[interesting.points]);
gene.five.cas.rnai <-  c("CTDSP2", "FGFR1", "FOXP4", "WIPF2", "MECOM","TNFSF10", "CCNE1");
# gene.five.cas.rnai <- unique(effect.05.box$gene)
rnai.05.box.4 <- rnai.05.box[rnai.05.box$gene %in% gene.five.cas.rnai, ];
rnai.05.box.4$label <- rep('RNAi', nrow(rnai.05.box.4));
effect.05.box.4 <- effect.05.box[effect.05.box$gene %in% gene.five.cas.rnai, ];
effect.05.box.4$label <- rep('Cas', nrow(effect.05.box.4));
rnai.effect.05.box.4 <- rbind(
    rnai.05.box.4,
    effect.05.box.4
    );
rnai.effect.05.box.4$order <- paste(rnai.effect.05.box.4$gene, rnai.effect.05.box.4$label, sep = '');
rnai.effect.05.box.4.order <- rnai.effect.05.box.4[order(rnai.effect.05.box.4$order), ];



dot.colours <- vector(length = nrow(rnai.effect.05.box.4.order));
dot.colours <- rep('grey70', nrow(rnai.effect.05.box.4.order));
dot.colours[rnai.effect.05.box.4.order$status == 1] <- 'dodgerblue2';
# dot.colours[rnai.05.box.part$status == 1][13:14] <- 'red2';

gene.order.rnai <- unique(rnai.effect.05.box.4.order$gene);

key <- list(
    text = list(
        lab = gene.order.rnai[1],
        cex = 1
        ),
    # x = 0.05,
    y = 0.93,
    text = list(
        lab = gene.order.rnai[2],
        cex = 1
        ),
    text = list(
        lab = gene.order.rnai[3],
        cex = 1
        ),
    # x = 0.9,
    y = 0.93,
    text = list(
        lab = gene.order.rnai[4],
        cex = 1
        ),
    # x = 0.9,
    y = 0.93,
    text = list(
        lab = gene.order.rnai[5],
        cex = 1
        ),
    # x = 0.9,
    y = 0.93,
    text = list(
        lab = gene.order.rnai[6],
        cex = 1
        ),
    # x = 0.9,
    y = 0.93,
    text = list(
        lab = gene.order.rnai[7],
        cex = 1
        ),
    # x = 0.9,
    y = 0.93
    );

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure4j')));

cas.rnai.example.box <- BoutrosLab.plotting.general::create.boxplot(
    formula = as.numeric(score) ~ order,
    data = rnai.effect.05.box.4.order,
    main = expression('Gene effect score of outlier genes'),
    outlier = TRUE,
    add.stripplot = TRUE,
    add.rectangle = TRUE,
    xleft.rectangle = seq(2.5, 14.5, 4),
    xright.rectangle = seq(4.5, 16.5, 4),
    ybottom.rectangle = -3,
    ytop.rectangle = 5,
    # set rectangle colour
    col.rectangle = 'grey',
    # set rectangle alpha (transparency)
    alpha.rectangle = 0.25,
    main.cex = 1.5,
    xaxis.lab = rep(c('CRISPR', 'RNAi'), 9),
    # xaxis.lab = gene.five.cas.rnai,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Gene effect score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.rot = 90,
    ylimits = c(-1.8, 0.9),
    key = key,
    sample.order = 'none',
    points.pch = 16,
    points.cex = 1,
    points.col = dot.colours,
    lwd = 1.2,
    # col = c('gold2'),
    col = rep(c('red3', 'dodgerblue3'), 9),
    alpha = 0.25
    );
cas.rnai.example.box;


save.outlier.figure(
    cas.rnai.example.box,
    c('Figure4j', 'cas', 'rnai', 'example', 'box'),
    width = 6.7,
    height = 6.7
    );

save.session.profile(file.path('output', 'Figure4j.txt'));
