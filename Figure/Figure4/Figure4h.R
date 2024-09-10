### HISTORY #####################################################################
# This script analyzes RNAi gene dependency scores for outlier genes identified 
# in CCLE. 
# Date: 2024-08-16


# Score of the overlap outliers
gene.rnai.diff.matrix.05.overlap.minus.05 <- gene.rnai.diff.matrix.05.overlap[gene.rnai.diff.matrix.05.overlap$diff < -0.4,];
gene.rnai.diff.matrix.05.overlap.minus.05 <- na.omit(gene.rnai.diff.matrix.05.overlap.minus.05);
gene.rnai.diff.matrix.05.overlap.minus.05.order <- gene.rnai.diff.matrix.05.overlap.minus.05[order(gene.rnai.diff.matrix.05.overlap.minus.05$diff),];


rnai.score.05.overlap.minus.05 <- gene.rnai.breast.t.num.match.05.na[rownames(gene.rnai.diff.matrix.05.overlap.minus.05.order),];
outlier.status.05.overlap.minus.05 <- ccle.sample.outlier.status.overlap.na[rownames(gene.rnai.diff.matrix.05.overlap.minus.05.order), colnames(rnai.score.05.overlap.minus.05)];

gene.name.box.05 <-NULL;
for(i in 1:nrow(outlier.status.05.overlap.minus.05)){
    gene.name <- rep(sub("\\..*","", rownames(rnai.score.05.overlap.minus.05))[i], ncol(rnai.score.05.overlap.minus.05));
    gene.name.box.05 <- c(gene.name.box.05, gene.name);
    }

rnai.05.box <- data.frame(cbind(
    score = as.numeric(unlist(t(rnai.score.05.overlap.minus.05))),
    gene = gene.name.box.05,
    status = as.numeric(unlist(t(outlier.status.05.overlap.minus.05)))));
rnai.05.box$score <- as.numeric(rnai.05.box$score);
rnai.05.box$status <- as.numeric(rnai.05.box$status);

dot.colours <- vector(length = nrow(rnai.05.box));
dot.colours <- rep('grey70', nrow(rnai.05.box));
dot.colours[rnai.05.box$status == 1]<- 'dodgerblue2';

rnai.05.box.plot <- BoutrosLab.plotting.general::create.boxplot(
    formula = score ~ gene,
    data = rnai.05.box,
    main =expression('Gene effect score of outlier genes'),
    outlier =TRUE,
    add.stripplot =TRUE,
    add.rectangle =TRUE,
    xleft.rectangle = seq(0.5,14.5,2),
    xright.rectangle = seq(1.5,15.5,2),
    ybottom.rectangle =-3,
    ytop.rectangle =5,
    col.rectangle ="grey",
    alpha.rectangle =0.25,
    main.cex =1.5,
    xlab.label =NULL,
    xlab.cex =0,
    ylab.label =expression('Gene effect score'),
    ylab.cex =1.3,
    yaxis.cex =1.1,
    xaxis.cex =1.1,
    xaxis.fontface =1,
    yaxis.fontface =1,
    yaxis.tck =c(0.2,0),
    xaxis.tck =c(0.2,0),
    xaxis.rot =90,
    sample.order ="decreasing",
    points.pch =16,
    points.cex =1,
    points.col = dot.colours,
    lwd =1.2,
    col =c('gold2'),
    alpha =0.25
    );




# Save the box plot as a PDF
pdf(
    file = generate.filename(
        'gene_effect_example_rnai', 
        'box', 
        'pdf'
        ), 
    width = 6, 
    height = 6
    );
rnai.05.box.plot;
dev.off();

# Save the box plot as a PNG
png(
    file = generate.filename(
        'gene_effect_example_rnai', 
        'box', 
        'png'
        ), 
    width = 6, 
    height = 6,
    unit = 'in', 
    res = 1200
    );
rnai.05.box.plot;
dev.off();




