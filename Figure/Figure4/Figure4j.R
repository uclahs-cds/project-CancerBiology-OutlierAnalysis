### HISTORY #####################################################################
# This script compares gene effect scores of outlier genes from both RNAi and 
# Cas-CRISPR datasets, focusing on four specific genes: MECOM, FGFR2, FOXP4, WIPF2.
# Date: 2024-08-16



# Get RNAi data
# score of the overlap outliers
gene.rnai.diff.matrix.05.overlap.minus.05 <- gene.rnai.diff.matrix.05.overlap[gene.rnai.diff.matrix.05.overlap$diff < -0.4,];
gene.rnai.diff.matrix.05.overlap.minus.05 <- na.omit(gene.rnai.diff.matrix.05.overlap.minus.05);
gene.rnai.diff.matrix.05.overlap.minus.05.order <- gene.rnai.diff.matrix.05.overlap.minus.05[order(gene.rnai.diff.matrix.05.overlap.minus.05$diff),];

rnai.score.05.overlap.minus.05 <- gene.rnai.breast.t.num.match.05.na[rownames(gene.rnai.diff.matrix.05.overlap.minus.05.order),];
outlier.status.05.overlap.minus.05 <- ccle.sample.outlier.status.overlap.na[rownames(gene.rnai.diff.matrix.05.overlap.minus.05.order),];

gene.name.box.05 <- NULL;
for (i in 1:nrow(outlier.status.05.overlap.minus.05)) {
    gene.name <- rep(sub("\\..*", "", rownames(rnai.score.05.overlap.minus.05))[i], ncol(rnai.score.05.overlap.minus.05));
    gene.name.box.05 <- c(gene.name.box.05, gene.name);
    }

rnai.05.box <- data.frame(cbind(
    score = as.numeric(unlist(t(rnai.score.05.overlap.minus.05))),
    gene = gene.name.box.05,
    status = as.numeric(unlist(t(outlier.status.05.overlap.minus.05)))
    ));
rnai.05.box$score <- as.numeric(rnai.05.box$score);
rnai.05.box$status <- as.numeric(rnai.05.box$status);



# Get Cas-CRISPR data
gene.effect.diff.matrix.05.overlap.minus.05 <- gene.effect.diff.matrix.05.overlap[gene.effect.diff.matrix.05.overlap$diff < -0.5,];
gene.effect.diff.matrix.05.overlap.minus.05.order <- gene.effect.diff.matrix.05.overlap.minus.05[order(gene.effect.diff.matrix.05.overlap.minus.05$diff),];

effect.score.05.overlap.minus.05 <- cas.effect.breast.05.na[rownames(gene.effect.diff.matrix.05.overlap.minus.05.order),];
outlier.status.05.overlap.minus.05 <- sample.outlier.05.overlap.na[rownames(gene.effect.diff.matrix.05.overlap.minus.05.order),];

gene.name.box.05 <- NULL;
for (i in 1:nrow(outlier.status.05.overlap.minus.05)) {
    gene.name <- rep(sub("\\..*", "", rownames(effect.score.05.overlap.minus.05))[i], ncol(effect.score.05.overlap.minus.05));
    gene.name.box.05 <- c(gene.name.box.05, gene.name);
    }

effect.05.box <- data.frame(cbind(
    score = as.numeric(unlist(t(effect.score.05.overlap.minus.05))),
    gene = gene.name.box.05,
    status = as.numeric(unlist(t(outlier.status.05.overlap.minus.05)))
    ));
effect.05.box$score <- as.numeric(effect.05.box$score);
effect.05.box$status <- as.numeric(effect.05.box$status);




# Select four specific genes from both datasets
gene.five.cas.rnai <-c('MECOM','FGFR2','FOXP4','WIPF2');
rnai.05.box.4 <- rnai.05.box[rnai.05.box$gene %in% gene.five.cas.rnai,];
rnai.05.box.4$label <-rep('RNAi', nrow(rnai.05.box.4));

effect.05.box.4 <- effect.05.box[effect.05.box$gene %in% gene.five.cas.rnai,];
effect.05.box.4$label <-rep('Cas', nrow(effect.05.box.4));

# Combine the data from RNAi and Cas-CRISPR datasets
rnai.effect.05.box.4 <- rbind(
    rnai.05.box.4,
    effect.05.box.4
    );
rnai.effect.05.box.4$order <- paste(rnai.effect.05.box.4$gene, rnai.effect.05.box.4$label, sep ='');
rnai.effect.05.box.4.order <- rnai.effect.05.box.4[order(rnai.effect.05.box.4$order),];

# Set dot colors
dot.colours <- vector(length=nrow(rnai.effect.05.box.4.order));
dot.colours <-rep('grey70', nrow(rnai.effect.05.box.4.order));
dot.colours[rnai.effect.05.box.4.order$status ==1]<-'dodgerblue2';


key <-list(
    text =list(
        lab ='MECOM', 
        cex =1
        ), 
    x =0.05,
    y =0.93,
    text =list(
        lab ='FGFR2', 
        cex =1
        ), 
    text =list(
        lab ='FOXP4', 
        cex =1
        ), 
    x =0.9,
    y =0.93,
    text =list(
        lab ='WIPF2', 
        cex =1
        ), 
    x =0.9,
    y =0.93
    );

# Create the boxplot comparing RNAi and Cas-CRISPR
cas.rnai.example.box <- BoutrosLab.plotting.general::create.boxplot(
    formula = score ~ order,
    data = rnai.effect.05.box.4.order,
    main =expression('Gene effect score of outlier genes'),
    outlier =TRUE,
    add.stripplot =TRUE,
    add.rectangle =TRUE,
    xleft.rectangle = seq(2.5,10.5,4),
    xright.rectangle = seq(4.5,12.5,4),
    ybottom.rectangle =-3,
    ytop.rectangle =5,
    col.rectangle ="grey",
    alpha.rectangle =0.25,
    main.cex =1.5,
    xaxis.lab =rep(c('Cas-CRISPR','RNAi'),4),
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
    ylimits =c(-1.8,0.9),
    key = key,
    sample.order ="none",
    points.pch =16,
    points.cex =1,
    points.col = dot.colours,
    lwd =1.2,
    col =rep(c('red3','dodgerblue3'),4),
    alpha =0.25
    );


# Save the plot as a PNG
output.directory <- get0('output.directory', ifnotfound = 'figures');

# Save the box plot as a PDF
pdf(
    file = file.path(output.directory, 'Figure_4_j.pdf'),
    width = 5, 
    height = 6.5
    );
cas.rnai.example.box;
dev.off();

# Save the box plot as a PNG
png(
    file = file.path(output.directory, 'Figure_4_j.png'),
    width = 5, 
    height = 6.5,
    unit = 'in', 
    res = 1200
    );
cas.rnai.example.box;
dev.off();




