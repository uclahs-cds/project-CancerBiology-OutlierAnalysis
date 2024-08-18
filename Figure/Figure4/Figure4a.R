### HISTORY #####################################################################
# This script visualizes the expression of outliers identified in CCLE across 
# multiple datasets, overlapped with TCGA-BRCA, METABRIC, I-SPY2, MATADOR, and 
# ICGC BRCA-EU.
# Date: 2024-08-16



# sample outlier status
sample.outlier.05.only <- sample.outlier.05[rownames(gene.rank.p.value.one.gene.p0.fdr.05),];
# overlapped with tissue outliers
sample.outlier.05.only.five <- sample.outlier.05.only[sub("\\..*", "", rownames(sample.outlier.05.only)) %in% five.data.outlier.symbol,];

# Tissue overall outlier status
sample.outlier.05.only.five.symbol <- sub("\\..*", "", rownames(sample.outlier.05.only.five));

brca.row.name <- rownames(fpkm.tumor.symbol.filter.brca)[match(sample.outlier.05.only.five.symbol, fpkm.tumor.symbol.filter.brca$Symbol)]; 
outlier.patient.tag.01.brca.sum.overlap <- outlier.patient.tag.01.brca.sum[brca.row.name];

meta.row.name <- rownames(fpkm.tumor.symbol.filter.meta.symbol)[match(sample.outlier.05.only.five.symbol, fpkm.tumor.symbol.filter.meta.symbol$Symbol)]; 
outlier.patient.tag.01.meta.sum <- apply(outlier.patient.tag.01.meta, 1, sum);
outlier.patient.tag.01.meta.sum.overlap <- outlier.patient.tag.01.meta.sum[meta.row.name];

outlier.patient.tag.01.ispy.sum <- apply(outlier.patient.tag.01.ispy, 1, sum);
outlier.patient.tag.01.ispy.sum.overlap <- outlier.patient.tag.01.ispy.sum[sample.outlier.05.only.five.symbol];

matador.row.name <- rownames(fpkm.tumor.symbol.filter.metador.symbol)[match(sample.outlier.05.only.five.symbol, fpkm.tumor.symbol.filter.metador.symbol$Symbol)]; 
outlier.patient.tag.01.metador.sum <- apply(outlier.patient.tag.01.metador, 1, sum);
outlier.patient.tag.01.metador.sum.overlap <- outlier.patient.tag.01.metador.sum[matador.row.name];

outlier.patient.tag.01.icgc.sum.overlap <- icgc.outlier.symbol[match(sample.outlier.05.only.five.symbol, icgc.outlier.symbol)]; 


ccle.overlap.outlier.05.five.tissue <- data.frame(cbind(
    brca = outlier.patient.tag.01.brca.sum.overlap,
    meta = outlier.patient.tag.01.meta.sum.overlap,
    ispy = outlier.patient.tag.01.ispy.sum.overlap,
    matador = outlier.patient.tag.01.metador.sum.overlap,
    icgc = outlier.patient.tag.01.icgc.sum.overlap
    ));
rownames(ccle.overlap.outlier.05.five.tissue) <- sample.outlier.05.only.five.symbol;

ccle.overlap.outlier.05.five.tissue.0.1 <- apply(ccle.overlap.outlier.05.five.tissue, 1, function(x) { 
    x[is.na(x)] <- 0;
    ifelse(x > 0, 1, 0);});
ccle.overlap.outlier.05.five.tissue.0.1 <- t(ccle.overlap.outlier.05.five.tissue.0.1);

data.mean.outlier.05.only.five <- data.mean[rownames(sample.outlier.05.only.five),];

distance.matrix.t <- dist(t(data.mean.outlier.05.only.five), method = "euclidean");
fit.t <- hclust(distance.matrix.t, method = "ward.D2");
data.mean.outlier.05.only.five.t.order <- data.mean.outlier.05.only.five[,fit.t$order];

distance.matrix <- dist(data.mean.outlier.05.only.five.t.order, method = "euclidean");
fit <- hclust(distance.matrix, method = "ward.D2");
data.mean.outlier.05.only.five.t.p.order <- data.mean.outlier.05.only.five.t.order[fit$order,];
colnames(data.mean.outlier.05.only.five.t.p.order) <- gsub("\\.", "-", colnames(data.mean.outlier.05.only.five.t.p.order));

ccle.overlap.outlier.05.five.tissue.0.1.order.zscore <- ccle.overlap.outlier.05.five.tissue.0.1[sub("\\..*", "", rownames(data.mean.outlier.05.only.five.t.p.order)),];


library("RColorBrewer");
down.col <- brewer.pal(11,"RdBu")[11];
up.col <- brewer.pal(11,"RdBu")[1];

data.mean.outlier.05.only.five.t.p.order[data.mean.outlier.05.only.five.t.p.order > 5] <- 5;

# Row/col color - show the outlier status
sample.outlier.05.only.five.match <- sample.outlier.05.only.five[match(rownames(data.mean.outlier.05.only.five.t.p.order), rownames(sample.outlier.05.only.five)), match(colnames(data.mean.outlier.05.only.five.t.p.order), colnames(sample.outlier.05.only.five))];
sample.outlier.05.only.five.match.row.sum <- rowSums(sample.outlier.05.only.five.match);
sample.outlier.05.only.five.match.col.sum <- rowSums(sample.outlier.05.only.five.match);

main.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(data.mean.outlier.05.only.five.t.p.order),
    clustering.method = 'none',
    cluster.dimensions = 'none',
    plot.dendrograms = 'none',
    yaxis.cex = 0.2,
    xaxis.cex = 2,
    main = expression('Outlier status'),
    main.cex = 1.5,
    grid.col = TRUE,
    print.colour.key = FALSE,
    ylab.label = expression('The overlap outlier gene set'),
    ylab.cex = 1.3,
    xlab.label = expression('CCLE'),
    force.grid.row = TRUE,
    force.grid.col = TRUE,
    grid.colour = 'white',
    xlab.cex = 1.3,
    yaxis.tck = 0,
    xaxis.tck = 0,
    row.pos = which(t(sample.outlier.05.only.five.match) > 0, arr.ind = TRUE)[,2] + 1.5,
    col.pos = which(t(sample.outlier.05.only.five.match) > 0, arr.ind = TRUE)[,1],
    cell.text = rep(".", times = sum(t(sample.outlier.05.only.five.match) > 0)),
    text.col = "white",
    text.cex = 1.5,
    axes.lwd = 0.8,
    colour.scheme = c(down.col, 'white', up.col),
    colour.centering.value = 0,
    at = seq(-5, 5, 0.001),
    colourkey.cex = 1,
    covariate.legend = legend.sample,
    legend.side = 'right',
    legend.title.just = 'left',
    legend.cex = 1,
    resolution = 1000
    );  


ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore <- data.frame(ccle.overlap.outlier.05.five.tissue.0.1.order.zscore);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$brca <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$brca > 0, 2, 1);  
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$meta <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$meta > 0, 3, 1);  
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$ispy <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$ispy > 0, 4, 1);  
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$matador <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$matador > 0, 5, 1);  
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$icgc <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$icgc > 0, 6, 1);  

all.col <- c('grey95', brca.col, meta.col, ispy.col, matador.col, ccle.col);

sub.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore),
    clustering.method = 'none',
    colour.scheme = all.col, 
    total.colours = 7,
    row.colour = 'black',
    col.colour = 'black',
    grid.row = TRUE, 
    grid.col = TRUE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    xaxis.lab = c('TCGA-BRCA', 'METABRIC', 'I-SPY2', 'MATADOR', 'ICGC BRCA-EU'),
    xaxis.fontface = 1,
    xaxis.rot = 90,
    xaxis.cex = 1,
    print.colour.key = FALSE
    );


legend.sample.grob <- BoutrosLab.plotting.general:::legend.grob(
    list(
        legend = list(
            title = expression(underline('z-score')),
            continuous = TRUE,
            colours = c(down.col,'white', up.col),
            total.colours = 100,
            labels = c('-5', '0', '5'),
            cex = 0.9,
            height = 3
            ),
        legend = list(
            colours = c(brca.col, meta.col, ispy.col, matador.col, ccle.col),
            title = expression(underline('Datasets')),
            labels = c('TCGA-BRCA', 'METABRIC', 'I-SPY2', 'MATADOR', 'ICGC'),
            size = 2,
            label.cex = 1,
            continuous = FALSE,
            height = 3,
            padding.text = 2
            ),
        legend = list(
            colours = c(white),
            title = expression(underline('Status')),
            labels = c('Outlier'),
            size = 2,
            label.cex = 1,
            continuous = FALSE,
            height = 3,
            padding.text = 2
            )
        ),
    title.just = 'left',
    title.fontface = 'plain'
    );

heat.all <- BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(main.hetmap, sub.hetmap),
    x.relation = "sliced",
    y.relation = "sliced",
    main = expression('Outlier status'),
    xlab.label = expression('Cell line (CCLE)'),
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(2, 1),
    panel.widths = c(1, 0.22),
    ylab.padding = 2,
    xlab.to.xaxis.padding = -1.5,
    x.spacing = 0.3,
    main.cex = 1.5,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 1.3,
    xlab.cex = 1.1,
    xaxis.rot = 90,
    xaxis.tck = 0,
    legend = list(bottom = list(fun = legend.sample.grob)),
    print.new.legend = TRUE,
    resolution = 500
    );




# Save the multi plot as a PDF
pdf(
    file = generate.filename(
        'CCLE_outlier', 
        'multi', 
        'pdf'
        ), 
    width = 7.5, 
    height = 8.15
    );
heat.all;
dev.off();

# Save the multi plot as a PNG
png(
    file = generate.filename(
        'CCLE_outlier', 
        'multi', 
        'png'
        ), 
    width = 7.5, 
    height = 8.15,
    unit = 'in', 
    res = 1200
    );
heat.all;
dev.off();
