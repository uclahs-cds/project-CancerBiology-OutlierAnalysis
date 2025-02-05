### HISTORY #####################################################################
# This script visualizes the expression of outliers identified in CCLE across
# multiple datasets, overlapped with TCGA-BRCA, METABRIC, I-SPY2, MATADOR, and
# ICGC BRCA-EU.
# Date: 2024-08-16

### DESCRIPTION #################################################################
# This script processes and visualizes the expression of outliers identified in
# the Cancer Cell Line Encyclopedia (CCLE) dataset. It compares these outliers
# across multiple breast cancer datasets including TCGA-BRCA, METABRIC, I-SPY2,
# MATADOR, and ICGC BRCA-EU. The script generates a heatmap to visualize the
# outlier status and expression levels across different datasets.

### PREAMBLE ####################################################################
library(BoutrosLab.utilities);
library(RColorBrewer);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'ccle.sample.outlier.status.fdr.05.five',
    'ccle.sample.outlier.status.fdr.05.five.symbol',
    'outlier.symbol'
    ));

# overlapped with tissue outliers
outlier.patient.tag.01.brca.sum.overlap <- outlier.symbol$brca[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$brca)];
outlier.patient.tag.01.meta.sum.overlap <- outlier.symbol$metabric[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$metabric)];
outlier.patient.tag.01.ispy.sum.overlap <- outlier.symbol$ispy[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$ispy)];
outlier.patient.tag.01.matador.sum.overlap <- outlier.symbol$matador[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$matador)];
outlier.patient.tag.01.icgc.sum.overlap <- outlier.symbol$icgc[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$icgc)];
outlier.patient.tag.01.cheng.sum.overlap <- outlier.symbol$cheng[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$cheng)];
outlier.patient.tag.01.kao.sum.overlap <- outlier.symbol$kao[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$kao)];
outlier.patient.tag.01.hatzis.sum.overlap <- outlier.symbol$hatzis[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$hatzis)];
outlier.patient.tag.01.sjostrom.sum.overlap <- outlier.symbol$sjostrom[match(ccle.sample.outlier.status.fdr.05.five.symbol, outlier.symbol$sjostrom)];

ccle.overlap.outlier.05.five.tissue <- data.frame(cbind(
    meta = outlier.patient.tag.01.meta.sum.overlap,
    brca = outlier.patient.tag.01.brca.sum.overlap,
    ispy = outlier.patient.tag.01.ispy.sum.overlap,
    sjostrom = outlier.patient.tag.01.sjostrom.sum.overlap,
    cheng = outlier.patient.tag.01.cheng.sum.overlap,
    matador = outlier.patient.tag.01.matador.sum.overlap,
    icgc = outlier.patient.tag.01.icgc.sum.overlap,
    kao = outlier.patient.tag.01.kao.sum.overlap,
    hatzis = outlier.patient.tag.01.hatzis.sum.overlap
    ));
rownames(ccle.overlap.outlier.05.five.tissue) <- ccle.sample.outlier.status.fdr.05.five.symbol;

ccle.overlap.outlier.05.five.tissue.0.1 <- apply(ccle.overlap.outlier.05.five.tissue, 1, function(x) {
    x[is.na(x)] <- 0;
    ifelse(x > 0, 1, 0);
    });
ccle.overlap.outlier.05.five.tissue.0.1 <- t(ccle.overlap.outlier.05.five.tissue.0.1);


ccle.mean.zscore <- apply(fpkm.tumor.symbol.filter.ccle, 1, scale)
ccle.mean.zscore <- data.frame(t(ccle.mean.zscore));
colnames(ccle.mean.zscore) <- colnames(fpkm.tumor.symbol.filter.ccle);
ccle.mean.zscore.outlier.fdr.05.five <- ccle.mean.zscore[rownames(ccle.sample.outlier.status.fdr.05.five), ];

distance.matrix.t <- dist(t(ccle.mean.zscore.outlier.fdr.05.five), method = 'euclidean');
fit.t <- hclust(distance.matrix.t, method = 'ward.D2');
ccle.mean.zscore.outlier.fdr.05.five.t.order <- ccle.mean.zscore.outlier.fdr.05.five[, fit.t$order];

distance.matrix <- dist(ccle.mean.zscore.outlier.fdr.05.five.t.order, method = 'euclidean');
fit <- hclust(distance.matrix, method = 'ward.D2');
ccle.mean.zscore.outlier.fdr.05.five.t.p.order <- ccle.mean.zscore.outlier.fdr.05.five.t.order[fit$order, ];
colnames(ccle.mean.zscore.outlier.fdr.05.five.t.p.order) <- gsub('\\.', '-', colnames(ccle.mean.zscore.outlier.fdr.05.five.t.p.order));

ccle.overlap.outlier.05.five.tissue.0.1.order.zscore <- ccle.overlap.outlier.05.five.tissue.0.1[sub('\\..*', '', rownames(ccle.mean.zscore.outlier.fdr.05.five.t.p.order)), ];


library('RColorBrewer');
down.col <- brewer.pal(11, 'RdBu')[11];
up.col <- brewer.pal(11, 'RdBu')[1];

ccle.mean.zscore.outlier.fdr.05.five.t.p.order[ccle.mean.zscore.outlier.fdr.05.five.t.p.order > 5] <- 5;

# Row/col color - show the outlier status
ccle.sample.outlier.status.fdr.05.five.match <- ccle.sample.outlier.status.fdr.05.five[match(rownames(ccle.mean.zscore.outlier.fdr.05.five.t.p.order), rownames(ccle.sample.outlier.status.fdr.05.five)), match(colnames(ccle.mean.zscore.outlier.fdr.05.five.t.p.order), colnames(ccle.sample.outlier.status.fdr.05.five))];

main.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(ccle.mean.zscore.outlier.fdr.05.five.t.p.order),
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
    row.pos = which(t(ccle.sample.outlier.status.fdr.05.five.match) > 0, arr.ind = TRUE)[, 2] + 1.5,
    col.pos = which(t(ccle.sample.outlier.status.fdr.05.five.match) > 0, arr.ind = TRUE)[, 1],
    cell.text = rep('.', times = sum(t(ccle.sample.outlier.status.fdr.05.five.match) > 0)),
    text.col = 'white',
    text.cex = 1.5,
    axes.lwd = 0.8,
    colour.scheme = c(down.col, 'white', up.col),
    colour.centering.value = 0,
    at = seq(-5, 5, 0.001),
    colourkey.cex = 1,
    # covariate.legend = legend.sample,
    # legend.side = 'right',
    # legend.title.just = 'left',
    # legend.cex = 1,
    resolution = 1000
    );


ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore <- data.frame(ccle.overlap.outlier.05.five.tissue.0.1.order.zscore);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$brca <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$brca > 0, 3, 1);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$meta <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$meta > 0, 2, 1);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$ispy <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$ispy > 0, 4, 1);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$sjostrom <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$sjostrom > 0, 5, 1);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$cheng <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$cheng > 0, 6, 1);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$matador <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$matador > 0, 7, 1);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$icgc <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$icgc > 0, 8, 1);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$kao <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$kao > 0, 9, 1);
ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$hatzis <- ifelse(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore$hatzis > 0, 10, 1);


nine.col <- c(
    grDevices::adjustcolor('deepskyblue4', alpha.f = 0.7),
    grDevices::adjustcolor('firebrick3', alpha.f = 0.7),
    grDevices::adjustcolor('gold2', alpha.f = 0.7),
    grDevices::adjustcolor('darkgreen', alpha.f = 0.7),
    grDevices::adjustcolor('mediumpurple3', alpha.f = 0.7),
    grDevices::adjustcolor(c('darkorange2'), alpha.f = 0.7), 
    grDevices::adjustcolor(c('darkblue'), alpha.f = 0.7), 
    grDevices::adjustcolor(c('lightpink2'), alpha.f = 0.8),
    grDevices::adjustcolor(c('yellowgreen'), alpha.f = 0.8)
    );
all.col <- c('grey95', nine.col);

sub.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(ccle.overlap.outlier.05.five.tissue.0.1.order.col.zscore),
    clustering.method = 'none',
    colour.scheme = all.col,
    total.colours = 11,
    row.colour = 'black',
    col.colour = 'black',
    grid.row = TRUE,
    grid.col = TRUE,
    yaxis.tck = 0,
    xaxis.tck = 0,
    xaxis.lab = rev(c( 'Hatzis', 'Kao', 'ICGC\n BRCA-EU', 'MATADOR', 'Cheng', 'Sjostrom', 'ISPY-2', 'TCGA-BRCA', 'METABRIC')),
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
            colours = c(down.col, 'white', up.col),
            total.colours = 100,
            labels = c('-5', '0', '5'),
            cex = 0.9,
            height = 3
            ),
        legend = list(
            colours = nine.col,
            title = expression(underline('Datasets')),
            labels = rev(c( 'Hatzis', 'Kao', 'ICGC BRCA-EU', 'MATADOR', 'Cheng', 'Sjostrom', 'ISPY-2', 'TCGA-BRCA', 'METABRIC')),
            size = 2,
            label.cex = 1,
            continuous = FALSE,
            height = 3,
            padding.text = 2
            ),
        legend = list(
            colours = c('white'),
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
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = expression('Outlier status'),
    xlab.label = expression('Cell line (CCLE)'),
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(2, 1),
    panel.widths = c(1, 0.3),
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
    legend = list(right = list(fun = legend.sample.grob)),
    print.new.legend = TRUE,
    resolution = 500
    );


save.outlier.figure(
    heat.all,
    c('Figure4a', 'CCLE', 'outlier', 'multi'),
    width = 7.5,
    height = 8.15
    );

save.session.profile(file.path('output', 'Figure4a.txt'));
