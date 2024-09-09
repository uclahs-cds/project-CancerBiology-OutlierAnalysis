### HISTORY #####################################################################
# This script analyzes the DNA methylation on the promoter region of the PDXNL 
# gene in outlier and non-outlier patients in TCGA-BRCA data. It generates heatmaps to visualize 
# the methylation levels across these patient groups.
# Date: 2024-08-14



# Load necessary library
library(BoutrosLab.plotting.general);


# Check PXDNL gene
i <- 'PXDNL';
i.me <- brca.outlier.promoter.symbol.sample.match.brca[
    brca.outlier.promoter.symbol.sample.match.brca$Symbol == "PXDNL",
    1:ncol(brca.outlier.promoter.symbol.sample.match.brca)
    ];


i.me.normal <- brca.outlier.promoter.symbol.normal.match.filter.brca[
    brca.outlier.promoter.symbol.sample.match.brca[rownames(brca.outlier.promoter.symbol.normal.match.filter.brca),]$Symbol == "PXDNL",
    ];



# i.me <- brca.outlier.promoter.symbol.sample.match.merge.filter[
#     rownames(brca.outlier.promoter.symbol.sample.match.merge.filter) == i, 
#     1:ncol(outlier.patient.tag.01.t.p.order.me.sample.match.gene.sum.filter.brca)
#     ];

i.patient <- outlier.patient.tag.01.t.p.order.me.sample.match.gene.sum.filter.brca[
    rownames(fpkm.tumor.symbol.filter.brca)[fpkm.tumor.symbol.filter.brca$Symbol == i], 
    ];

i.me.mean <- apply(
    i.me[, 1:ncol(outlier.patient.tag.01.t.p.order.me.sample.match.gene.sum.filter.brca)], 
    1, 
    function(x) {mean(na.omit(x));}
    );

i.me.patient <- i.me[, which(i.patient == 1), drop = FALSE];
i.me.patient.non <- i.me[, which(i.patient == 0), drop = FALSE];


overlap.patient <- colnames(brca.outlier.promoter.symbol.normal.match.filter.brca)[
    substr(colnames(brca.outlier.promoter.symbol.normal.match.filter.brca), 1, 12) %in% 
        substr(names(i.patient)[i.patient == 1], 1, 12)
    ];


i.me.patient.normal <- i.me.normal[, overlap.patient, drop = FALSE];
i.me.patient.non.normal <- i.me.normal[, !(colnames(i.me.normal) %in% overlap.patient), drop = FALSE];

i.me.patient.non.mean <- apply(
    i.me.patient.non, 
    1, 
    function(x) {mean(na.omit(x));}
    );

# names(i.me.patient) <- names(i.me.mean);
promoters.i <- promoters.info[[i]];



promoters.i.order <- promoters.i[order(promoters.i$pos), ];
i.me.mean.order <- i.me.patient[rownames(promoters.i[order(promoters.i$pos), 1:7]),, drop = FALSE];
i.me.patient.order <- i.me.patient.non.mean[rownames(promoters.i[order(promoters.i$pos), 1:7])];

fpkm.i <- fpkm.tumor.symbol.filter.brca[
    fpkm.tumor.symbol.filter.brca$Symbol == i, 
    , 
    drop = FALSE
    ];

fpkm.i.order <- fpkm.i[
    , 
    colnames(outlier.patient.tag.01.t.p.order.me.sample.match.gene.sum.filter.brca), 
    drop = FALSE
    ];

fpkm.i.order <- fpkm.i.order[
    , 
    order(as.numeric(fpkm.i.order[1,]), decreasing = TRUE), 
    drop = FALSE
    ];

i.me.order <- i.me[, colnames(fpkm.i.order)];

# Location
promoters.i.order <- promoters.i[order(promoters.i$pos), ];

# Heatmap
i.me.beta.order <- i.me.patient[
    rownames(promoters.i.order[order(promoters.i.order$pos), 1:7]),, drop = FALSE
    ];

i.me.beta.order.non <- i.me.patient.non.mean[
    rownames(promoters.i.order[order(promoters.i.order$pos), 1:7])
    ];

i.me.beta.order.normal <- i.me.patient.normal[
    rownames(promoters.i.order[order(promoters.i.order$pos), 1:7]),, drop = FALSE
    ];

i.me.beta.order.non.normal <- i.me.patient.non.normal[
    rownames(promoters.i.order[order(promoters.i.order$pos), 1:7]),
    ];


i.me.patient.non.order <- i.me.patient.non[
    rownames(promoters.i.order[order(promoters.i.order$pos), 1:7]), 
    ];

col.key <- c("#b2182b", "white", "#2166ac");

i.heat.out <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(i.me.beta.order),
    clustering.method = 'none',
    colour.scheme = col.key, 
    grid.row = FALSE, 
    grid.col = FALSE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.cex = 0,
    yaxis.rot = 0,
    ylab.cex = 0,
    at = seq(0, 1, 0.001),
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );
i.heat.out.normal <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(i.me.beta.order.normal),
    clustering.method = 'none',
    colour.scheme = col.key, 
    grid.row = FALSE, 
    grid.col = FALSE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.cex = 0,
    yaxis.rot = 0,
    ylab.cex = 0,
    at = seq(0, 1, 0.001),
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );


# All non-outlier patients
i.heat.out.non <- BoutrosLab.plotting.general:::create.heatmap(
    x = i.me.patient.non.order,
    clustering.method = 'ward.D2',
    cluster.dimensions = 'row',
    plot.dendrograms = FALSE,
    colour.scheme = col.key, 
    grid.row = FALSE, 
    grid.col = FALSE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.cex = 0,
    yaxis.rot = 0,
    ylab.cex = 0,
    at = seq(0, 1, 0.001),
    colourkey.cex = 1.3,
    print.colour.key = TRUE
    );
i.heat.out.non.normal <- BoutrosLab.plotting.general:::create.heatmap(
    x = i.me.beta.order.non.normal,
    clustering.method = 'ward.D2',
    cluster.dimensions = 'row',
    plot.dendrograms = FALSE,
    colour.scheme = col.key, 
    grid.row = FALSE, 
    grid.col = FALSE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.cex = 0,
    yaxis.rot = 0,
    ylab.cex = 0,
    at = seq(0, 1, 0.001),
    colourkey.cex = 1.3,
    print.colour.key = TRUE
    );



i.heat  <- BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(i.heat.out.non.normal, i.heat.out.non, i.heat.out.normal, i.heat.out),
    x.relation = "sliced",
    y.relation = "sliced",
    main = i,
    xlab.label = expression('Beta value'),
    ylab.label = c(expression('Outlier patient'), "", "", "", expression('Non-outlier patients'), "", "", "", ""),
    yaxis.fontface = 1,
    plot.layout = c(1, 4),
    main.key.padding = 3,
    panel.heights = c(0.08, 0.08, 1, 0.2),
    ylab.padding = 1,
    y.spacing = -0.7,
    main.cex = 1.6,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    xlab.padding = -10, 
    xlab.to.xaxis.padding = -1, 
    bottom.padding = 3, 
    right.padding = 3, 
    yaxis.cex = 1.2,
    yaxis.tck = 0,
    ylab.cex = 1.15,
    xlab.cex = 1.3,
    xaxis.rot = 90,
    xaxis.tck = 0,
    xlab.key.padding = 4,
    resolution = 500
    );



# Save the heatmap as a PDF
pdf(
    file = generate.filename(
        i, 
        'heatmap_me', 
        'pdf'
        ), 
    width = 7.5, 
    height = 8.5
    );
i.heat;
dev.off();

# Save the heatmap as a PNG
png(
    file = generate.filename(
        i, 
        'heatmap_me', 
        'png'
        ), 
    width = 7.5, 
    height = 8.5,
    unit = 'in', 
    res = 1200
    );
i.heat;
dev.off();

