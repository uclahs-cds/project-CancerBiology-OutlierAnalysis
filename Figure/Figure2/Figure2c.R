### HISTORY ######################################################################
# This script processes CNA data for chromosome 10 across multiple datasets.
# It compares outlier and non-outlier samples, generates heatmaps, and
# combines these plots into a single multiplot figure.
# Date: 2024-08-13

### PREAMBLE ####################################################################
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-08-23_Figure2a-d.rda'));



### DATA EXTRACTION #############################################################

# 1. TCGA-BRCA
# Filtering data for chromosome 10
brca.cnv.chr.new.gis.fpkm.order.match.chr10 <- brca.cnv.chr.new.gis.fpkm.order.match[
    brca.cnv.chr.new.gis.fpkm.order.match.chr$chromosome == 'chr10',
    ];

# Extracting outlier patient samples for FGFR2
brca.out.sample <- colnames(outlier.patient.tag.01.brca.cnv.match)[
    outlier.patient.tag.01.brca.cnv.match[
        rownames(fpkm.tumor.symbol.filter.brca[
            fpkm.tumor.symbol.filter.brca$Symbol %in% 'FGFR2',
            ]),
        ] == 1
    ];

# Extracting non-outlier patient samples for FGFR2
brca.non.out.sample <- colnames(outlier.patient.tag.01.brca.cnv.match)[
    outlier.patient.tag.01.brca.cnv.match[
        rownames(fpkm.tumor.symbol.filter.brca[
            fpkm.tumor.symbol.filter.brca$Symbol %in% 'FGFR2',
            ]),
        ] == 0
    ];

# Filtering data for non-outlier patient samples in chromosome 10
brca.cnv.chr.new.gis.fpkm.order.match.chr10.out.non <- brca.cnv.chr.new.gis.fpkm.order.match.chr10[
    , 2:ncol(brca.cnv.chr.new.gis.fpkm.order.match.chr10)
    ][
    , colnames(brca.cnv.chr.new.gis.fpkm.order.match.chr10)[2:ncol(brca.cnv.chr.new.gis.fpkm.order.match.chr10)] %in% substr(brca.non.out.sample, 1, 15)
    ];



# 2. METABRIC
# Filtering data for chromosome 10
meta.cnv.chr.new.gis.fpkm.order.match.chr10 <- meta.cnv.chr.new.gis.fpkm.order.match[
    meta.cnv.chr.new.gis.fpkm.order.match.chr$chromosome == 'chr10',
    ];

# Extracting outlier patient samples for FGFR2
meta.out.sample <- colnames(outlier.patient.tag.01.meta.cnv.match)[
    outlier.patient.tag.01.meta.cnv.match[
        rownames(fpkm.tumor.symbol.filter.meta.symbol[
            fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% 'FGFR2',
            ]),
        ] == 1
    ];

# Extracting non-outlier patient samples for FGFR2
meta.non.out.sample <- colnames(outlier.patient.tag.01.meta.cnv.match)[
    outlier.patient.tag.01.meta.cnv.match[
        rownames(fpkm.tumor.symbol.filter.meta.symbol[
            fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% 'FGFR2',
            ]),
        ] == 0
    ];

# Filtering data for non-outlier patient samples in chromosome 10
meta.cnv.chr.new.gis.fpkm.order.match.chr10.out.non <- meta.cnv.chr.new.gis.fpkm.order.match.chr10[
    , 2:ncol(meta.cnv.chr.new.gis.fpkm.order.match.chr10)
    ][
    , colnames(meta.cnv.chr.new.gis.fpkm.order.match.chr10)[
        2:ncol(meta.cnv.chr.new.gis.fpkm.order.match.chr10)
        ] %in% substr(meta.non.out.sample, 1, 15)
    ];


# 3. ICGC-BRCA_EU
# Extracting gene symbols from raw data
icgc.cnv.all.symbol <- sub('\\|.*', '', icgc.cnv.chr.new.gis.raw$Gene.Symbol);

# Filtering data for chromosome 10
icgc.cnv.chr.new.gis.fpkm.order.match.chr10 <- icgc.cnv.chr.new.gis.fpkm.order.match[
    icgc.cnv.chr.new.gis.fpkm.order.match.chr == '10',
    ];

# Filtering gene symbols for chromosome 10
icgc.cnv.all.symbol.10 <- icgc.cnv.all.symbol[
    icgc.cnv.chr.new.gis.fpkm.order.match.chr == '10'
    ];

# Extracting outlier patient samples for FGFR2 in ICGC data
icgc.out.sample <- colnames(outlier.patient.tag.01.icgc)[
    outlier.patient.tag.01.icgc[
        rownames(fpkm.data.icgc)[
            fpkm.data.icgc$Name %in% 'FGFR2'
            ],
        ] == 1
    ];

# Filtering data for non-outlier patient samples in chromosome 10
icgc.cnv.chr.new.gis.fpkm.order.match.chr10.out.non <- icgc.cnv.chr.new.gis.fpkm.order.match.chr10[
    , !(colnames(icgc.cnv.chr.new.gis.fpkm.order.match.chr10) %in% icgc.out.sample)
    ];

# 1. Outlier Data
brca.chr10.out <- brca.cnv.chr.new.gis.fpkm.order.match.chr10[
    , substr(brca.out.sample, 1, 15),
    drop = FALSE
    ];
brca.chr10.out.symbol <- brca.cnv.chr.new.gis.fpkm.order.match.chr10$Hugo_Symbol;

meta.chr10.out <- meta.cnv.chr.new.gis.fpkm.order.match.chr10[
    , substr(na.omit(meta.out.sample), 1, 15)
    ];
meta.chr10.out.symbol <- meta.cnv.chr.new.gis.fpkm.order.match.chr10$Hugo_Symbol;

icgc.chr10.out <- icgc.cnv.chr.new.gis.fpkm.order.match.chr10[
    , icgc.out.sample,
    drop = FALSE
    ];
icgc.chr10.out.symbol <- icgc.cnv.all.symbol.10;

unique.chr10.symbol <- intersect(
    brca.chr10.out.symbol,
    intersect(meta.chr10.out.symbol, icgc.chr10.out.symbol)
    );

# Matching Outlier Data
brca.chr10.out.match <- brca.chr10.out[
    brca.chr10.out.symbol %in% unique.chr10.symbol,
    ];
brca.chr10.out.match <- brca.chr10.out.match[
    !duplicated(brca.chr10.out.symbol[brca.chr10.out.symbol %in% unique.chr10.symbol]),
    ];
meta.chr10.out.match <- meta.chr10.out[
    meta.chr10.out.symbol %in% unique.chr10.symbol,
    ];
meta.chr10.out.match <- meta.chr10.out.match[
    !duplicated(meta.chr10.out.symbol[meta.chr10.out.symbol %in% unique.chr10.symbol]),
    ];
icgc.chr10.out.match <- icgc.chr10.out[
    icgc.chr10.out.symbol %in% unique.chr10.symbol,
    ];

all.chr10.out.match <- cbind(
    brca.chr10.out.match,
    meta.chr10.out.match,
    icgc.chr10.out.match
    );


# 2. Non-Outlier Data
brca.chr10.non.out.match <- brca.cnv.chr.new.gis.fpkm.order.match.chr10.out.non[
    brca.chr10.out.symbol %in% unique.chr10.symbol,
    ];
brca.chr10.non.out.match <- brca.chr10.non.out.match[
    !duplicated(brca.chr10.out.symbol[brca.chr10.out.symbol %in% unique.chr10.symbol]),
    ];
meta.chr10.non.out.match <- meta.cnv.chr.new.gis.fpkm.order.match.chr10.out.non[
    meta.chr10.out.symbol %in% unique.chr10.symbol,
    ];
meta.chr10.non.out.match <- meta.chr10.non.out.match[
    !duplicated(meta.chr10.out.symbol[meta.chr10.out.symbol %in% unique.chr10.symbol]),
    ];
icgc.chr10.non.out.match <- icgc.cnv.chr.new.gis.fpkm.order.match.chr10.out.non[
    icgc.chr10.out.symbol %in% unique.chr10.symbol,
    ];

all.chr10.non.out.match <- cbind(
    brca.chr10.non.out.match,
    meta.chr10.non.out.match,
    icgc.chr10.non.out.match
    );

### CLUSTERING AND ORDERING
distance.matrix.t.all.chr10.non.out.match <- dist(t(all.chr10.non.out.match), method = 'euclidean');
fit.t.all.chr10.non.out.match <- hclust(distance.matrix.t.all.chr10.non.out.match, method = 'ward.D2');
all.chr10.non.out.match.order <- all.chr10.non.out.match[, fit.t.all.chr10.non.out.match$order];

### HEATMAP
chr10.out <- create.heatmap(
    x = all.chr10.out.match,
    clustering.method = 'none',
    colour.scheme = c('#2166ac', 'white', '#b2182b'),
    grid.row = FALSE,
    grid.col = FALSE,
    yaxis.tck = 0,
    xaxis.tck = 0,
    yaxis.cex = 1.5,
    yaxis.rot = 0,
    ylab.cex = 0,
    colour.centering.value = 0,
    at = seq(-2, 2, 0.1),
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );

chr10.non <- create.heatmap(
    x = all.chr10.non.out.match.order,
    clustering.method = 'none',
    colour.scheme = c('#2166ac', 'white', '#b2182b'),
    grid.row = FALSE,
    grid.col = FALSE,
    yaxis.tck = 0,
    xaxis.tck = 0,
    yaxis.cex = 1.5,
    yaxis.rot = 0,
    ylab.cex = 0,
    colour.centering.value = 0,
    at = seq(-2, 2, 0.1),
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );

### MULTIPLOT
cnv.col <- c('#2166ac', 'white', '#b2182b');
cnv.col.ramp <- colorRampPalette(cnv.col);
cnv.col.ramp.5 <- cnv.col.ramp(5);

cna.multi <- create.multiplot(
    plot.objects = list(chr10.non, chr10.out),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = expression('CNA of chromosome 10'),
    xlab.label = expression('Genes on chromosome 10'),
    ylab.label = c(
        expression('Outlier patients'), '', '', '',
        expression('Non-outlier patients'), '', '', '', ''
        ),
    yaxis.fontface = 1,
    plot.layout = c(1, 2),
    main.key.padding = 3,
    panel.heights = c(0.12, 1),
    ylab.padding = 1,
    y.spacing = -0.7,
    main.cex = 1.6,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    xlab.padding = -10,
    xlab.to.xaxis.padding = -1,
    right.padding = 3,
    bottom.padding = 10,
    # Setting groups
    legend = list(
        bottom = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = 'black',
                        pch = 22,
                        cex = 2.5,
                        fill = cnv.col.ramp.5
                        ),
                    text = list(
                        lab = c('Homozygous deletion', 'Hemizygous deletion', 'Neutral', 'Gain', 'High level amplification')
                        ),
                    padding.text = 3,
                    cex = 1,
                    just = 'left'
                    )
                )
            )
        ),
    print.new.legend = TRUE,
    yaxis.cex = 1.2,
    yaxis.tck = 0,
    ylab.cex = 1.15,
    xlab.cex = 1.3,
    xaxis.rot = 90,
    xaxis.tck = 0,
    xlab.key.padding = 9,
    resolution = 500
    );

save.outlier.figure(
    cna.multi,
    c('Figure2c', 'CNA', 'chr10', 'multipanel'),
    width = 10.4,
    height = 4.5
    );

save.session.profile(file.path('output', 'Figure2c.txt'));
