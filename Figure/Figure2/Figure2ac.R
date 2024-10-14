### HISTORY ######################################################################
# This script performs CNA analysis for METABRIC/TCGA-BRCA/ICGC BRCA-EU data.
# Date: 2024-08-23

### DESCRIPTION ##################################################################
# This script conducts Copy Number Aberration (CNA) analysis across three breast
# cancer datasets: METABRIC, TCGA-BRCA, and ICGC BRCA-EU. It processes CNA data
# for outlier and non-outlier samples, performs meta-analysis, and creates
# visualizations to represent the results.

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'all.weighted.mean',
    'brca.cnv.chr.new.gis.fpkm.order.match',
    'brca.cnv.chr.new.gis.match',
    'brca.gis.table.num',
    'icgc.cnv.chr.new.gis.fpkm.order.match',
    'icgc.cnv.chr.new.gis.fpkm.order.match.chr',
    'icgc.gis.table.num',
    'meta.cnv.chr.new.gis.match',
    'meta.gis.table.num',
    'metafor.cnv.lower.all',
    'metafor.cnv.odd.all',
    'metafor.cnv.p.all',
    'metafor.cnv.p.all.fdr',
    'metafor.cnv.upper.all',
    'outlier.patient.tag.01.brca.cnv.match',
    'outlier.patient.tag.01.meta.cnv.match',
    'outlier.symbol'
    ));






### PLOTTING RESULTS ############################################################
# - Use weighted mean
# Initialize an empty vector to store the weighted means
# Initialize vectors to store the weighted sums
brca.sum.vector <- sapply(unique(brca.gis.table.num$new.gis.value), function(x) {
    sum(brca.gis.table.num$sample[brca.gis.table.num$new.gis.value == x])
    })
brca.sum.vector <- rep(brca.sum.vector, each = 2);

meta.sum.vector <- sapply(unique(meta.gis.table.num$new.gis.value), function(x) {
    sum(meta.gis.table.num$sample[meta.gis.table.num$new.gis.value == x])
    })
meta.sum.vector <- rep(meta.sum.vector, each = 2);

icgc.sum.vector <- sapply(unique(icgc.gis.table.num$new.gis.value), function(x) {
    sum(icgc.gis.table.num$sample[icgc.gis.table.num$new.gis.value == x])
    })
icgc.sum.vector <- rep(icgc.sum.vector, each = 2);

# Initialize vector to store the weighted means
all.weighted.mean <- numeric(10);

# Compute weighted means for each status and GIS value combination
for (status in c('non', 'out')) {
    for (i in 1:5) {
        gis.value <- c(-2, -1, 0, 1, 2)[i];

        weights <- c(
            sum(brca.gis.table.num$sample[brca.gis.table.num$new.gis.value == gis.value & brca.gis.table.num$status == status]),
            sum(meta.gis.table.num$sample[meta.gis.table.num$new.gis.value == gis.value & meta.gis.table.num$status == status]),
            sum(icgc.gis.table.num$sample[icgc.gis.table.num$new.gis.value == gis.value & icgc.gis.table.num$status == status])
            );

        values <- c(
            brca.gis.table.num$Freq[brca.gis.table.num$new.gis.value == gis.value & brca.gis.table.num$status == status],
            meta.gis.table.num$Freq[meta.gis.table.num$new.gis.value == gis.value & meta.gis.table.num$status == status],
            icgc.gis.table.num$Freq[icgc.gis.table.num$new.gis.value == gis.value & icgc.gis.table.num$status == status]
            );

        # Compute the weighted average and store it
        index <- ifelse(status == 'non', i, i + 5);
        all.weighted.mean[index] <- weighted.mean(values, weights);
        }
    }



# 1. Create non-outlier fraction bar plot
metafor.cnv.odd.all.df.fraction.non <- data.frame(
    group = rev(c('a', 'b', 'c', 'd', 'e')),
    sample = rep('a', 5),
    odd = all.weighted.mean[1:5]
    );

cnv.col.ramp.5.seg <- c('#2166AC', '#90B2D5', 'grey60', '#D88B95', '#B2182B')
cnv.col.ramp.5.bar <- c('#2166AC', '#90B2D5', 'white', '#D88B95', '#B2182B')

non.fraction.bar <- BoutrosLab.plotting.general:::create.barplot(
    formula = group ~ odd,
    data = metafor.cnv.odd.all.df.fraction.non,
    main = NULL,
    plot.horizontal = TRUE,
    col = cnv.col.ramp.5.bar,
    stack = FALSE,
    main.cex = 0,
    border.col = 'black',
    border.lwd = 1,
    ylab.label = NULL,
    xlab.cex = 1.1,
    xlab.label = expression('Fraction'),
    xlimits = c(-0.07, 0.82),
    ylab.cex = 1.3,
    xaxis.cex = 1.1,
    yaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.rectangle = TRUE,
    xleft.rectangle = -5,
    xright.rectangle = 10,
    ybottom.rectangle = c(1.5, 3.5),
    ytop.rectangle = c(2.5, 4.5),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    );
non.fraction.bar;



# 2. outlier fraction
metafor.cnv.odd.all.df.fraction.out <- data.frame(
    group = rev(c('a', 'b', 'c', 'd', 'e')),
    sample = rep('a', 5),
    odd = all.weighted.mean[6:10]
    );

out.fraction.bar <- BoutrosLab.plotting.general:::create.barplot(
    formula = group ~ odd,
    data = metafor.cnv.odd.all.df.fraction.out,
    main = NULL,
    plot.horizontal = TRUE,
    col = cnv.col.ramp.5.bar,
    stack = FALSE,
    main.cex = 0,
    border.col = 'black',
    border.lwd = 1,
    # main.just = 'right',
    ylab.label = NULL,
    xlab.cex = 1.1,
    xlab.label = expression('Fraction'),
    xlimits = c(-0.07, 0.82),
    ylab.cex = 1.3,
    xaxis.cex = 1.1,
    yaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    # xaxis.rot = 45,
    add.rectangle = TRUE,
    xleft.rectangle = -5,
    xright.rectangle = 10,
    ybottom.rectangle = c(1.5, 3.5),
    ytop.rectangle = c(2.5, 4.5),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    );
out.fraction.bar;



# 3. segment plot
metafor.cnv.odd.ci.p.data <- data.frame(cbind(
    p.value = metafor.cnv.p.all,
    odd = metafor.cnv.odd.all,
    ci.min = metafor.cnv.lower.all,
    ci.max = metafor.cnv.upper.all,
    fdr = metafor.cnv.p.all.fdr
    ));

metafor.cnv.odd.ci.p.data$label <- rev(c('a', 'b', 'c', 'd', 'e'));
metafor.cnv.odd.ci.p.data$label <- as.factor(metafor.cnv.odd.ci.p.data$label);

# Change the order from lowest to higest odds from meta analysis
metafor.cnv.odd.ci.p.data.rev <- metafor.cnv.odd.ci.p.data[rev(1:5), ];

metafor.all.segplot <- BoutrosLab.plotting.general::create.segplot(
    formula = label ~ log2(ci.min) + log2(ci.max),
    data = metafor.cnv.odd.ci.p.data.rev,
    # add middle dots
    centers = log2(metafor.cnv.odd.ci.p.data.rev$odd),
    main.cex = 0,
    yaxis.fontface = 1,
    xlab.cex = 1.3,
    xlab.label = expression('Odds Ratio'),
    symbol.cex = 1.2,
    ylab.cex = 0,
    yaxis.cex = 0,
    xaxis.cex = 1,
    xaxis.lab = c(0, 0.25, 0.5, 1, 2, 4, 8, 16),
    xaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    segments.col = rev(cnv.col.ramp.5.seg),
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -3,
    xright.rectangle = 13,
    ybottom.rectangle = seq(1.5, 23.5, 2),
    ytop.rectangle = seq(2.5, 24.5, 2),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
    );
metafor.all.segplot;


# 4. FDR
metafor.cnv.fdr <- data.frame(
    group = rev(c('a', 'b', 'c', 'd', 'e')),
    sample = rep('a', 5),
    odd = metafor.cnv.odd.ci.p.data$fdr
    );

fdr.bar <- BoutrosLab.plotting.general:::create.barplot(
    formula = group ~ -log10(odd),
    data = metafor.cnv.fdr,
    main = NULL,
    plot.horizontal = TRUE,
    col = 'grey20',
    stack = FALSE,
    main.cex = 0,
    border.col = 'black',
    border.lwd = 1,
    ylab.label = NULL,
    xlab.cex = 1.1,
    xlab.label = expression('-log'[10] * '(FDR)'),
    xlimits = c(-1, 13),
    ylab.cex = 1.3,
    xaxis.cex = 1.1,
    yaxis.cex = 0,
    abline.v = 2,
    abline.lty = 3,
    abline.lwd = 1.5,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    # xaxis.rot = 45,
    add.rectangle = TRUE,
    xleft.rectangle = -5,
    xright.rectangle = 15,
    ybottom.rectangle = c(1.5, 3.5),
    ytop.rectangle = c(2.5, 4.5),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    );
fdr.bar;

multi.gene <- create.multipanelplot(
    list(non.fraction.bar, out.fraction.bar, metafor.all.segplot, fdr.bar),
    main.cex = 0,
    main.y = 0.5,
    resolution = 300,
    layout.height = 1,
    layout.width = 4,
    layout.skip = c(FALSE, FALSE, FALSE, FALSE),
    plot.objects.widths = c(1, 1, 1, 1),
    ylab.axis.padding = -5,
    xlab.cex = 0,
    xlab.axis.padding = 1,
    x.spacing = 1,
    right.legend.padding = 0
    );

multi.gene;

save.outlier.figure(
    multi.gene,
    c('Figure2ac', 'CNA', 'multipanel'),
    width = 10.4,
    height = 4.5
    );

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
meta.cnv.chr.new.gis.fpkm <- meta.cnv.chr.new.gis[
    meta.cnv.chr.new.gis$Hugo_Symbol %in%
        fpkm.tumor.symbol.filter.meta.symbol$Symbol,
    ];
meta.cnv.chr.new.gis.fpkm.order <- meta.cnv.chr.new.gis.fpkm[
    na.omit(match(
        meta.cnv.chr.new.gis.fpkm.order.match.chr$gene_name,
        meta.cnv.chr.new.gis.fpkm$Hugo_Symbol
        )),
    ];
meta.cnv.chr.new.gis.fpkm.order.match <- meta.cnv.chr.new.gis.fpkm.order[
    ,
    c('Hugo_Symbol', colnames(meta.cnv.chr.new.gis.match)[
        colnames(meta.cnv.chr.new.gis.match) %in%
            substr(colnames(outlier.patient.tag.01.meta), 1, 15)
        ])
    ];
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
    c('Figure2ac', 'CNA', 'chr10', 'multipanel'),
    width = 10.4,
    height = 4.5
    );
save.session.profile(file.path('output', 'Figure2ac.txt'));
