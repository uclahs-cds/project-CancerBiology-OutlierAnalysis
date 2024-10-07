### HISTORY #####################################################################
# This script generates a multi-panel plot for RNA, protein, CNV, and gene
# effect scores (RNAi and CRISPR-Cas) for a specific gene.
# Date: 2024-08-16

### DESCRIPTION #################################################################
# The script processes various data types (RNA, protein, CNV, RNAi, and CRISPR-Cas)
# for a specific gene of interest. It creates individual plots for
# each data type: a barplot for RNA abundance, a barplot for protein abundance with
# handling for missing values, a heatmap for CNV data, and barplots for gene effect
# scores from RNAi and CRISPR-Cas datasets.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-10-03_Figure1_2_3_4_min_input.rda'));

load.multiple.computed.variables(c(
    'protein.info.breast.num.symbol',
    'ccle.sample.outlier.status.overlap.na',
    'rnai.05.box',
    'rnai.score.05.overlap.minus.05',
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


i <- 'FOXP4';

# Prepare RNA abundance data
i.fpkm <- fpkm.tumor.symbol.filter.ccle[i, ];
i.fpkm.data <- data.frame(
    gene = as.numeric(i.fpkm),
    sample = colnames(fpkm.tumor.symbol.filter.ccle)
    );
i.fpkm.data.order <- i.fpkm.data[order(i.fpkm.data$gene, decreasing = TRUE), ];

# Prepare protein abundance data
i.protein <- protein.info.breast.num[protein.info.breast.num.symbol %in% i, match(i.fpkm.data.order$sample, colnames(protein.info.breast.num))];


# Create RNA abundance barplot
bar.fpkm <- BoutrosLab.plotting.general::create.barplot(
    formula = gene ~ sample,
    main = NULL,
    ylab.label = expression('RNA abundance (TPM)'),
    xlab.label = NULL,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0, 0),
    xaxis.cex = 0,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.2,
    main.cex = 1.4,
    col = c('red3', rep('grey20', (nrow(i.fpkm.data) - 1))),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    data = i.fpkm.data,
    sample.order = 'decreasing'
    );

# Prepare protein data with NA values handled
i.protein.na <- i.protein;
i.protein.na[is.na(i.protein.na)] <- 0;
i.protein.na.df <- data.frame(
    gene = i.protein.na,
    sample = i.fpkm.data.order$sample,
    order = c(1:length(i.protein.na))
    );

# Create protein abundance barplot
na.start <- as.numeric(which(is.na(i.protein))) - 0.5;
na.end <- na.start + 1;

bar.protein.na <- BoutrosLab.plotting.general::create.barplot(
    formula = gene ~ order,
    main = NULL,
    ylab.label = expression('Protein abundance'),
    xlab.label = NULL,
    xaxis.lab = i.protein.na.df$sample,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0, 0),
    xaxis.cex = 0,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.2,
    main.cex = 1.4,
    col = c('red3', rep('grey20', (nrow(i.fpkm.data) - 1))),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    data = i.protein.na.df,
    sample.order = 'none',
    add.rectangle = TRUE,
    xleft.rectangle = na.start,
    xright.rectangle = na.end,
    ybottom.rectangle = -6,
    ytop.rectangle = 10,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25
    );

# Prepare CNV data for heatmap
i.cnv <- cnv.info.breast.t.num[rownames(i.fpkm), i.fpkm.data.order$sample];
max.lim.cnv <- ceiling(max(na.omit(i.cnv)));

# Create CNV heatmap
cnv.plot <- BoutrosLab.plotting.general::create.heatmap(
    x = t(data.frame(i.cnv)),
    clustering.method = 'none',
    colour.scheme = c('dodgerblue4', 'white', 'red4'),
    grid.row = FALSE,
    grid.col = FALSE,
    yaxis.tck = 0,
    xaxis.tck = 0,
    axes.lwd = 0.5,
    yaxis.lab = c('', '', '', 'CNV', '', '', ''),
    xlab.label = expression('Breast cancer cell lines (CCLE)'),
    yaxis.fontface = 1,
    xaxis.fontface = 1,
    yaxis.cex = 1,
    xaxis.cex = 1,
    yaxis.rot = 0,
    ylab.cex = 0,
    xlab.cex = 1,
    colour.centering.value = 1,
    at = seq(0, 2.1, 0.001),
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );

# Prepare RNAi data for barplot
rnai.05.box.4.FOXP4 <- rnai.05.box[rnai.05.box$gene %in% i, ];
rnai.05.box.4.FOXP4$label <- rep('RNAi', nrow(rnai.05.box.4.FOXP4));
rownames(rnai.05.box.4.FOXP4) <- colnames(rnai.score.05.overlap.minus.05);
rnai.05.box.4.FOXP4.order <- rnai.05.box.4.FOXP4[match(i.protein.na.df$sample, rownames(rnai.05.box.4.FOXP4)), ];
rnai.05.box.4.FOXP4.order$order <- 1:nrow(rnai.05.box.4.FOXP4.order);

na.start.rnai <- as.numeric(which(is.na(rnai.05.box.4.FOXP4.order$score))) - 0.5;
na.end.rnai <- na.start.rnai + 1;

# Create RNAi gene effect score barplot
bar.rnai.na <- BoutrosLab.plotting.general::create.barplot(
    formula = score ~ order,
    main = NULL,
    ylab.label = expression('Gene effect score (RNAi)'),
    xlab.label = NULL,
    xaxis.lab = rownames(rnai.05.box.4.FOXP4.order),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0, 0),
    xaxis.cex = 0,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.2,
    main.cex = 1.4,
    col = c('red3', rep('grey20', (nrow(i.fpkm.data) - 1))),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    data = rnai.05.box.4.FOXP4.order,
    sample.order = 'none',
    add.rectangle = TRUE,
    xleft.rectangle = na.start.rnai,
    xright.rectangle = na.end.rnai,
    ybottom.rectangle = -6,
    ytop.rectangle = 10,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25
    );

# Prepare CRISPR-Cas data for barplot
effect.05.box.4.FOXP4 <- effect.05.box[effect.05.box$gene %in% i, ];
effect.05.box.4.FOXP4$label <- rep('Cas', nrow(effect.05.box.4.FOXP4));
rownames(effect.05.box.4.FOXP4) <- colnames(effect.score.05.overlap.minus.05);
effect.05.box.4.FOXP4.order <- effect.05.box.4.FOXP4[match(i.protein.na.df$sample, rownames(effect.05.box.4.FOXP4)), ];
effect.05.box.4.FOXP4.order$order <- 1:nrow(effect.05.box.4.FOXP4.order);

# Create CRISPR-Cas gene effect score barplot
bar.cas.na <- BoutrosLab.plotting.general::create.barplot(
    formula = score ~ order,
    main = NULL,
    ylab.label = expression('Gene effect score (CRISPR-Cas)'),
    xlab.label = NULL,
    xaxis.lab = rownames(effect.05.box.4.FOXP4.order),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0, 0),
    xaxis.cex = 0,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.2,
    main.cex = 1.4,
    col = c('red3', rep('grey20', (nrow(i.fpkm.data) - 1))),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    data = effect.05.box.4.FOXP4.order,
    sample.order = 'none',
    add.rectangle = FALSE
    );

# Create legend for CNV plot
legendG <- legend.grob(
    list(
        legend = list(
            colours = c('dodgerblue4', 'white', 'red4'),
            title = expression(underline('Copy number')),
            labels = c('0', max.lim.cnv),
            size = 1,
            label.cex = 1,
            continuous = TRUE,
            height = 2
            )
        ),
    label.cex = 1,
    title.cex = 1,
    title.just = 'left',
    x = 0.2,
    y = 1.3
    );

# Combine all plots into a multi-panel plot
multi.gene.protein.bar <- create.multipanelplot(
    list(bar.fpkm, bar.protein.na, bar.rnai.na, bar.cas.na, cnv.plot),
    main = as.expression(substitute(paste(var), list(var = i))),
    main.cex = 1.5,
    main.y = 0.1,
    resolution = 300,
    layout.height = 5,
    layout.width = 1,
    layout.skip = c(FALSE, FALSE, FALSE, FALSE, FALSE),
    plot.objects.heights = c(9, 8, 8, 8, 3),
    legend = list(bottom = list(fun = legendG)),
    y.spacing = c(-9, -9, -9, -2, -9),
    right.legend.padding = 0
    );

save.outlier.figure(
    multi.gene.protein.bar,
    c('Figure4k', i, 'multi', 'bar'),
    width = 9,
    height = 12
    );

save.session.profile(file.path('output', 'Figure4k.txt'));
