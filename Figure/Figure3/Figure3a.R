### HISTORY ######################################################################
# This script analyzes protein abundance data (z-scores) for outlier and
# non-outlier genes in the TCGA-BRCA dataset using CPTAC data.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script processes and analyzes protein abundance data for outlier and 
# non-outlier genes in breast cancer samples from the TCGA-BRCA dataset, using 
# CPTAC (Clinical Proteomic Tumor Analysis Consortium) data. It performs the 
# following main tasks:
# 1. Identifies outlier genes with available protein data
# 2. Compares protein abundance between outlier and non-outlier patients
# 3. Performs statistical analysis (Wilcoxon test) on the differences
# 4. Creates a boxplot visualization of the protein abundance distribution

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
args <- commandArgs();
source(file.path(
    dirname(dirname(normalizePath(sub('^--file=', '', args[grep('^--file=', args)])))),
    'common_functions.R'
    ));
# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-09-10_Figure3a-d.rda'));

# Outlier symbol
outlier.symbol <- fpkm.tumor.symbol.filter.brca[rownames(brca.outlier.patient.tag.01.t.p.order), 'Symbol'];

# Protein CPTAC z-score list
protein.cptac.zscore.gene <- rownames(brca.protein.cptac.zscore);

# Outlier genes with protein CPTAC z-score data
outlier.protein.cptac.zscore.gene <- outlier.symbol[outlier.symbol %in% protein.cptac.zscore.gene];

brca.protein.cptac.zscore.outlier.match <- brca.protein.cptac.zscore[
    ,
    colnames(brca.protein.cptac.zscore) %in% substr(colnames(outlier.patient.tag.01.brca), 1, 15)
    ];

# Only outlier gene's FPKM
fpkm.tumor.symbol.filter.brca.outlier <- fpkm.tumor.symbol.filter.brca[
    rownames(outlier.gene.fdr.01.brca),
    ];

outlier.protein.cptac.zscore.list <- list();
non.outlier.protein.cptac.zscore.list <- list();
target.gene.cptac.zscore.list <- NULL;

for (i in 1:length(outlier.protein.cptac.zscore.gene)) {
    target.gene.name.protein <- brca.protein.cptac.zscore.outlier.match[
        rownames(brca.protein.cptac.zscore.outlier.match) %in% outlier.protein.cptac.zscore.gene[i],
        ];
    row.name.target <- rownames(fpkm.tumor.symbol.filter.brca.outlier)[
        fpkm.tumor.symbol.filter.brca.outlier$Symbol %in% outlier.protein.cptac.zscore.gene[i]
        ];
    target.col <- colnames(outlier.patient.tag.01.brca.protein.cptac.zscore.match)[
        outlier.patient.tag.01.brca.protein.cptac.zscore.match[row.name.target, ] == 1
        ];
    non.target.col <- colnames(outlier.patient.tag.01.brca.protein.cptac.zscore.match)[
        outlier.patient.tag.01.brca.protein.cptac.zscore.match[row.name.target, ] == 0
        ];
    target.gene.cptac.zscore.list <- c(target.gene.cptac.zscore.list, outlier.protein.cptac.zscore.gene[i]);
    outlier.protein.cptac.zscore.list[[i]] <- target.gene.name.protein[, substr(target.col, 1, 15)];
    non.outlier.protein.cptac.zscore.list[[i]] <- target.gene.name.protein[, substr(non.target.col, 1, 15)];
    }

# Box plot - compare the values between patients
# Exclude the genes with no outlier patient info

names(outlier.protein.cptac.zscore.list) <- outlier.protein.cptac.zscore.gene;
names(non.outlier.protein.cptac.zscore.list) <- outlier.protein.cptac.zscore.gene;

protein.cptac.na.value <- data.frame(
    protein.cptac.na.value = c(
        as.numeric(unlist(non.outlier.protein.cptac.list.no.p.na)),
        as.numeric(unlist(outlier.protein.cptac.list.no.p.na))
        )
    );

protein.cptac.na.value.box <- data.frame(
    cbind(
        protein.cptac.na.value$protein.cptac.na.value,
        c(
            rep('non', length(as.numeric(unlist(non.outlier.protein.cptac.list.no.p.na)))),
            rep('out', length(as.numeric(unlist(outlier.protein.cptac.list.no.p.na))))
            )
        )
    );

colnames(protein.cptac.na.value.box) <- c('protein.cptac.value', 'status');
protein.cptac.na.value.box[, 1] <- as.numeric(protein.cptac.na.value.box[, 1]);

wilcox.result.protein.na <- wilcox.test(
    as.numeric(unlist(outlier.protein.cptac.list.no.p.na)),
    as.numeric(unlist(non.outlier.protein.cptac.list.no.p.na)),
    alternative = 'two.sided',
    conf.int = TRUE
    );

text.pvalue.protein.na <- display.statistical.result(
    x = wilcox.result.protein.na$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

key.protein.na <- list(
    text = list(
        lab = text.pvalue.protein.na,
        cex = 1
        ),
    x = 0.25,
    y = 0.95
    );

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure3a')));

cptac.box <- BoutrosLab.plotting.general::create.boxplot(
    formula = protein.cptac.value ~ status,
    data = protein.cptac.na.value.box,
    main = expression('Protein abundance of outlier genes'),
    main.cex = 1.3,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Protein abundance (z-score)'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.lab = c('Non-outlier\n patients', 'Outlier\n patients'),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.rot = 90,
    outliers = FALSE,
    key = key.protein.na,
    ylimits = c(-5, 8.5),
    add.stripplot = TRUE,
    points.pch = 1,
    points.cex = 0.8,
    points.col = 'grey60',
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 4),
    xright.rectangle = c(4, 5),
    ybottom.rectangle = -6,
    ytop.rectangle = 10,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    lwd = 1.2,
    col = c('red2', 'dodgerblue3'),
    alpha = 0.3
    );


save.outlier.figure(
    cptac.box,
    c('Figure3a', 'cptac', 'box'),
    width = 3.5,
    height = 6.5
    );

save.session.profile(file.path('output', 'Figure3a.txt'));
