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
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'outlier.symbol',
    'outlier.gene.fdr.01'
    ));

# Protein CPTAC z-score list
protein.cptac.zscore.gene <- rownames(brca.protein.cptac.zscore);

# Outlier genes with protein CPTAC z-score data
outlier.protein.cptac.zscore.gene <- outlier.symbol$brca[outlier.symbol$brca %in% protein.cptac.zscore.gene];

brca.protein.cptac.zscore.outlier.match <- brca.protein.cptac.zscore[
    ,
    colnames(brca.protein.cptac.zscore) %in% substr(colnames(outlier.patient.tag.01.brca), 1, 15)
    ];

outlier.patient.tag.01.brca.protein.cptac.zscore.match <- outlier.patient.tag.01.brca[
    rownames(fpkm.tumor.symbol.filter.brca)[
        fpkm.tumor.symbol.filter.brca$Symbol %in% unique(outlier.protein.cptac.zscore.gene)
        ],
    substr(colnames(outlier.patient.tag.01.brca), 1, 15) %in% colnames(brca.protein.cptac.zscore)
    ];

# Only outlier gene's FPKM
fpkm.tumor.symbol.filter.brca.outlier <- fpkm.tumor.symbol.filter.brca[
    rownames(outlier.gene.fdr.01$brca),
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
    # target.col <- colnames(outlier.patient.tag.01.brca.protein.cptac.zscore.match)[
    #     outlier.patient.tag.01.brca.protein.cptac.zscore.match[row.name.target, ][1,] == 1
    #     ];
    # non.target.col <- colnames(outlier.patient.tag.01.brca.protein.cptac.zscore.match)[
    #     outlier.patient.tag.01.brca.protein.cptac.zscore.match[row.name.target, ][1,] == 0
    #     ];
    target.col <- colnames(outlier.patient.tag.01.brca.protein.cptac.zscore.match)[
        apply(outlier.patient.tag.01.brca.protein.cptac.zscore.match[row.name.target, ], 2, sum) == 1
        ];
    non.target.col <- colnames(outlier.patient.tag.01.brca.protein.cptac.zscore.match)[
        apply(outlier.patient.tag.01.brca.protein.cptac.zscore.match[row.name.target, ], 2, sum) == 0
        ];
    target.gene.cptac.zscore.list <- c(target.gene.cptac.zscore.list, outlier.protein.cptac.zscore.gene[i]);
    outlier.protein.cptac.zscore.list[[i]] <- target.gene.name.protein[, substr(target.col, 1, 15)];
    non.outlier.protein.cptac.zscore.list[[i]] <- target.gene.name.protein[, substr(non.target.col, 1, 15)];
    }

# Box plot - compare the values between patients
# Exclude the genes with no outlier patient info

names(outlier.protein.cptac.zscore.list) <- outlier.protein.cptac.zscore.gene;
names(non.outlier.protein.cptac.zscore.list) <- outlier.protein.cptac.zscore.gene;


outlier.protein.cptac.list.no.p.na <- na.omit(unlist(outlier.protein.cptac.zscore.list))
non.outlier.protein.cptac.list.no.p.na <- non.outlier.protein.cptac.zscore.list[names(outlier.protein.cptac.list.no.p.na)];

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
    c('Figure1i', 'cptac', 'box'),
    width = 3.5,
    height = 6.5
    );



# Calculate quantiles of protein abundance for outlier genes
percent.protein.cptac.quantile <- NULL;
for (i in 1:length(outlier.protein.cptac.zscore.list)) {
    unequal.quan <- rev(seq(0, 0.9, 0.1));
    value.vector <- na.omit(as.numeric(unlist(non.outlier.protein.cptac.zscore.list[i])));
    non.value <- quantile(value.vector, p = unequal.quan);
    out.value <- mean(as.numeric(unlist(outlier.protein.cptac.zscore.list[i])));
    all.value <- c(out.value, non.value);
    percent.protein.cptac.quantile <- rbind(percent.protein.cptac.quantile, all.value);
    }
percent.protein.cptac.quantile <- na.omit(percent.protein.cptac.quantile);

# Prepare data for heatmap
heat.df <- t(data.frame(percent.protein.cptac.quantile));
rownames(heat.df) <- c('Outliers', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100');
heat.df.rev <- heat.df[rev(seq(nrow(heat.df))), ];

# Define legend for the heatmap
legend.col <- list(
    legend = list(
        colours = c('black', 'white', '#b2182b'),
        title = expression(underline('z-score')),
        labels = c(-3, 0, 3),
        size = 3,
        label.cex = 1,
        continuous = TRUE,
        height = 3
        )
    );

# Generate the heatmap
heat.out <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(heat.df.rev),
    clustering.method = 'none',
    colour.scheme = c('black', 'white', '#b2182b'),
    col.colour = 'white',
    grid.row = FALSE,
    grid.col = TRUE,
    yaxis.tck = 0,
    xaxis.tck = 0,
    xaxis.lab = NULL,
    yaxis.lab = rev(c('Outliers', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100')),
    xlab.label = expression('Outlier Genes'),
    yaxis.cex = 1.2,
    xaxis.cex = 1.2,
    yaxis.rot = 0,
    xaxis.rot = 90,
    xlab.cex = 1.2,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    colour.centering.value = 0,
    at = seq(-3, 3, 0.01),
    covariate.legend = legend.col,
    legend.cex = 1,
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );


save.outlier.figure(
    heat.out,
    c('Figure1k', 'cptac', 'heatmap'),
    width = 6,
    height = 4.5
    );



### Example gene - CPTAC
# - CLU gene

# Use rank for protein
brca.protein.cptac.outlier.match.CLU <- brca.protein.cptac.zscore.outlier.match['CLU',];
fpkm.protein.cptac.match.CLU.row <- fpkm.tumor.symbol.filter.brca[fpkm.tumor.symbol.filter.brca$Symbol == 'CLU', match( colnames(brca.protein.cptac.outlier.match.CLU),substr(colnames(fpkm.tumor.symbol.filter.brca), 1, 15))];

fpkm.protein.cptac.match.CLU <- scale(as.numeric(fpkm.tumor.symbol.filter.brca[fpkm.tumor.symbol.filter.brca$Symbol == 'CLU', patient.part.brca]))[match(colnames(brca.protein.cptac.outlier.match.CLU),substr(colnames(fpkm.tumor.symbol.filter.brca), 1, 15)),]

outlier.patient.tag.01.brca.protein.cptac.zscore.match.CLU <- colnames(outlier.patient.tag.01.brca.protein.cptac.zscore.match)[outlier.patient.tag.01.brca.protein.cptac.zscore.match[rownames(fpkm.protein.cptac.match.CLU.row),] == 1];
CLU.fpkm <- fpkm.protein.cptac.match.CLU[outlier.patient.tag.01.brca.protein.cptac.zscore.match.CLU];
CLU.protein.cptac <- brca.protein.cptac.outlier.match.CLU[substr(outlier.patient.tag.01.brca.protein.cptac.zscore.match.CLU, 1, 15)];

protein.cptac.rna.CLU.comparison <- data.frame(cbind(rank(-(as.numeric(brca.protein.cptac.outlier.match.CLU[1,]))),
                                         as.numeric(fpkm.protein.cptac.match.CLU)));
colnames(protein.cptac.rna.CLU.comparison) <- c('non', 'out');

dot.colours <- vector(length= nrow(protein.cptac.rna.CLU.comparison));
dot.colours <- rep('black',nrow(protein.cptac.rna.CLU.comparison));
dot.colours[which(colnames(fpkm.protein.cptac.match.CLU.row) == outlier.patient.tag.01.brca.protein.cptac.zscore.match.CLU)] <- 'red2';


CLU.scatter <- create.scatterplot(
    formula = -non ~ out,
    data = protein.cptac.rna.CLU.comparison,
    col = dot.colours,
    alpha = .75,
    yat = c(-100, -50, -1),
    ylimits = c(-110, 15),
    yaxis.lab = c(100, 50, 1),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    add.grid = TRUE,
	grid.colour = 'grey80',
    cex = 1,
    xaxis.cex = 1.3,
    yaxis.cex = 1.3,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.3,
    main = expression('CLU'),
    xlab.label = expression('mRNA abundance (z-score)'),
    ylab.label = expression('Protein abundance (rank)'),
    text.x = log2(as.numeric(CLU.fpkm)),
    text.y = as.numeric(CLU.protein.cptac),
    text.labels = '*outlier patient',
    text.guess.labels = TRUE,
    text.guess.label.position = 180,
    text.guess.radius.factor = 1.5,
    text.fontface = 1,
    text.col = 'red2',
    add.text = TRUE,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = -protein.cptac.rna.CLU.comparison$non,
                    y = protein.cptac.rna.CLU.comparison$out,
                    label.items = c('spearman'),
                    alpha.background = 0,
                    key.cex = 1.1
                    )
                ),
            x = 0.03,
            y = 0.95,
            corner = c(0,1)
            )
        ),
    );



save.outlier.figure(
    CLU.scatter,
    c('Figure1l', 'CPTAC_CLU', 'box'),
    width = 3.5,
    height = 6.5
    );




save.session.profile(file.path('output', 'Figure1ikl_9.txt'));
