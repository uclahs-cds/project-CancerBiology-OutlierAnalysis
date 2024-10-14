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
    c('Figure3acd', 'cptac', 'box'),
    width = 3.5,
    height = 6.5
    );

# Calculate quantiles of protein abundance for outlier genes
percent.protein.cptac.quantile <- NULL;
for (i in 1:length(outlier.protein.cptac.list.no.p.na)) {
    unequal.quan <- rev(seq(0, 0.9, 0.1));
    value.vector <- na.omit(as.numeric(unlist(non.outlier.protein.cptac.list.no.p.na[i])));
    non.value <- quantile(value.vector, p = unequal.quan);
    out.value <- as.numeric(unlist(outlier.protein.cptac.list.no.p.na[i]));
    all.value <- c(out.value, non.value);
    percent.protein.cptac.quantile <- rbind(percent.protein.cptac.quantile, all.value);
    }

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
    c('Figure3acd', 'cptac', 'heatmap'),
    width = 6,
    height = 4.5
    );

### 1. CPTAC example
# Two inputs: IGF1R, CLU - run separately
do.plot.3d <- function(i) {
    # Extract matching protein data for a specific gene (i)
    brca.protein.cptac.outlier.match <- brca.protein.cptac[, colnames(brca.protein.cptac) %in% substr(colnames(outlier.patient.tag.01.brca), 1, 15)];
    brca.protein.cptac.outlier.match.i <- brca.protein.cptac.outlier.match[i, ];

    # Match corresponding FPKM data
    fpkm.protein.cptac.match.i <- fpkm.tumor.symbol.filter.brca[
        fpkm.tumor.symbol.filter.brca$Symbol == i,
        match(colnames(brca.protein.cptac.outlier.match.i), substr(colnames(fpkm.tumor.symbol.filter.brca), 1, 15))
        ];

    # Identify outlier patients for the specific gene (i)
    outlier.patient.tag.01.brca.protein.cptac.zscore.match.i <- colnames(
        outlier.patient.tag.01.brca.protein.cptac.zscore.match
        )[outlier.patient.tag.01.brca.protein.cptac.zscore.match[
        rownames(fpkm.protein.cptac.match.i),
        ] == 1];

    # Extract FPKM and protein CPTAC data for outlier patients
    i.fpkm <- fpkm.protein.cptac.match.i[outlier.patient.tag.01.brca.protein.cptac.zscore.match.i];
    i.protein.cptac <- brca.protein.cptac.outlier.match.i[
        substr(outlier.patient.tag.01.brca.protein.cptac.zscore.match.i, 1, 15)
        ];

    # Combine protein and mRNA data into a comparison data frame
    protein.cptac.rna.i.comparison <- data.frame(
        cbind(
            as.numeric(brca.protein.cptac.outlier.match.i[1, ]),
            as.numeric(fpkm.protein.cptac.match.i)
            )
        );

    # Rename columns for clarity
    colnames(protein.cptac.rna.i.comparison) <- c('non', 'out');

    # Define colors for scatter plot based on outlier status
    dot.colours <- rep('black', nrow(protein.cptac.rna.i.comparison));
    dot.colours[
        which(colnames(fpkm.protein.cptac.match.i) == outlier.patient.tag.01.brca.protein.cptac.zscore.match.i)
        ] <- 'red2';

    cptac.gene.scatter <- create.scatterplot(
        formula = non ~ log2(out),
        data = protein.cptac.rna.i.comparison,
        col = dot.colours,
        alpha = .55,
        xaxis.fontface = 1,
        yaxis.fontface = 1,
        yaxis.tck = c(0.2, 0),
        xaxis.tck = c(0.2, 0),
        add.grid = TRUE,
        grid.colour = 'grey80',
        cex = 0.75,
        xaxis.cex = 1,
        yaxis.cex = 1,
        xlab.cex = 1.3,
        ylab.cex = 1.3,
        main.cex = 1.5,
        main = i,
        xlab.label = expression('mRNA abundance (log'[2] * '(FPKM))'),
        ylab.label = expression('Protein abundance (z-score)'),
        text.x = log2(as.numeric(i.fpkm)),
        text.y = as.numeric(i.protein.cptac),
        text.labels = '*outlier patient',
        text.guess.labels = TRUE,
        text.guess.label.position = 180,
        text.guess.radius.factor = 1.5,
        text.fontface = 1,
        text.col = 'red2',
        add.text = TRUE,
        type = c('p', 'r', 'g'),
        legend = list(
            inside = list(
                fun = draw.key,
                args = list(
                    key = get.corr.key(
                        x = protein.cptac.rna.i.comparison$non,
                        y = protein.cptac.rna.i.comparison$out,
                        label.items = c('spearman'),
                        alpha.background = 0,
                        key.cex = 1.1
                        )
                    ),
                x = 0.03,
                y = 0.95,
                corner = c(0, 1)
                )
            ),
        );


    save.outlier.figure(
        cptac.gene.scatter,
        c('Figure3acd', i, 'scatter', 'cptac'),
        width = 6,
        height = 6
        );
    }


do.plot.3d('IGF1R');
do.plot.3d('CLU');

### 2. RPPA example
outlier.symbol <- fpkm.tumor.symbol.filter.brca[rownames(outlier.patient.tag.01.brca), 'Symbol'];
protein.gene <- unlist(strsplit(protein.antibody$gene_name, '/'));

# Outlier genes with protein data
outlier.protein.gene <- outlier.symbol[outlier.symbol %in% protein.gene];
protein.antibody.outlier <- NULL;

for (i in 1:nrow(protein.antibody)) {
    if (sum(unlist(strsplit(protein.antibody$gene_name[i], '/')) %in% outlier.protein.gene) > 0) {
        protein.antibody.outlier <- rbind(protein.antibody.outlier, protein.antibody[i, ]);
        }
    }

protein.antibody.outlier.id <- rownames(protein.antibody.outlier);
brca.protein.outlier <- brca.protein[protein.antibody.outlier.id, 5:ncol(brca.protein)];
brca.protein.outlier.match <- brca.protein.outlier[
    ,
    colnames(brca.protein.outlier) %in% colnames(outlier.patient.tag.01.brca)
    ];


# Exclude phosphorylated protein
protein.antibody.outlier.no.p <- protein.antibody.outlier[
    -(grep('_p', protein.antibody.outlier$peptide_target)),
    ];

protein.antibody.outlier.id.no.p <- rownames(protein.antibody.outlier.no.p);
brca.protein.outlier.no.p <- brca.protein[protein.antibody.outlier.id.no.p, 5:ncol(brca.protein)];
brca.protein.outlier.match.no.p <- brca.protein.outlier.no.p[
    ,
    colnames(brca.protein.outlier.no.p) %in% colnames(outlier.patient.tag.01.brca)
    ];

outlier.patient.tag.01.brca.protein.match.no.p <- outlier.patient.tag.01.brca[
    rownames(fpkm.tumor.symbol.filter.brca)[fpkm.tumor.symbol.filter.brca$Symbol %in% unique(protein.antibody.outlier.no.p$gene_name)],
    colnames(brca.protein.outlier.match.no.p)
    ];

i <- 'NRAS';
# Extract i protein data from the matched dataset
brca.protein.outlier.match.i <- brca.protein.outlier.match[
    rownames(protein.antibody.outlier)[protein.antibody.outlier$gene_name == i],
    ];

# Match corresponding FPKM data for i
fpkm.protein.match.i <- fpkm.tumor.symbol.filter.brca[
    fpkm.tumor.symbol.filter.brca$Symbol == i,
    colnames(brca.protein.outlier.match.i)
    ];

# Combine protein and mRNA data into a comparison data frame for i
protein.rna.i.comparison <- data.frame(
    cbind(
        as.numeric(brca.protein.outlier.match.i[1, ]),
        as.numeric(fpkm.protein.match.i)
        )
    );

# Rename columns for clarity
colnames(protein.rna.i.comparison) <- c('non', 'out');

# Identify outlier patients for i
outlier.patient.tag.01.brca.protein.match.no.p.i <- colnames(
    outlier.patient.tag.01.brca.protein.match.no.p
    )[outlier.patient.tag.01.brca.protein.match.no.p[
    rownames(fpkm.protein.match.i),
    ] == 1];

# Extract FPKM and protein data for outlier patients
i.fpkm <- fpkm.protein.match.i[
    outlier.patient.tag.01.brca.protein.match.no.p.i
    ];

i.protein <- brca.protein.outlier.match.i[
    1, outlier.patient.tag.01.brca.protein.match.no.p.i
    ];

# Define colors for scatter plot based on outlier status
dot.colours <- rep('black', nrow(protein.rna.i.comparison));
dot.colours[
    which(colnames(fpkm.protein.match.i) %in% outlier.patient.tag.01.brca.protein.match.no.p.i)
    ] <- 'red2';


rppa.gene.scatter <- BoutrosLab.plotting.general::create.scatterplot(
    formula = non ~ log2(out),
    data = protein.rna.i.comparison[nrow(protein.rna.i.comparison):1, ],
    col = rev(dot.colours),
    alpha = .55,
    grid.colour = 'grey80',
    cex = 0.75,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.5,
    main = i,
    xlab.label = expression('mRNA abundance (log'[2] * '(FPKM))'),
    ylab.label = expression('Protein abundance (RPPA)'),
    text.x = log2(as.numeric(i.fpkm)),
    text.y = as.numeric(i.protein),
    text.labels = rep('*outlier patient', length(i.protein)),
    text.guess.labels = TRUE,
    text.guess.label.position = 160,
    text.guess.radius.factor = 1.5,
    text.fontface = 1,
    text.col = 'red2',
    add.text = TRUE,
    type = c('p', 'r', 'g'),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = protein.rna.i.comparison$non,
                    y = protein.rna.i.comparison$out,
                    label.items = c('spearman'),
                    alpha.background = 0,
                    key.cex = 1.1
                    )
                ),
            x = 0.03,
            y = 0.95,
            corner = c(0, 1)
            )
        )
    );

save.outlier.figure(
    rppa.gene.scatter,
    c('Figure3acd', i, 'scatter', 'rppa'),
    width = 6,
    height = 6
    );

save.session.profile(file.path('output', 'Figure3acd.txt'));
