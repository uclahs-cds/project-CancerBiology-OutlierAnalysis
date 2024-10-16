### HISTORY ######################################################################
# This script analyzes the protein expression data (RPPA) for outlier and
# non-outlier genes in the TCGA-BRCA dataset.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script processes and analyzes Reverse Phase Protein Array (RPPA) data
# for outlier and non-outlier genes in breast cancer samples from the TCGA-BRCA
# dataset. It performs the following main tasks:
# 1. Identifies outlier genes with available RPPA data
# 2. Processes RPPA data, excluding phosphorylated proteins
# 3. Compares protein abundance between outlier and non-outlier patients
# 4. Performs statistical analysis (Wilcoxon test) on the differences
# 5. Creates a boxplot visualization of the protein abundance distribution

### PREAMBLE #####################################################################
# Load necessary library
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'outlier.symbol'
    ));

# Protein gene list from antibody data
protein.gene <- unlist(strsplit(protein.antibody$gene_name, '/'));

# Outlier genes with protein data
outlier.protein.gene <- outlier.symbol$brca[outlier.symbol$brca %in% protein.gene];

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

outlier.patient.tag.01.brca.protein.match <- outlier.patient.tag.01.brca[
    rownames(fpkm.tumor.symbol.filter.brca)[fpkm.tumor.symbol.filter.brca$Symbol %in% unique(outlier.protein.gene)],
    colnames(brca.protein.outlier.match)
    ];

outlier.protein.list <- list();
non.outlier.protein.list <- list();
target.gene.list <- NULL;

for (i in 1:nrow(brca.protein.outlier.match)) {
    target.gene.name <- protein.antibody[rownames(brca.protein.outlier.match), 'gene_name'][i];
    target.gene.name.split <- unlist(strsplit(target.gene.name, '/'));
    target.gene.name.single <- outlier.protein.gene[outlier.protein.gene %in% target.gene.name.split];
    row.name.target <- rownames(fpkm.tumor.symbol.filter.brca)[fpkm.tumor.symbol.filter.brca$Symbol %in% target.gene.name.single];
    target.col <- colnames(outlier.patient.tag.01.brca.protein.match)[outlier.patient.tag.01.brca.protein.match[row.name.target, ] == 1];
    non.target.col <- colnames(outlier.patient.tag.01.brca.protein.match)[outlier.patient.tag.01.brca.protein.match[row.name.target, ] == 0];
    target.gene.list <- c(target.gene.list, target.gene.name.single);
    outlier.protein.list[[i]] <- brca.protein.outlier.match[i, target.col];
    non.outlier.protein.list[[i]] <- brca.protein.outlier.match[i, non.target.col];
    }

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

outlier.protein.list.no.p <- list();
non.outlier.protein.list.no.p <- list();
target.gene.list.no.p <- NULL;

for (i in 1:nrow(brca.protein.outlier.match.no.p)) {
    target.gene.name <- protein.antibody[rownames(brca.protein.outlier.match.no.p), 'gene_name'][i];
    target.gene.name.split <- unlist(strsplit(target.gene.name, '/'));
    target.gene.name.single <- outlier.protein.gene[outlier.protein.gene %in% target.gene.name.split];
    row.name.target <- rownames(fpkm.tumor.symbol.filter.brca)[fpkm.tumor.symbol.filter.brca$Symbol %in% target.gene.name.single];
    target.col <- colnames(outlier.patient.tag.01.brca.protein.match.no.p)[outlier.patient.tag.01.brca.protein.match.no.p[row.name.target, ] == 1];
    non.target.col <- colnames(outlier.patient.tag.01.brca.protein.match.no.p)[outlier.patient.tag.01.brca.protein.match.no.p[row.name.target, ] == 0];
    target.gene.list.no.p <- c(target.gene.list.no.p, target.gene.name.single);
    outlier.protein.list.no.p[[i]] <- brca.protein.outlier.match.no.p[i, target.col];
    non.outlier.protein.list.no.p[[i]] <- brca.protein.outlier.match.no.p[i, non.target.col];
    }

# Box plot
# Exclude genes with no outlier patient info

names(outlier.protein.list.no.p) <- rownames(brca.protein.outlier.match.no.p);
outlier.protein.list.no.p.na <- na.omit(unlist(outlier.protein.list.no.p));
names(non.outlier.protein.list.no.p) <- rownames(brca.protein.outlier.match.no.p);
non.outlier.protein.list.no.p.na <- non.outlier.protein.list.no.p[names(outlier.protein.list.no.p.na)];

protein.na.value <- data.frame(
    protein.na.value = c(
        as.numeric(unlist(non.outlier.protein.list.no.p.na)),
        as.numeric(unlist(outlier.protein.list.no.p.na))
        )
    );

protein.na.value.box <- data.frame(
    cbind(
        protein.na.value$protein.na.value,
        c(
            rep('non', length(as.numeric(unlist(non.outlier.protein.list.no.p.na)))),
            rep('out', length(as.numeric(unlist(outlier.protein.list.no.p.na))))
            )
        )
    );

colnames(protein.na.value.box) <- c('protein.value', 'status');
protein.na.value.box[, 1] <- as.numeric(protein.na.value.box[, 1]);

wilcox.result.protein.na <- wilcox.test(
    as.numeric(unlist(outlier.protein.list.no.p.na)),
    as.numeric(unlist(non.outlier.protein.list.no.p.na)),
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
set.seed(sum(utf8ToInt('Figure3b')));

rppa.box <- BoutrosLab.plotting.general::create.boxplot(
    formula = protein.value ~ status,
    data = protein.na.value.box,
    main = expression('Protein abundance of outlier genes'),
    main.cex = 1.3,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Protein abundance (RPPA)'),
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
    ylimits = c(-4, 8),
    add.stripplot = TRUE,
    points.pch = 1,
    points.cex = 0.8,
    points.col = 'grey60',
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 4),
    xright.rectangle = c(4, 5),
    ybottom.rectangle = -4,
    ytop.rectangle = 8,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    lwd = 1.2,
    col = c('red2', 'dodgerblue3'),
    alpha = 0.3
    );


save.outlier.figure(
    rppa.box,
    c('Figure3b', 'rppa', 'box'),
    width = 3.5,
    height = 6.5
    );

save.session.profile(file.path('output', 'Figure3b.txt'));
