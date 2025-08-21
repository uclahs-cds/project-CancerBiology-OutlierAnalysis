### HISTORY ######################################################################
# This script analyzes the protein expression data (RPPA) for outlier and
# non-outlier genes in the TCGA-BRCA and I-SPY2 datasets.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script processes and analyzes Reverse Phase Protein Array (RPPA) data
# for outlier and non-outlier genes in breast cancer samples from the TCGA-BRCA
# and I-SPY2 datasets. It performs the following main tasks:
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
    # target.col <- colnames(outlier.patient.tag.01.brca.protein.match)[outlier.patient.tag.01.brca.protein.match[row.name.target, ][1,] == 1];
    # non.target.col <- colnames(outlier.patient.tag.01.brca.protein.match)[outlier.patient.tag.01.brca.protein.match[row.name.target, ][1,] == 0];
    target.col <- colnames(outlier.patient.tag.01.brca.protein.match)[apply(outlier.patient.tag.01.brca.protein.match[row.name.target, ], 2, sum) == 1];
    non.target.col <- colnames(outlier.patient.tag.01.brca.protein.match)[apply(outlier.patient.tag.01.brca.protein.match[row.name.target, ], 2, sum) == 0];
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

outlier.patient.tag.01.protein.match.no.p.brca <- outlier.patient.tag.01.brca[
    rownames(fpkm.tumor.symbol.filter.brca)[fpkm.tumor.symbol.filter.brca$Symbol %in% unique(protein.antibody.outlier.no.p$gene_name)],
    colnames(brca.protein.outlier.match.no.p)
    ];


# scale the data
brca.protein.outlier.match.no.p.scale <- apply(brca.protein.outlier.match.no.p, 1, function(x) {scale(as.numeric(x))});
brca.protein.outlier.match.no.p.scale <- data.frame(t(brca.protein.outlier.match.no.p.scale));
colnames(brca.protein.outlier.match.no.p.scale) <- colnames(brca.protein.outlier.match.no.p);
rownames(brca.protein.outlier.match.no.p.scale) <- protein.antibody[rownames(brca.protein.outlier.match.no.p.scale), 'gene_name'];

outlier.patient.tag.01.protein.match.no.p.name.brca <- outlier.patient.tag.01.protein.match.no.p.brca;
rownames(outlier.patient.tag.01.protein.match.no.p.name.brca) <- fpkm.tumor.symbol.filter.brca[rownames(outlier.patient.tag.01.protein.match.no.p.brca),]$Symbol;




# 2. I-SPY2
ispy.protein.outlier.match <- ispy.protein.outlier[,colnames(ispy.protein.outlier) %in% substr(colnames(outlier.patient.tag.01.ispy), 1, 7)];

ispy.protein.outlier.match.total <- ispy.protein.outlier.match[which(substr(rownames(ispy.protein.outlier.match), nchar(rownames(ispy.protein.outlier.match))-4, nchar(rownames(ispy.protein.outlier.match))) == 'total'),];
rownames(ispy.protein.outlier.match.total) <- substr(rownames(ispy.protein.outlier.match.total), 1, nchar(rownames(ispy.protein.outlier.match.total))-6);

outlier.patient.tag.01.protein.match.ispy.total <- outlier.patient.tag.01.ispy[rownames(ispy.protein.outlier.match.total), substr(colnames(outlier.patient.tag.01.ispy), 1, 7) %in% colnames(ispy.protein.outlier)[1:(ncol(ispy.protein.outlier)-1)]];

# scale data
ispy.protein.outlier.match.total.scale <- apply(ispy.protein.outlier.match.total, 1, function(x) {scale(as.numeric(x))});
ispy.protein.outlier.match.total.scale <- data.frame(t(ispy.protein.outlier.match.total.scale));
colnames(ispy.protein.outlier.match.total.scale) <- colnames(ispy.protein.outlier.match);




# Merge
rppa.tcga.ispy.unique <- unique(c(rownames(brca.protein.outlier.match.no.p.scale), rownames(ispy.protein.outlier.match.total)));
brca.protein.outlier.match.no.p.scale.merge <- brca.protein.outlier.match.no.p.scale[rppa.tcga.ispy.unique,];
outlier.patient.tag.01.protein.match.no.p.name.brca.merge <- outlier.patient.tag.01.protein.match.no.p.name.brca[rppa.tcga.ispy.unique,];
ispy.protein.outlier.match.total.scale.merge <- ispy.protein.outlier.match.total.scale[rppa.tcga.ispy.unique,];
outlier.patient.tag.01.protein.match.ispy.total.merge <- outlier.patient.tag.01.protein.match.ispy.total[rppa.tcga.ispy.unique,];

rownames(brca.protein.outlier.match.no.p.scale.merge) <- rppa.tcga.ispy.unique;
rownames(outlier.patient.tag.01.protein.match.no.p.name.brca.merge) <- rppa.tcga.ispy.unique;
rownames(ispy.protein.outlier.match.total.scale.merge) <- rppa.tcga.ispy.unique;
rownames(outlier.patient.tag.01.protein.match.ispy.total.merge) <- rppa.tcga.ispy.unique;

merge.rppa.data <- cbind(brca.protein.outlier.match.no.p.scale.merge, ispy.protein.outlier.match.total.scale.merge);
merge.rppa.patient <- cbind(outlier.patient.tag.01.protein.match.no.p.name.brca.merge, outlier.patient.tag.01.protein.match.ispy.total.merge);



outlier.protein.list.no.p.scale.merge <- list();
non.outlier.protein.list.no.p.scale.merge <- list();
for (i in 1:nrow(merge.rppa.data)) {
    
    target.col <- colnames(merge.rppa.patient)[merge.rppa.patient[i,] %in% 1];
    non.target.col <- colnames(merge.rppa.patient)[!(merge.rppa.patient[i,] %in% 1)];
    outlier.protein.list.no.p.scale.merge[[i]] <- merge.rppa.data[i,target.col];
    non.outlier.protein.list.no.p.scale.merge[[i]] <- merge.rppa.data[i,non.target.col]
    }


names(outlier.protein.list.no.p.scale.merge) <- rownames(merge.rppa.patient);
outlier.protein.value.no.p.scale.merge <- na.omit(unlist(outlier.protein.list.no.p.scale.merge));
names(non.outlier.protein.list.no.p.scale.merge) <- rownames(merge.rppa.patient);
non.outlier.protein.value.no.p.scale.merge <- non.outlier.protein.list.no.p.scale.merge[names(non.outlier.protein.list.no.p.scale.merge) %in% unique(sub("\\..*", "", names(outlier.protein.value.no.p.scale.merge)))];

outlier.protein.list.no.p.scale.merge.value <- na.omit(as.numeric(unlist(outlier.protein.list.no.p.scale.merge)));
non.outlier.protein.list.no.p.scale.merge.value <- na.omit(as.numeric(unlist(non.outlier.protein.list.no.p.scale.merge)));


protein.na.value <- data.frame(
    protein.na.value = c(
        non.outlier.protein.list.no.p.scale.merge.value,
        outlier.protein.list.no.p.scale.merge.value
        )
    );

protein.na.value.box <- data.frame(
    cbind(
        protein.na.value$protein.na.value,
        c(
            rep('non', length(non.outlier.protein.list.no.p.scale.merge.value)),
            rep('out', length(outlier.protein.list.no.p.scale.merge.value))
            )
        )
    );

colnames(protein.na.value.box) <- c('protein.value', 'status');
protein.na.value.box[, 1] <- as.numeric(protein.na.value.box[, 1]);

wilcox.result.protein.na <- wilcox.test(
    non.outlier.protein.list.no.p.scale.merge.value,
    outlier.protein.list.no.p.scale.merge.value,
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
    ylimits = c(-9, 16.5),
    add.stripplot = TRUE,
    points.pch = 1,
    points.cex = 0.8,
    points.col = 'grey60',
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 4),
    xright.rectangle = c(4, 5),
    ybottom.rectangle = -10,
    ytop.rectangle = 17,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    lwd = 1.2,
    col = c('red2', 'dodgerblue3'),
    alpha = 0.3
    );


save.outlier.figure(
    rppa.box,
    c('Figure1j', 'rppa', 'box'),
    width = 3.5,
    height = 6.5
    );

save.session.profile(file.path('output', 'Figure1j.txt'));
