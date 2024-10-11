### HISTORY #####################################################################
# This script analyzes the protein abundance of outlier genes identified in CCLE
# overlapped with tissue datasets.
# Date: 2024-08-16

### DESCRIPTION #################################################################
# This script processes data from the Cancer Cell Line Encyclopedia (CCLE) and
# tissue datasets to analyze the protein abundance of outlier genes. It compares
# the protein abundance between outlier and non-outlier samples and visualizes
# the results using a boxplot.

### PREAMBLE ####################################################################
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'five.data.outlier.symbol',
    'ccle.sample.outlier.status.fdr.05'
    ));

ccle.sample.outlier.status.fdr.05.five <- ccle.sample.outlier.status.fdr.05[sub('\\..*', '', rownames(ccle.sample.outlier.status.fdr.05)) %in% five.data.outlier.symbol, ];
ccle.sample.outlier.status.fdr.05.five.symbol <- sub('\\..*', '', rownames(ccle.sample.outlier.status.fdr.05.five));


protein.info.breast.num.symbol <- sapply(strsplit(rownames(protein.info.breast.num), '\\|'), function(x) x[3]);
protein.info.breast.num.symbol <- sub('_HUMAN', '', protein.info.breast.num.symbol);
protein.info.breast.num.match <- protein.info.breast.num[, colnames(protein.info.breast.num) %in% colnames(fpkm.tumor.symbol.filter.ccle)];
rownames(protein.info.breast.num.match) <- protein.info.breast.num.symbol;
protein.info.breast.num.match.05 <- protein.info.breast.num.match[rownames(protein.info.breast.num.match) %in% ccle.sample.outlier.status.fdr.05.five.symbol, ];

cache.multiple.computed.variables(c(
    'protein.info.breast.num.match',
    'protein.info.breast.num.symbol'
    ));

# Use the mean of the duplicated gene and patient
#   - duplicated gene
uni.gene.protein <- unique(rownames(protein.info.breast.num.match.05));
protein.info.breast.num.match.05.uni <- NULL;
for (i in 1:length(uni.gene.protein)) {
    if (1 == sum(rownames(protein.info.breast.num.match.05) %in% uni.gene.protein[i])) {
        dup.gene.protein.num.mean <- protein.info.breast.num.match.05[rownames(protein.info.breast.num.match.05) %in% uni.gene.protein[i], ];
        } else {
        dup.gene.protein.num <- protein.info.breast.num.match.05[rownames(protein.info.breast.num.match.05) %in% uni.gene.protein[i], ];
        dup.gene.protein.num.mean <- apply(dup.gene.protein.num, 2, function(x) {
            mean(na.omit(x))
            });
        }
    protein.info.breast.num.match.05.uni <- rbind(protein.info.breast.num.match.05.uni, dup.gene.protein.num.mean);
    }
rownames(protein.info.breast.num.match.05.uni) <- uni.gene.protein;

#   - duplicated sample
uni.sample.protein <- unique(colnames(protein.info.breast.num.match.05));
protein.info.breast.num.match.05.uni.gene.sample <- NULL;
for (i in 1:length(uni.sample.protein)) {
    if (1 == sum(colnames(protein.info.breast.num.match.05) %in% uni.sample.protein[i])) {
        dup.gene.protein.num.mean <- protein.info.breast.num.match.05.uni[, colnames(protein.info.breast.num.match.05) %in% uni.sample.protein[i]];
        } else {
        dup.gene.protein.num <- protein.info.breast.num.match.05.uni[, colnames(protein.info.breast.num.match.05) %in% uni.sample.protein[i]];
        dup.gene.protein.num.mean <- apply(dup.gene.protein.num, 1, function(x) {
            mean(na.omit(x))
            });
        }
    protein.info.breast.num.match.05.uni.gene.sample <- cbind(protein.info.breast.num.match.05.uni.gene.sample, dup.gene.protein.num.mean);
    }
colnames(protein.info.breast.num.match.05.uni.gene.sample) <- uni.sample.protein;



# only outlier gene's fpkm
fpkm.tumor.symbol.filter.ccle.outlier <- fpkm.tumor.symbol.filter.ccle[sub('\\..*', '', rownames(fpkm.tumor.symbol.filter.ccle)) %in% rownames(protein.info.breast.num.match.05.uni.gene.sample), colnames(protein.info.breast.num.match.05.uni.gene.sample)];

# only outlier gene's status
ccle.sample.outlier.status.protein.match <- ccle.sample.outlier.status[rownames(fpkm.tumor.symbol.filter.ccle.outlier), colnames(protein.info.breast.num.match.05.uni.gene.sample)];

rownames(fpkm.tumor.symbol.filter.ccle.outlier) <- sub('\\..*', '', rownames(fpkm.tumor.symbol.filter.ccle.outlier));
rownames(ccle.sample.outlier.status.protein.match) <- sub('\\..*', '', rownames(ccle.sample.outlier.status.protein.match));


outlier.protein.ccle.zscore.list <- list();
non.outlier.protein.ccle.zscore.list <- list();
target.gene.ccle.zscore.list <- NULL;
for (i in 1:nrow(protein.info.breast.num.match.05.uni.gene.sample)) {
    # protein target name
    target.gene.name.protein <- rownames(protein.info.breast.num.match.05.uni.gene.sample)[i];

    # outlier patient column
    target.col <- colnames(ccle.sample.outlier.status.protein.match)[ccle.sample.outlier.status.protein.match[target.gene.name.protein, ] == 1];
    non.target.col <- colnames(ccle.sample.outlier.status.protein.match)[ccle.sample.outlier.status.protein.match[target.gene.name.protein, ] == 0];

    #

    target.gene.ccle.zscore.list <- c(target.gene.ccle.zscore.list, target.gene.name.protein);
    outlier.protein.ccle.zscore.list[[i]] <- protein.info.breast.num.match.05.uni.gene.sample[i, target.col];
    non.outlier.protein.ccle.zscore.list[[i]] <- protein.info.breast.num.match.05.uni.gene.sample[i, non.target.col];
    }


# box plot - compare the values between patients
# - excluding the genes having no outlier patient info

names(outlier.protein.ccle.zscore.list) <- target.gene.ccle.zscore.list;
outlier.protein.ccle.list.no.p.na <- na.omit(unlist(outlier.protein.ccle.zscore.list))
names(non.outlier.protein.ccle.zscore.list) <- target.gene.ccle.zscore.list;
non.outlier.protein.ccle.list.no.p.na <- non.outlier.protein.ccle.zscore.list[names(outlier.protein.ccle.list.no.p.na)];

protein.ccle.na.value <- data.frame(protein.ccle.na.value = c(as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na)), as.numeric(unlist(outlier.protein.ccle.list.no.p.na))));
protein.ccle.na.value.box <- data.frame(cbind(protein.ccle.na.value$protein.ccle.na.value, c(rep('non', length(as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na)))), rep('out', length(as.numeric(unlist(outlier.protein.ccle.list.no.p.na)))))));
colnames(protein.ccle.na.value.box) <- c('protein.ccle.value', 'status');
protein.ccle.na.value.box[, 1] <- as.numeric(protein.ccle.na.value.box[, 1]);

wilcox.result.protien.na <- wilcox.test(as.numeric(unlist(outlier.protein.ccle.list.no.p.na)),
    as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na)),
    alternative = 'two.sided', conf.int = TRUE
    );


text.pvalue.protien.na <- display.statistical.result(
    x = wilcox.result.protien.na$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

key.protien.na <- list(
    text = list(
        lab = text.pvalue.protien.na,
        cex = 1
        ),
    x = 0.25,
    y = 0.95
    );

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure4c')));

ccle.protein.box <- BoutrosLab.plotting.general::create.boxplot(
    formula = protein.ccle.value ~ status,
    data = protein.ccle.na.value.box,
    main = expression('Protein abundance of outlier genes'),
    main.cex = 1.3,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Protein abundance (z-score)'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.lab = c('Non-outlier\n samples', 'Outlier\n samples'),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.rot = 0,
    outliers = FALSE,
    key = key.protien.na,
    ylimits = c(-5, 8.5),
    # yat = seq(-110, 110, 20),
    # add.text = TRUE,
    # text.x = 1.5,
    # text.y = 6.2,
    # text.fontface = 1,
    # text.labels = paste('p =', sprintf("%.1e",wilcox.result.protien$p.value)),
    add.stripplot = TRUE,
    points.pch = 1,
    points.cex = 0.8,
    points.col = 'grey60',
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 4),
    xright.rectangle = c(4, 5),
    ybottom.rectangle = -6,
    ytop.rectangle = 10,
    # set rectangle colour
    col.rectangle = 'grey',
    # set rectangle alpha (transparency)
    alpha.rectangle = 0.25,
    lwd = 1.2,
    col = c('red2', 'dodgerblue3'),
    alpha = 0.3
    );


save.outlier.figure(
    ccle.protein.box,
    c('Figure4c', 'CCLE', 'outlier', 'protein', 'box'),
    width = 3.5,
    height = 6.5
    );

save.session.profile(file.path('output', 'Figure4c.txt'));
