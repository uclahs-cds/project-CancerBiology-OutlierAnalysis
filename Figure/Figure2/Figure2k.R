### HISTORY #####################################################################
# This script processes gene expression data to compare outlier and non-outlier
# patients across tumor and normal samples. It generates a boxplot to visualize
# the mean beta values and performs a Kruskal-Wallis test for statistical analysis.
# Date: 2024-08-14

# Load necessary library
library(BoutrosLab.plotting.general);

# Source the helper library
args <- commandArgs();
source(file.path(
    dirname(dirname(normalizePath(sub('^--file=', '', args[grep('^--file=', args)])))),
    'common_functions.R'
    ));
# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-08-26_Figure2h-l_input.rda'));


# Main analysis
# divide into outlier and non-outlier


two.out.non.tumor.normal.gene.500 <- mean.minus.ma.merge.two.500$Symbol[mean.minus.ma.merge.two.500$Symbol %in% rownames(normal.tumor.beta.comparison.two.minus.order.500)];


two.out.non.tumor.normal.gene.value.mean.500 <- NULL;
for (i in 1:length(two.out.non.tumor.normal.gene.500)) {
    target.symbol <- two.out.non.tumor.normal.gene.500[i];

    # 1. outlier patient - tumor
    target.out.tumor.patient <- colnames(two.outlier.patient.status.merge.filter.500)[which(two.outlier.patient.status.merge.filter.500[target.symbol, ] == 1)];
    target.out.tumor.mean <- mean(na.omit(as.numeric(two.outlier.promoter.symbol.sample.match.merge.filter.500[target.symbol, target.out.tumor.patient])));

    # 2. outlier patient - normal
    #    - METABRIC
    target.out.normal.patient.meta <- sample.info.norm[sample.info.norm$matched_tumor %in% gsub('\\.', '_', target.out.tumor.patient), ]$samp;
    target.out.normal.patient.meta <- gsub('\\_', '.', target.out.normal.patient.meta);
    target.out.normal.meta <- na.omit(as.numeric(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500[target.symbol, target.out.normal.patient.meta]))
    #    - TCGA-BRCA
    target.out.normal.patient.brca <- colnames(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500)[substr(colnames(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500), 1, 12) %in% substr(target.out.tumor.patient, 1, 12)];
    target.out.normal.brca <- na.omit(as.numeric(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500[target.symbol, target.out.normal.patient.brca]))
    #   - combine
    target.out.normal.mean <- mean(na.omit(c(target.out.normal.meta, target.out.normal.brca)));

    # 3. non-outlier patient - tumor
    target.non.tumor.mean <- mean(na.omit(as.numeric(two.outlier.promoter.symbol.sample.match.merge.filter.500[target.symbol, -which(colnames(two.outlier.promoter.symbol.sample.match.merge.filter.500) %in% target.out.tumor.patient)])));

    # 4. non-outlier patient - normal
    target.non.normal.mean <- mean(na.omit(as.numeric(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500[target.symbol, -which(colnames(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500) %in% c(target.out.normal.patient.meta, target.out.normal.patient.brca))])));

    # Combine
    all.target.mean <- c(target.out.tumor.mean, target.out.normal.mean, target.non.tumor.mean, target.non.normal.mean);

    two.out.non.tumor.normal.gene.value.mean.500 <- rbind(two.out.non.tumor.normal.gene.value.mean.500, all.target.mean);
    }


rownames(two.out.non.tumor.normal.gene.value.mean.500) <- two.out.non.tumor.normal.gene.500;
two.out.non.tumor.normal.gene.value.mean.na.500 <- na.omit(two.out.non.tumor.normal.gene.value.mean.500);



# Use only differntially methylated genes
colnames(two.out.non.tumor.normal.gene.value.mean.na.500) <- c('outlier_tumor', 'outlier_normal', 'non_outlier_tumor', 'non_outlier_normal');
diff.outlier.normal.tumor <- two.out.non.tumor.normal.gene.value.mean.na.500[, 'outlier_normal'] - two.out.non.tumor.normal.gene.value.mean.na.500[, 'outlier_tumor'];
diff.non.outlier.tumor.outlier.tumor <- two.out.non.tumor.normal.gene.value.mean.na.500[, 'non_outlier_tumor'] - two.out.non.tumor.normal.gene.value.mean.na.500[, 'outlier_tumor'];


threshold <- 0.2;
two.out.non.tumor.normal.gene.value.mean.na.500.02 <- two.out.non.tumor.normal.gene.value.mean.na.500[
    diff.outlier.normal.tumor > threshold &
        diff.non.outlier.tumor.outlier.tumor > threshold,
    ];


# boxplot data
box.data <- data.frame(
    value = as.vector(two.out.non.tumor.normal.gene.value.mean.na.500.02),
    sample = rep(c('a', 'b', 'c', 'd'), each = nrow(two.out.non.tumor.normal.gene.value.mean.na.500.02))
    );

# Kruskal-Wallis test
kruskal.result.tumor.na <- kruskal.test(value ~ sample, data = box.data);
text.pvalue.tumor.na <- display.statistical.result(
    x = kruskal.result.tumor.na$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure2k')));

# Create boxplot
tumor.normal.box.plot <- BoutrosLab.plotting.general::create.boxplot(
    formula = value ~ sample,
    data = box.data,
    xlab.label = NULL,
    ylab.label = expression(paste('Mean of ', beta, ' value')),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    outliers = FALSE,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    ylimits = c(-0.05, 1.1),
    yat = seq(0, 1, 0.2),
    key = list(
        text = list(lab = text.pvalue.tumor.na, cex = 1),
        x = 0.55,
        y = 0.95
        ),
    add.stripplot = TRUE,
    points.pch = 1,
    points.cex = 0.8,
    points.col = 'grey50',
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 3.5),
    xright.rectangle = c(2.5, 5),
    ybottom.rectangle = -2,
    ytop.rectangle = 3,
    # set rectangle colour
    col.rectangle = 'grey',
    # set rectangle alpha (transparency)
    alpha.rectangle = 0.25,
    lwd = 1.2,
    col = c('red3', 'gold3', 'dodgerblue4', 'darkgreen'),
    alpha = 0.4
    );


save.outlier.figure(
    tumor.normal.box.plot,
    c('Figure2k', 'merge', 'tumour', 'normal', 'me', 'box'),
    width = 3.7,
    height = 5
    );
