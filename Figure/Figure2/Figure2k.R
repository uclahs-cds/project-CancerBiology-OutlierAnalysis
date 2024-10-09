### HISTORY ######################################################################
# This script processes gene expression data to compare outlier and non-outlier
# patients across tumor and normal samples. It generates a boxplot to visualize
# the mean beta values and performs a Kruskal-Wallis test for statistical analysis.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script analyzes gene expression data, comparing outlier and non-outlier
# patients in both tumor and normal samples. It calculates mean beta values,
# identifies differentially methylated genes, and visualizes the results using
# a boxplot. The script also performs a Kruskal-Wallis test for statistical
# significance and includes data from multiple sources (METABRIC, TCGA-BRCA).

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-10-08_Figure1_2_3_4_min_input.rda'));

load.multiple.computed.variables(c(
    'meta.me.outlier.match',
    'outlier.patient.tag.01.brca.me.match',
    'outlier.patient.tag.01.meta.me.match',
    'two.outlier.patient.status.merge.filter.500',
    'two.outlier.promoter.symbol.sample.match.merge.filter.500'
    ));


### Analyze the normal data
# 1. TCGA-BRCA
process.outliers <- function(df, outlier.patients) {
    lapply(1:nrow(df), function(i) {
        symbol.name <- rownames(df)[i];
        outlier.patient.gene <- outlier.patients[
            rownames(fpkm.tumor.symbol.filter.brca)[
                fpkm.tumor.symbol.filter.brca$Symbol == symbol.name
                ],
            ];
        overlap.patient <- colnames(df)[
            substr(colnames(df), 1, 12) %in%
                substr(names(outlier.patient.gene)[outlier.patient.gene == 1], 1, 12)
            ];
        df[i, overlap.patient];
        });
    }


process.sample.methylation <- function(df, outlier.patients) {
    lapply(1:nrow(df), function(i) {
        symbol.name <- rownames(df)[i];
        outlier.patient.gene <- outlier.patients[
            rownames(fpkm.tumor.symbol.filter.brca)[
                fpkm.tumor.symbol.filter.brca$Symbol == symbol.name
                ],
            ];
        outlier.patient.gene.na <- na.omit(outlier.patient.gene);
        list(
            outlier = df[i, ][outlier.patient.gene.na == 1],
            non.outlier = df[i, ][outlier.patient.gene.na == 0]
            );
        });
    }


# Process the outliers for the normal sample methylation data
outlier.normal.me.merge.500 <- process.outliers(
    brca.outlier.promoter.symbol.normal.match.merge.500,
    outlier.patient.tag.01.brca.me.match
    );

# Process sample methylation data to separate outliers from non-outliers
sample.methylation <- process.sample.methylation(
    brca.me.outlier.match,
    outlier.patient.tag.01.brca.me.match
    );

# Extract outliers and non-outliers
outlier.sample.me.merge.beta.500 <- lapply(
    sample.methylation,
    `[[`, 'outlier'
    );
non.outlier.sample.me.merge.beta.500 <- lapply(
    sample.methylation,
    `[[`, 'non.outlier'
    );

# Calculate the means for outlier and non-outlier samples
outlier.sample.me.merge.unlist.beta.500 <- as.numeric(unlist(outlier.sample.me.merge.beta.500));
non.outlier.sample.me.merge.unlist.beta.500 <- as.numeric(unlist(non.outlier.sample.me.merge.beta.500));
outlier.sample.me.merge.unlist.beta.mean.500 <- sapply(
    outlier.sample.me.merge.beta.500,
    function(x) mean(as.numeric(x))
    );
non.outlier.sample.me.merge.unlist.beta.mean.500 <- sapply(
    non.outlier.sample.me.merge.beta.500,
    function(x) mean(na.omit(as.numeric(x)))
    );

# Process normal samples and filter for relevant genes
outlier.sample.normal.beta.merge.500 <- sapply(
    outlier.normal.me.merge.500,
    mean
    );

names(outlier.sample.normal.beta.merge.500) <- rownames(brca.me.outlier.match);
outlier.sample.normal.beta.merge.unlist.filter.500 <- outlier.sample.normal.beta.merge.500[
    rownames(brca.me.outlier.match)
    ];

# Prepare comparison data between normal and tumor methylation samples
normal.tumor.beta.comparison.merge.500 <- data.frame(
    normal = as.numeric(outlier.sample.normal.beta.merge.unlist.filter.500),
    tumor = as.numeric(outlier.sample.me.merge.unlist.beta.mean.500)
    );
rownames(normal.tumor.beta.comparison.merge.500) <- rownames(brca.me.outlier.match);
normal.tumor.beta.comparison.merge.500 <- na.omit(normal.tumor.beta.comparison.merge.500);

# Calculate the difference between normal and tumor samples
matched.me.outlier.normal.unlist.minus.name.na.brca.500 <- cbind(
    normal.tumor.beta.comparison.merge.500,
    minus = normal.tumor.beta.comparison.merge.500[, 2] - normal.tumor.beta.comparison.merge.500[, 1]
    );

# Order the data by the 'minus' column
matched.me.outlier.normal.unlist.minus.name.na.order.brca.500 <- matched.me.outlier.normal.unlist.minus.name.na.brca.500[
    order(matched.me.outlier.normal.unlist.minus.name.na.brca.500$minus),
    ];



# 2. METABRIC
# Filter the outlier genes having matched normal
meta.com.outlier.promoter.symbol.normal.match.name <- sub('^chr\\d+[.+-]', '', rownames(meta.com.outlier.promoter.symbol.normal.match));
meta.com.outlier.promoter.symbol.normal.match.symbol <- sub('^.*[+-]', '', meta.com.outlier.promoter.symbol.normal.match.name);
meta.com.outlier.promoter.symbol.normal.match.symbol <- sub('\\.\\d+$', '', meta.com.outlier.promoter.symbol.normal.match.symbol);
meta.com.outlier.promoter.symbol.normal.match.symbol <- sub('^\\.', '', meta.com.outlier.promoter.symbol.normal.match.symbol);


meta.com.outlier.promoter.symbol.normal.match.meta.name <- meta.com.outlier.promoter.symbol.normal.match;
meta.com.outlier.promoter.symbol.normal.match.meta.name$name <- meta.com.outlier.promoter.symbol.normal.match.symbol;


meta.com.outlier.promoter.symbol.normal.match.filter.meta.name.unique <- na.omit(unique(meta.com.outlier.promoter.symbol.normal.match.meta.name$name));
meta.com.outlier.promoter.symbol.normal.match.meta.unique <- NULL;
for (i in 1:length(meta.com.outlier.promoter.symbol.normal.match.filter.meta.name.unique)) {
    target.gene <- meta.com.outlier.promoter.symbol.normal.match.filter.meta.name.unique[i];
    target.beta <- meta.com.outlier.promoter.symbol.normal.match[
        meta.com.outlier.promoter.symbol.normal.match.meta.name$name %in% target.gene,
        ];
    if (1 == nrow(target.beta)) {
        target.beta.mean <- target.beta;
        } else {
        target.beta.mean <- apply(target.beta, 2, function(x) {
            mean(na.omit(x))
            });
        }
    meta.com.outlier.promoter.symbol.normal.match.meta.unique <- rbind(
        meta.com.outlier.promoter.symbol.normal.match.meta.unique,
        target.beta.mean
        );
    }
rownames(meta.com.outlier.promoter.symbol.normal.match.meta.unique) <- meta.com.outlier.promoter.symbol.normal.match.filter.meta.name.unique;


meta.com.matched.normal <- list();
meta.com.matched.outlier <- list();
meta.com.matched.symbol <- NULL;

for (i in 1:nrow(outlier.patient.tag.01.meta.me.match)) {
    target.outlier.patient <- colnames(outlier.patient.tag.01.meta.me.match)[
        outlier.patient.tag.01.meta.me.match[i, ] == 1
        ];

    target.symbol <- fpkm.tumor.symbol.filter.meta.symbol[
        rownames(outlier.patient.tag.01.meta.me.match[i, ]), 'Symbol'
        ];

    target.level <- meta.me.outlier.match[
        rownames(meta.me.outlier.match) %in% target.symbol,
        ];

    target.level.normal <- meta.com.outlier.promoter.symbol.normal.match[
        meta.com.outlier.promoter.symbol.normal.match.symbol %in% target.symbol,
        ];

    target.normal.patient <- sample.info.norm[
        sample.info.norm$matched_tumor %in% gsub('\\.', '_', target.outlier.patient),
        ]$samp;

    target.outlier.patient.match <- gsub(
        '_', '.', sample.info.norm[sample.info.norm$samp %in% target.normal.patient, 'matched_tumor']
        );

    meta.com.matched.normal[[i]] <- target.level.normal[
        , gsub('_', '.', target.normal.patient),
        drop = FALSE
        ];

    meta.com.matched.outlier[[i]] <- target.level[
        , target.outlier.patient.match,
        drop = FALSE
        ];

    meta.com.matched.symbol <- c(meta.com.matched.symbol, target.symbol);
    }

meta.com.matched.normal.unlist <- unlist(lapply(meta.com.matched.normal, function(x) {
    mean(na.omit(as.numeric(unlist(x))))
    }));
meta.com.matched.outlier.unlist <- unlist(lapply(meta.com.matched.outlier, function(x) {
    mean(na.omit(as.numeric(unlist(x))))
    }));

normal.tumor.com.beta.comparison <- data.frame(cbind(
    meta.com.matched.normal.unlist,
    meta.com.matched.outlier.unlist,
    meta.com.matched.symbol
    ));
colnames(normal.tumor.com.beta.comparison) <- c('normal', 'tumor', 'Symbol');
normal.tumor.com.beta.comparison[normal.tumor.com.beta.comparison == 'NaN'] <- NA;
normal.tumor.com.beta.comparison <- na.omit(normal.tumor.com.beta.comparison);




# divide into outlier and non-outlier
me.out.normal.symbol.two.500 <- unique(c(
    rownames(matched.me.outlier.normal.unlist.minus.name.na.order.brca.500),
    normal.tumor.com.beta.comparison$Symbol
    ));

two.outlier.promoter.symbol.sample.normal.match.merge.filter.500 <- list();
for (i in 1:length(me.out.normal.symbol.two.500)) {
    if (me.out.normal.symbol.two.500[i] %in% rownames(brca.outlier.promoter.symbol.normal.match.merge.500)) {
        target.gene.brca <- as.numeric(brca.outlier.promoter.symbol.normal.match.merge.500[me.out.normal.symbol.two.500[i], ]);
        } else {
        target.gene.brca <- rep('NA', ncol(brca.outlier.promoter.symbol.normal.match.merge.500))
        }

    if (me.out.normal.symbol.two.500[i] %in% rownames(meta.com.outlier.promoter.symbol.normal.match.meta.unique)) {
        target.gene.meta <- as.numeric(meta.com.outlier.promoter.symbol.normal.match.meta.unique[me.out.normal.symbol.two.500[i], ]);
        } else {
        target.gene.meta <- rep('NA', ncol(meta.com.outlier.promoter.symbol.normal.match.meta.unique))
        }

    both.target.gene <- c(target.gene.brca, target.gene.meta)

    two.outlier.promoter.symbol.sample.normal.match.merge.filter.500[[i]] <- both.target.gene;
    }


two.outlier.promoter.symbol.sample.normal.match.merge.filter.500 <- do.call(rbind, two.outlier.promoter.symbol.sample.normal.match.merge.filter.500);
rownames(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500) <- me.out.normal.symbol.two.500;
colnames(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500) <- c(colnames(brca.outlier.promoter.symbol.normal.match.merge.500), colnames(meta.com.outlier.promoter.symbol.normal.match.meta.unique));
two.outlier.promoter.symbol.sample.normal.match.merge.filter.500 <- data.frame(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500);


two.outlier.promoter.symbol.sample.normal.match.merge.filter.mean.500 <- apply(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500, 1, function(x) {
    mean(na.omit(as.numeric(x)))
    });
two.outlier.promoter.symbol.sample.normal.match.merge.filter.mean.na.500 <- na.omit(two.outlier.promoter.symbol.sample.normal.match.merge.filter.mean.500);



two.out.non.tumor.normal.gene.value.mean.500 <- NULL;
for (i in 1:length(me.out.normal.symbol.two.500)) {
    target.symbol <- me.out.normal.symbol.two.500[i];

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


rownames(two.out.non.tumor.normal.gene.value.mean.500) <- me.out.normal.symbol.two.500;
two.out.non.tumor.normal.gene.value.mean.na.500 <- na.omit(two.out.non.tumor.normal.gene.value.mean.500);



# Use only differntially methylated genes
colnames(two.out.non.tumor.normal.gene.value.mean.na.500) <- c('outlier_tumor', 'outlier_normal', 'non_outlier_tumor', 'non_outlier_normal');
diff.outlier.normal.tumor <- two.out.non.tumor.normal.gene.value.mean.na.500[, 'outlier_normal'] - two.out.non.tumor.normal.gene.value.mean.na.500[, 'outlier_tumor'];
diff.non.outlier.tumor.outlier.tumor <- two.out.non.tumor.normal.gene.value.mean.na.500[, 'non_outlier_tumor'] - two.out.non.tumor.normal.gene.value.mean.na.500[, 'outlier_tumor'];

threshold <- 0.2;
two.out.non.tumor.normal.gene.value.mean.na.500.02 <- two.out.non.tumor.normal.gene.value.mean.na.500[
    # diff.outlier.normal.tumor > threshold &
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

save.session.profile(file.path('output', 'Figure2k.txt'));
