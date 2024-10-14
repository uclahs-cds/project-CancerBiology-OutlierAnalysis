library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'outlier.symbol'
    ));

# Get DNA methylation data
# 1. TCGA-BRCA
brca.me.outlier.match <- brca.me.data[rownames(brca.me.data) %in% outlier.symbol$brca, ];
brca.outlier.non.promoter.symbol.sample.match.merge.500 <- brca.me.data[!(rownames(brca.me.data) %in% outlier.symbol$brca), ];
# 2. METABRIC
meta.me.outlier.match <- meta.me.data[rownames(meta.me.data) %in% outlier.symbol$metabric, ];
meta.outlier.non.promoter.symbol.sample.match.merge.500 <- meta.me.data[!(rownames(meta.me.data) %in% outlier.symbol$metabric), ];


# Use promoter region TSS ~ +500bp
me.out.symbol.two.500 <- unique(c(
    rownames(meta.me.outlier.match),
    rownames(brca.me.outlier.match)
    ));

# merge two dataset
#   - 1. merge methylation data
# Initialize the list to store merged data
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- lapply(me.out.symbol.two.500, function(symbol) {

    # Get BRCA target gene data or assign NA if not found
    target_gene_brca <- if (symbol %in% rownames(brca.me.outlier.match)) {
        as.numeric(brca.me.outlier.match[symbol, ])
    } else {
        rep(NA, ncol(brca.me.outlier.match))
    }

    # Get METABRIC target gene data or assign NA if not found
    target_gene_meta <- if (symbol %in% rownames(meta.me.outlier.match)) {
        as.numeric(meta.me.outlier.match[symbol, ])
    } else {
        rep(NA, ncol(meta.me.outlier.match))
    }

    # Combine BRCA and METABRIC data
    c(target_gene_brca, target_gene_meta)
})

# Combine the list into a data frame and set row and column names
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- do.call(rbind, two.outlier.promoter.symbol.sample.match.merge.filter.500)
rownames(two.outlier.promoter.symbol.sample.match.merge.filter.500) <- me.out.symbol.two.500
colnames(two.outlier.promoter.symbol.sample.match.merge.filter.500) <- c(colnames(brca.me.outlier.match), colnames(meta.me.outlier.match))

# Convert to data frame if necessary
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- as.data.frame(two.outlier.promoter.symbol.sample.match.merge.filter.500)

# 2. Merge outlier status data
outlier.patient.tag.01.brca.me.match <- outlier.patient.tag.01.brca[, colnames(brca.me.outlier.match)];
outlier.patient.tag.01.brca.me.match <- outlier.patient.tag.01.brca.me.match[
    rownames(outlier.patient.tag.01.brca.me.match) %in% rownames(fpkm.tumor.symbol.filter.brca)[
        fpkm.tumor.symbol.filter.brca$Symbol %in%
            rownames(brca.me.outlier.match)
        ],
    ];
outlier.patient.tag.01.meta.me.match <- outlier.patient.tag.01.meta[, colnames(meta.me.outlier.match)];
outlier.patient.tag.01.meta.me.match <- outlier.patient.tag.01.meta.me.match[
    rownames(outlier.patient.tag.01.meta.me.match) %in% rownames(fpkm.tumor.symbol.filter.meta.symbol)[
        fpkm.tumor.symbol.filter.meta.symbol$Symbol %in%
            rownames(meta.me.outlier.match)
        ],
    ];

two.outlier.patient.status.merge.filter.list.500 <- list();

for (i in 1:length(me.out.symbol.two.500)) {
    if (me.out.symbol.two.500[i] %in% rownames(brca.me.outlier.match)) {
        row.brca <- rownames(fpkm.tumor.symbol.filter.brca)[
            fpkm.tumor.symbol.filter.brca$Symbol == me.out.symbol.two.500[i]
            ];
        row.brca <- row.brca[row.brca %in% rownames(outlier.patient.tag.01.brca.me.match)];
        target.gene.brca <- as.numeric(
            outlier.patient.tag.01.brca.me.match[row.brca, ]
            );
        } else {
        target.gene.brca <- rep('NA', ncol(brca.me.outlier.match));
        }

    if (me.out.symbol.two.500[i] %in% rownames(meta.me.outlier.match)) {
        target.gene.meta <- as.numeric(
            outlier.patient.tag.01.meta.me.match[
                rownames(fpkm.tumor.symbol.filter.meta.symbol)[
                    fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% me.out.symbol.two.500[i]
                    ],
                ]
            );
        } else {
        target.gene.meta <- rep('NA', ncol(outlier.patient.tag.01.meta.me.match));
        }

    both.target.gene <- c(target.gene.brca, target.gene.meta);

    two.outlier.patient.status.merge.filter.list.500[[i]] <- both.target.gene;
    }

two.outlier.patient.status.merge.filter.500 <- do.call(rbind, two.outlier.patient.status.merge.filter.list.500);
two.outlier.patient.status.merge.filter.500 <- data.frame(two.outlier.patient.status.merge.filter.500);
rownames(two.outlier.patient.status.merge.filter.500) <- me.out.symbol.two.500;
colnames(two.outlier.patient.status.merge.filter.500) <- c(
    colnames(brca.me.data),
    colnames(meta.me.data)
    );

# Use lapply to split outliers and non-outliers efficiently
split_outliers <- function(gene_row, promoter_row) {
    list(
        outlier = promoter_row[gene_row == 1],
        non_outlier = promoter_row[gene_row == 0]
        )
    }

# Apply split_outliers function to each row
outlier_non_outlier_list <- Map(
    split_outliers,
    as.data.frame(t(two.outlier.patient.status.merge.filter.500)),
    as.data.frame(t(two.outlier.promoter.symbol.sample.match.merge.filter.500))
    )

# Split results into outliers and non-outliers
outlier.sample.me.two.500 <- lapply(outlier_non_outlier_list, `[[`, 'outlier')
non.outlier.sample.me.two.500 <- lapply(outlier_non_outlier_list, `[[`, 'non_outlier')

# Flatten lists into numeric vectors
outlier.sample.me.two.unlist.500 <- as.numeric(unlist(outlier.sample.me.two.500))
non.outlier.sample.me.two.unlist.500 <- as.numeric(unlist(non.outlier.sample.me.two.500))

outlier.sample.me.two.unlist.mean.500 <- lapply(outlier.sample.me.two.500, function(x) {
    mean(na.omit(as.numeric(x)))
    });
non.outlier.sample.me.two.unlist.mean.500 <- lapply(non.outlier.sample.me.two.500, function(x) {
    mean(na.omit(as.numeric(x)))
    });

mean.beta.merge.two.500 <- apply(
    two.outlier.promoter.symbol.sample.match.merge.filter.500,
    1,
    function(x) {
        mean(na.omit(as.numeric(x)))
        }
    );
minus.beta.merge.two.500 <- as.numeric(outlier.sample.me.two.unlist.mean.500) - as.numeric(non.outlier.sample.me.two.unlist.mean.500);
mean.minus.ma.merge.two.500 <- data.frame(cbind(
    as.numeric(mean.beta.merge.two.500),
    as.numeric(minus.beta.merge.two.500),
    rownames(two.outlier.promoter.symbol.sample.match.merge.filter.500)
    ));

mean.minus.ma.merge.two.500[, 1] <- as.numeric(mean.minus.ma.merge.two.500[, 1]);
mean.minus.ma.merge.two.500[, 2] <- as.numeric(mean.minus.ma.merge.two.500[, 2]);
colnames(mean.minus.ma.merge.two.500) <- c('mean.beta', 'minus.beta', 'Symbol');

p.me <- wilcox.test(
    outlier.sample.me.two.unlist.500,
    non.outlier.sample.me.two.unlist.500,
    alternative = 'two.sided',
    conf.int = TRUE
    );

cache.multiple.computed.variables(c(
    'p.me',
    'mean.minus.ma.merge.two.500',
    'brca.outlier.non.promoter.symbol.sample.match.merge.500',
    'me.out.symbol.two.500',
    'meta.me.outlier.match',
    'meta.outlier.non.promoter.symbol.sample.match.merge.500',
    'non.outlier.sample.me.two.500',
    'outlier.patient.tag.01.brca.me.match',
    'outlier.patient.tag.01.meta.me.match',
    'outlier.sample.me.two.500',
    'two.outlier.patient.status.merge.filter.500',
    'two.outlier.promoter.symbol.sample.match.merge.filter.500'
    ));


save.session.profile(file.path('output', '3.dna.methylation.analysis.txt'));
