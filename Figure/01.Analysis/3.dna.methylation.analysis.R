library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

# Load required variables from cache
load.multiple.computed.variables(c(
    'outlier.symbol'
    ));

# Get DNA methylation data
# 1. TCGA-BRCA
brca.me.outlier.match <- brca.me.data[rownames(brca.me.data) %in% outlier.symbol$brca, ]; # Filter BRCA methylation data for outliers
brca.outlier.non.promoter.symbol.sample.match.merge.500 <- brca.me.data[!(rownames(brca.me.data) %in% outlier.symbol$brca), ]; # Non-outlier BRCA data

# 2. METABRIC
meta.me.outlier.match <- meta.me.data[rownames(meta.me.data) %in% outlier.symbol$metabric, ] # Filter METABRIC methylation data for outliers
meta.outlier.non.promoter.symbol.sample.match.merge.500 <- meta.me.data[!(rownames(meta.me.data) %in% outlier.symbol$metabric), ] # Non-outlier METABRIC data

# Combine BRCA and METABRIC promoter data
me.out.symbol.two.500 <- unique(c(rownames(meta.me.outlier.match), rownames(brca.me.outlier.match))) # Unique symbols in both datasets

# Merge BRCA and METABRIC methylation data
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- lapply(me.out.symbol.two.500, function(symbol) {
    # Get BRCA data if available, otherwise fill with NAs
    target_gene_brca <- if (symbol %in% rownames(brca.me.outlier.match)) {
        as.numeric(brca.me.outlier.match[symbol, ])
        } else {
        rep(NA, ncol(brca.me.outlier.match))
        }

    # Get METABRIC data if available, otherwise fill with NAs
    target_gene_meta <- if (symbol %in% rownames(meta.me.outlier.match)) {
        as.numeric(meta.me.outlier.match[symbol, ])
        } else {
        rep(NA, ncol(meta.me.outlier.match))
        }

    # Combine BRCA and METABRIC data
    c(target_gene_brca, target_gene_meta)
    })

# Convert the list to a data frame and set row/column names
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- do.call(rbind, two.outlier.promoter.symbol.sample.match.merge.filter.500)
rownames(two.outlier.promoter.symbol.sample.match.merge.filter.500) <- me.out.symbol.two.500
colnames(two.outlier.promoter.symbol.sample.match.merge.filter.500) <- c(colnames(brca.me.outlier.match), colnames(meta.me.outlier.match))

# Convert to data frame if necessary
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- as.data.frame(two.outlier.promoter.symbol.sample.match.merge.filter.500)

# Merge outlier status data for BRCA and METABRIC
outlier.patient.tag.01.brca.me.match <- outlier.patient.tag.01.brca[, colnames(brca.me.outlier.match)]
outlier.patient.tag.01.brca.me.match <- outlier.patient.tag.01.brca.me.match[
    rownames(outlier.patient.tag.01.brca.me.match) %in% rownames(fpkm.tumor.symbol.filter.brca)[fpkm.tumor.symbol.filter.brca$Symbol %in% rownames(brca.me.outlier.match)],
    ]

outlier.patient.tag.01.meta.me.match <- outlier.patient.tag.01.meta[, colnames(meta.me.outlier.match)]
outlier.patient.tag.01.meta.me.match <- outlier.patient.tag.01.meta.me.match[
    rownames(outlier.patient.tag.01.meta.me.match) %in% rownames(fpkm.tumor.symbol.filter.meta.symbol)[fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% rownames(meta.me.outlier.match)],
    ]

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

# Convert the list into a data frame
two.outlier.patient.status.merge.filter.500 <- do.call(rbind, two.outlier.patient.status.merge.filter.list.500)
rownames(two.outlier.patient.status.merge.filter.500) <- me.out.symbol.two.500
colnames(two.outlier.patient.status.merge.filter.500) <- c(colnames(brca.me.data), colnames(meta.me.data))

# Function to split outliers and non-outliers
split_outliers <- function(gene_row, promoter_row) {
    list(
        outlier = promoter_row[gene_row == 1],
        non_outlier = promoter_row[gene_row == 0]
        )
    }

# Apply the split_outliers function to each row
outlier_non_outlier_list <- Map(split_outliers, as.data.frame(t(two.outlier.patient.status.merge.filter.500)), as.data.frame(t(two.outlier.promoter.symbol.sample.match.merge.filter.500)))

# Separate outliers and non-outliers into lists
outlier.sample.me.two.500 <- lapply(outlier_non_outlier_list, `[[`, 'outlier')
non.outlier.sample.me.two.500 <- lapply(outlier_non_outlier_list, `[[`, 'non_outlier')

# Convert lists to numeric vectors and calculate means
outlier.sample.me.two.unlist.500 <- as.numeric(unlist(outlier.sample.me.two.500))
non.outlier.sample.me.two.unlist.500 <- as.numeric(unlist(non.outlier.sample.me.two.500))

outlier.sample.me.two.unlist.mean.500 <- lapply(outlier.sample.me.two.500, function(x) mean(na.omit(as.numeric(x))))
non.outlier.sample.me.two.unlist.mean.500 <- lapply(non.outlier.sample.me.two.500, function(x) mean(na.omit(as.numeric(x))))

# Calculate mean beta values and their differences
mean.beta.merge.two.500 <- apply(two.outlier.promoter.symbol.sample.match.merge.filter.500, 1, function(x) mean(na.omit(as.numeric(x))))
minus.beta.merge.two.500 <- as.numeric(outlier.sample.me.two.unlist.mean.500) - as.numeric(non.outlier.sample.me.two.unlist.mean.500)

# Combine results into a data frame
mean.minus.ma.merge.two.500 <- data.frame(
    mean.beta = as.numeric(mean.beta.merge.two.500),
    minus.beta = as.numeric(minus.beta.merge.two.500),
    Symbol = rownames(two.outlier.promoter.symbol.sample.match.merge.filter.500)
    )

# Perform a Wilcoxon test between outliers and non-outliers
p.me <- wilcox.test(
    outlier.sample.me.two.unlist.500,
    non.outlier.sample.me.two.unlist.500,
    alternative = 'two.sided',
    conf.int = TRUE
    );

# Cache the computed variables for future use
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

# Save the session profile
save.session.profile(file.path('output', '3.dna.methylation.analysis.txt'));
