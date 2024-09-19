### HISTORY ######################################################################
# This script calculates the z-scores for the expression of the FGFR2 gene across different
# datasets. It identifies outlier and non-outlier samples and generates bar plots based on
# the calculated z-scores.
# Date: 2024-08-13

### PREAMBLE ####################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-08-23_Figure2a-d.rda'));

unique.datasets <- list(
    gene.dataset.name = c(
        'fpkm.tumor.symbol.filter.meta.symbol',
        'fpkm.tumor.symbol.filter.brca',
        'fpkm.tumor.symbol.filter.symbol.icgc'
        ),
    outlier.dataset.name = c(
        'outlier.patient.tag.01.meta',
        'outlier.patient.tag.01.brca',
        'outlier.patient.tag.01.icgc'
        )
    );

gene.alias <- "FGFR2"

subset.genes.and.outliers <- function(gene.dataset.name, outlier.dataset.name) {
    gene.dataset <- get(gene.dataset.name);
    outlier.dataset <- get(outlier.dataset.name);

    # Get a subset of the input gene data, eliminating all rows that don't
    # correspond to the gene of interest (or are NA) and the Symbol column.
    gene.subset <- gene.dataset[
        !is.na(gene.dataset$Symbol) & gene.alias == gene.dataset$Symbol,
        !('Symbol' == colnames(gene.dataset))
        ];

    # Get the outlier status of each patient for this gene, treating
    # missing data as not outlying.
    outlier.status <- as.numeric(outlier.dataset[rownames(gene.subset), ]);
    outlier.status[is.na(outlier.status)] <- 0;

    # Compute the mean and standard deviation for the non-outliers
    non.outliers.subset <- as.numeric(gene.subset[outlier.status == 0]);
    non.outliers.mean <- mean(non.outliers.subset);
    non.outliers.sd <- sd(non.outliers.subset);

    # Compute the z-score for all patients
    z.score <- as.numeric((gene.subset - non.outliers.mean) / non.outliers.sd);

    # Return a single-gene-single-dataset dataframe with one row per
    # patient
    data.frame(
        gene.alias = gene.alias,
        z.score = z.score,
        outlier.status = outlier.status
        # outlier.status = outlier.status,
        # order = letters[gene.index]
        );
    }

# Build a single-gene-all-dataset dataframe
gene.data <- do.call(rbind, mapply(
    subset.genes.and.outliers,
    unique.datasets$gene.dataset.name,
    unique.datasets$outlier.dataset.name,
    SIMPLIFY = FALSE
    ));

# Split data into outliers and non-outliers
outlier_z_scores <- gene.data[gene.data$outlier.status == 1, 'z.score'];
non_outlier_z_scores <- gene.data[gene.data$outlier.status == 0, 'z.score'];

# Quantile calculation and plotting data
# FIXME Wny are we not including the 100% quantile?
quantiles <- seq(0, 0.9, 0.1);

plot_data <- rbind(
        data.frame(
            quant = quantile(non_outlier_z_scores, p = quantiles),
            color = 'grey10'
            ),
        data.frame(
            quant = setNames(mean(outlier_z_scores), 'Outlier samples'),
            color = 'darkred'
        ));

# Add a factor column to the dataframe to control the order of the bars.
# Reverse the order of the levels so that the outlier samples are plotted first.
plot_data$label <- ordered(
    rownames(plot_data),
    levels = rev(rownames(plot_data))
    );

i.fpkm.merge.data.order.quan.unequal.three <- create.barplot(
    formula = quant ~ label,
    main = as.expression(substitute(paste(var), list(var = gene.alias))),
    ylab.label = expression('z-score'),
    xlab.label = NULL,
    xaxis.rot = 90,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.2,
    main.cex = 1.4,
    col = plot_data$color,
    border.col = 'black',
    border.lwd = 0.25,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    data = plot_data,
    ylimits = c(-12, 140)
    );

# Add points to plot
your.points <- xyplot(
    outlier_z_scores ~ c(0.9, 1, 1, 1.1, 1),
    pch = 23,
    col = 'black',
    fill = 'red2',
    cex = 1.3
    );
i.fpkm.merge.data.order.quan.unequal.three.dot <- i.fpkm.merge.data.order.quan.unequal.three + as.layer(your.points);

save.outlier.figure(
    i.fpkm.merge.data.order.quan.unequal.three.dot,
    c('Figure2b', gene.alias, 'barplot'),
    width = 4.5,
    height = 5
    );

save.session.profile(file.path('output', 'Figure2b.txt'));
