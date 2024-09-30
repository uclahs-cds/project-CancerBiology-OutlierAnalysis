### HISTORY #####################################################################
# Original script was adapted for better readability and maintainability according to the
# Boutros Lab R coding standards.
# Date: 2024-08-12

### DESCRIPTION #################################################################
# This script combines data from multiple genes into a single data frame and
# creates a strip plot showing the RNA abundance distribution across different genes.
# The plot distinguishes between outlier and non-outlier patients using different
# colors and point shapes.

### PREAMBLE ####################################################################
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-09-10_Figure1.rda'));

### DATA PREPARATION ############################################################


genes <- c('IGF2', 'TMEM30A', 'NRAS', 'IGF2R', 'GAPDH', 'B2M');

# fpkm.tumor.symbol.filter.XXX: RNA abundance matrix with rows of genes and
# columns of patients. Row names are Ensembl IDs.

# fpkm.tumor.symbol.filter.XXX.symbol is the same, except it has one
# additional `Symbol` column of gene aliases.

# All of the RNA abundance matricies have a Symbol column _except_ ispy, which
# has the symbols as its row names. Standardize that now.
fpkm.tumor.symbol.filter.ispy$Symbol <- rownames(fpkm.tumor.symbol.filter.ispy);

# outlier.patient.tag.01.XXX: Outlier status matrix of XXX dataset with
# rows of genes and columns of patients. Values are 1 for outlier events
# and 0 otherwise.

unique.datasets <- list(
    gene.dataset.name = c(
        'fpkm.tumor.symbol.filter.metador.symbol',
        'fpkm.tumor.symbol.filter.meta.symbol',
        'fpkm.tumor.symbol.filter.brca',
        'fpkm.tumor.symbol.filter.ispy',
        'fpkm.tumor.symbol.filter.symbol.icgc'
        ),
    outlier.dataset.name = c(
        'outlier.patient.tag.01.metador',
        'outlier.patient.tag.01.meta',
        'outlier.patient.tag.01.brca',
        'outlier.patient.tag.01.ispy',
        'outlier.patient.tag.01.icgc'
        )
    );

gene.data <- do.call(rbind, lapply(seq_along(genes), function(gene.index) {
    # Extract data for the specified gene across different datasets
    gene.alias <- genes[gene.index];

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
            outlier.status = outlier.status,
            order = letters[gene.index]
            );
        }

    # Build a matrix of the single-gene-single-dataset dataframes, one for each
    # dataset
    all.dataset.gene.results <- mapply(
        subset.genes.and.outliers,
        unique.datasets$gene.dataset.name,
        unique.datasets$outlier.dataset.name,
        SIMPLIFY = FALSE
        );

    # Return a single-gene-all-dataset dataframe
    do.call(rbind, all.dataset.gene.results);
    }));

# gene.data is now a single dataframe with one row for each patient-gene
# combination, across all datasets.
#
# # Define point colors, shapes, and sizes based on outlier status
gene.data$color <- ifelse(gene.data$outlier.status == 1, 'red2', 'black');
gene.data$shape <- ifelse(gene.data$outlier.status == 1, 23, 21);
gene.data$size <- ifelse(gene.data$outlier.status == 1, 0.9, 0.75);

# DATA ANALYSIS ###############################################################

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure1b')));

# Create the strip plot
stripplot.gene.z.scores <- BoutrosLab.plotting.general::create.stripplot(
    formula = as.numeric(z.score) ~ order,
    data = gene.data,
    xaxis.cex = 1.1,
    yaxis.cex = 1,
    xaxis.lab = genes,
    yat = seq(0, 200, 20),
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    cex = gene.data$size,
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 3.5, 5.5),
    xright.rectangle = c(2.5, 4.5, 6.5),
    ybottom.rectangle = -50,
    ytop.rectangle = 1000,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    ylab.label = expression('z-score'),
    xlab.label = NULL,
    col = gene.data$color,
    col.border = 'black',
    fill = 'transparent',
    xaxis.rot = 90,
    colour.alpha = 0.95,
    main.cex = 1.4,
    jitter.data = TRUE,
    jitter.factor = 1.1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    pch = gene.data$shape,
    main = expression('Distribution of RNA abundance'),
    bottom.padding = 6.5,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        fill = c('red2', 'grey5'),
                        pch = c(23, 21)
                        ),
                    text = list(
                        lab = c('Outlier Patient', 'Non-outlier Patient')
                        ),
                    padding.text = 3,
                    cex = 1.1
                    )
                ),
            x = 0.04,
            y = -0.34
            )
        )
    );


stripplot.gene.z.scores;

### OUTPUT ######################################################################

save.outlier.figure(
    stripplot.gene.z.scores,
    c('Figure1b', 'six_gene', 'scatter'),
    width = 4.7,
    height = 5.3
    );

save.session.profile(file.path('output', 'Figure1b.txt'));
