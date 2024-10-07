### HISTORY ######################################################################
# This script processes NGF and LRP4 gene expression and DNA methylation data, and generates
# a scatter plot to visualize the relationship between RNA abundance and DNA
# methylation levels across different patients.
# The codes are connected to Figure 2h.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script analyzes the relationship between gene expression (RNA abundance) and
# DNA methylation for NGF and LRP4 genes.

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'))

# Load required data
load(file.path(get.outlier.data.dir(), '2024-10-03_Figure1_2_3_4_min_input.rda'))

### DESCRIPTION #################################################################
# The function `do.plot.2i` processes data for RNA abundance (FPKM) and DNA methylation
# for two example genes (NGF, LRP4) and creates scatter plots for visualization.

# Helper function to process FPKM data
process.fpkm.data <- function(fpkm.brca, fpkm.meta, gene) {
    fpkm.brca.gene <- fpkm.brca[fpkm.brca$Symbol %in% gene, -ncol(fpkm.brca), drop = FALSE]
    fpkm.meta.gene <- fpkm.meta[fpkm.meta$Symbol %in% gene, -ncol(fpkm.meta), drop = FALSE]

    fpkm.combined <- c(scale(as.numeric(fpkm.brca.gene)), scale(as.numeric(fpkm.meta.gene)))
    fpkm.df <- data.frame(t(fpkm.combined))
    rownames(fpkm.df) <- gene
    colnames(fpkm.df) <- c(colnames(fpkm.brca.gene), colnames(fpkm.meta.gene))

    return(fpkm.df)
    }

# Function to create the scatter plot for given gene
do.plot.2i <- function(gene) {
    # Process methylation and patient data
    gene.methyl <- two.outlier.promoter.symbol.sample.match.merge.filter.500[gene, , drop = FALSE]
    gene.patient <- two.outlier.patient.status.merge.filter.500[gene, ]

    # Process FPKM data for the gene
    fpkm.gene <- process.fpkm.data(fpkm.tumor.symbol.filter.brca, fpkm.tumor.symbol.filter.meta.symbol, gene)
    fpkm.gene.ordered <- fpkm.gene[, colnames(two.outlier.patient.status.merge.filter.500), drop = FALSE]
    fpkm.gene.ordered <- fpkm.gene.ordered[, order(as.numeric(fpkm.gene.ordered[1, ]), decreasing = TRUE), drop = FALSE]

    # Order methylation data
    gene.methyl.ordered <- gene.methyl[, colnames(fpkm.gene.ordered)]

    # Prepare scatter plot data
    scatter.data <- data.frame(fpkm = as.numeric(fpkm.gene.ordered), me = as.numeric(gene.methyl.ordered))
    rownames(scatter.data) <- colnames(gene.methyl.ordered)
    gene.patient.ordered <- gene.patient[colnames(gene.methyl.ordered)]

    # Define colors for plot points
    dot.colors <- ifelse(gene.patient.ordered == 1, 'red2', 'black')

    # Reverse order for scatter plot
    scatter.data.rev <- scatter.data[rev(seq(nrow(scatter.data))), ]
    dot.colors.rev <- rev(dot.colors)

    # Create scatter plot
    scatter.plot <- create.scatterplot(
        formula = fpkm ~ me,
        data = scatter.data.rev,
        col = dot.colors.rev,
        alpha = .6,
        xlimits = c(-0.06, 1.07),
        xaxis.fontface = 1,
        yaxis.fontface = 1,
        yaxis.tck = c(0.2, 0),
        xaxis.tck = c(0.2, 0),
        add.grid = TRUE,
        grid.colour = 'grey80',
        cex = 0.9,
        main.cex = 1.6,
        xaxis.cex = 1,
        yaxis.cex = 1,
        main = gene,
        xlab.cex = 1.3,
        ylab.cex = 1.3,
        ylab.label = expression(paste('RNA abundance (z-score)')),
        xlab.label = expression(paste('DNA methylation (', beta, ' value)')),
        type = c('p', 'r', 'g'),
        legend = list(
            inside = list(
                fun = draw.key,
                args = list(
                    key = get.corr.key(
                        x = scatter.data$fpkm,
                        y = scatter.data$me,
                        label.items = c('spearman'),
                        alpha.background = 0,
                        key.cex = 1.1
                        )
                    ),
                x = 0.75,
                y = 0.95,
                corner = c(0, 1)
                )
            )
        )

    # Save scatter plot
    save.outlier.figure(
        scatter.plot,
        c('Figure2i', gene, 'scatter'),
        width = 6,
        height = 6
        )
    }

# Generate scatter plots for NGF and LRP4 genes
do.plot.2i('NGF')
do.plot.2i('LRP4')

# Save session profile
save.session.profile(file.path('output', 'Figure2i.txt'))
