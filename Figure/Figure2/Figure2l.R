### HISTORY ######################################################################
# This script analyzes the DNA methylation on the promoter region of the PXDNL
# gene in outlier and non-outlier patients in TCGA-BRCA data.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script focuses on analyzing DNA methylation patterns in the promoter region
# of the PXDNL gene. It compares methylation levels between outlier and non-outlier
# patients using TCGA-BRCA data. The script processes methylation data for tumor
# and normal samples, orders the data based on gene expression levels, and creates
# heatmaps to visualize methylation patterns across different patient groups and
# sample types.

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
library(metafor);

# Source helper library
source(here::here('common_functions.R'))

# Load the data file
load(file.path(get.outlier.data.dir(), '2024-10-03_Figure1_2_3_4_min_input.rda'))


load.multiple.computed.variables(c(
    'outlier.patient.tag.01.brca.me.match'
    ));

### DESCRIPTION #################################################################
# Function to process DNA methylation and FPKM data for the PDXNL gene and generate heatmaps.

# Helper function to process and order methylation data
process.methylation.data <- function(gene, outlier.data, normal.data, fpkm.data, tag.data) {
    # Methylation data for outliers
    gene.me <- outlier.data[outlier.data$Symbol == gene, !('Symbol' == colnames(outlier.data))];
    gene.me.normal <- normal.data[rownames(gene.me), !('Symbol' == colnames(normal.data))];

    # Patient grouping
    gene.patient <- tag.data[rownames(fpkm.data)[fpkm.data$Symbol == gene], ]
    gene.me.patient <- gene.me[, gene.patient == 1, drop = FALSE]
    gene.me.patient.non <- gene.me[, gene.patient == 0, drop = FALSE]

    overlap_patient <- colnames(normal.data)[substr(colnames(normal.data), 1, 12) %in% substr(names(gene.patient)[gene.patient == 1], 1, 12)]

    gene.me.patient_normal <- gene.me.normal[, overlap_patient, drop = FALSE]
    gene.me.patient.non_normal <- gene.me.normal[, !colnames(gene.me.normal) %in% overlap_patient, drop = FALSE]

    # Order the data by promoter position
    gene.promoters <- promoters.info[[gene]];
    gene.promoters.order <- rownames(gene.promoters[order(gene.promoters$pos), ]);

    list(
        gene.me.patient = gene.me.patient[gene.promoters.order, ],
        gene.me.patient.non = gene.me.patient.non[gene.promoters.order, ],
        gene.me.patient_normal = gene.me.patient_normal[gene.promoters.order, ],
        gene.me.patient.non_normal = gene.me.patient.non_normal[gene.promoters.order, ]
        )
    }

# Function to create a heatmap for a given methylation data set
create_gene_heatmap <- function(methylation_data, clustering = 'none', cluster.dimension = NULL, show.color.key = FALSE) {
    BoutrosLab.plotting.general:::create.heatmap(
        x = methylation_data,
        clustering.method = clustering,
        cluster.dimensions = cluster.dimension,
        plot.dendrograms = FALSE,
        colour.scheme = c('#b2182b', 'white', '#2166ac'),
        grid.row = FALSE,
        grid.col = FALSE,
        yaxis.tck = 0,
        xaxis.tck = 0,
        yaxis.cex = 0,
        yaxis.rot = 0,
        ylab.cex = 0,
        at = seq(0, 1, 0.001),
        colourkey.cex = 1.3,
        print.colour.key = show.color.key
        )
    }

# Process PXDNL data
gene <- 'PXDNL'
methylation_data <- process.methylation.data(
    gene = gene,
    outlier.data = brca.outlier.promoter.symbol.sample.match.brca,
    normal.data = brca.outlier.promoter.symbol.normal.match.filter.brca,
    fpkm.data = fpkm.tumor.symbol.filter.brca,
    tag.data = outlier.patient.tag.01.brca.me.match
    )

# Create heatmaps for different patient groups
heatmap.outlier <- create_gene_heatmap(t(methylation_data$gene.me.patient));
heatmap.outlier.normal <- create_gene_heatmap(t(methylation_data$gene.me.patient_normal));

heatmap.non.outlier <- create_gene_heatmap(
    methylation_data$gene.me.patient.non,
    clustering = 'ward.D2',
    cluster.dimension = 'row',
    show.color.key = TRUE
    );
heatmap.non.outlier.normal <- create_gene_heatmap(
    methylation_data$gene.me.patient.non_normal,
    clustering = 'ward.D2',
    cluster.dimension = 'row',
    show.color.key = TRUE
    );

# Combine heatmaps into a multiplot
combined_heatmap <- BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(heatmap.non.outlier.normal, heatmap.non.outlier, heatmap.outlier.normal, heatmap.outlier),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = gene,
    xlab.label = expression('Beta value'),
    ylab.label = c(expression('Outlier patient'), '', '', '', expression('Non-outlier patients'), '', '', '', ''),
    yaxis.fontface = 1,
    plot.layout = c(1, 4),
    main.key.padding = 3,
    panel.heights = c(0.08, 0.08, 1, 0.2),
    ylab.padding = 1,
    y.spacing = -0.7,
    main.cex = 1.6,
    xaxis.cex = 0,
    xlab.padding = -10,
    xlab.to.xaxis.padding = -1,
    bottom.padding = 3,
    right.padding = 3,
    yaxis.cex = 1.2,
    yaxis.tck = 0,
    ylab.cex = 1.15,
    xlab.cex = 1.3,
    xaxis.rot = 90,
    xaxis.tck = 0,
    xlab.key.padding = 4,
    resolution = 500
    )

# Save the heatmap plot
save.outlier.figure(
    combined_heatmap,
    c('Figure2l', gene, 'heatmap', 'me'),
    width = 7.5,
    height = 8.5
    )

# Save session profile
save.session.profile(file.path('output', 'Figure2l.txt'))
