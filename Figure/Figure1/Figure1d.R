### HISTORY #####################################################################
# This script calculates the sum of outlier genes per patient and per gene,
# and then generates a histogram that displays the distribution of the number of outlier genes
# per patient on a logarithmic scale.
# Date: 2024-08-12

### DESCRIPTION #################################################################
# The script processes a matrix of outlier genes across multiple patients. It sums
# the outlier genes per patient and per gene, and then plots a histogram to visualize
# the distribution of outlier genes per patient.

library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-09-10_Figure1.rda'));

### PREAMBLE ####################################################################

five.data.outlier.symbol.na <- na.omit(five.data.outlier.symbol);

### 1. MATADOR
# Find the position of "_" in Metador data row names
outlier.patient.tag.01.metador.pos <- which(
    strsplit(rownames(outlier.patient.tag.01.metador), '')[[1]] == '_'
    );
# Extract the substring after "_"
outlier.patient.tag.01.metador.symbol <- substring(
    rownames(outlier.patient.tag.01.metador),
    outlier.patient.tag.01.metador.pos + 1
    );
# Match and extract rows from Metador data based on 'five.data.outlier.symbol.na'
outlier.patient.tag.01.metador.match.five <- outlier.patient.tag.01.metador[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.metador.symbol),
    ];

### 2. TCGA-BRCA
outlier.patient.tag.01.brca.symbol <- fpkm.tumor.symbol.filter.brca[
    rownames(outlier.patient.tag.01.brca),
    ]$Symbol;
outlier.patient.tag.01.brca.match.five <- outlier.patient.tag.01.brca[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.brca.symbol),
    ];


### 3. METABRIC
outlier.patient.tag.01.meta.symbol <- fpkm.tumor.symbol.filter.meta.symbol[
    rownames(outlier.patient.tag.01.meta),
    ]$Symbol;
outlier.patient.tag.01.meta.match.five <- outlier.patient.tag.01.meta[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.meta.symbol),
    ];


### 4. ISPY-2
outlier.patient.tag.01.ispy.symbol <- rownames(outlier.patient.tag.01.ispy);
outlier.patient.tag.01.ispy.match.five <- outlier.patient.tag.01.ispy[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.ispy.symbol),
    ];


### 5. ICGC BRCA-EU
outlier.patient.tag.01.icgc.symbol <- fpkm.tumor.symbol.filter.symbol.icgc[
    rownames(outlier.patient.tag.01.icgc),
    ]$Symbol;
# Match and extract rows from ICGC data based on 'five.data.outlier.symbol.na'
outlier.patient.tag.01.icgc.match.five <- outlier.patient.tag.01.icgc[
    match(five.data.outlier.symbol.na, outlier.patient.tag.01.icgc.symbol),
    ];



# Combine matched data from all sources into a single data frame
outlier.patient.all.five.01 <- data.frame(
    cbind(
        outlier.patient.tag.01.brca.match.five,
        outlier.patient.tag.01.meta.match.five,
        outlier.patient.tag.01.ispy.match.five,
        outlier.patient.tag.01.metador.match.five,
        outlier.patient.tag.01.icgc.match.five
        )
    );

# Set row names of the combined data frame
rownames(outlier.patient.all.five.01) <- five.data.outlier.symbol.na;


# Calculate the sum of outlier genes per patient and per gene
outlier.patient.all.five.01.sum.gene <- apply(
    outlier.patient.all.five.01,
    2,
    function(x) {
        sum(na.omit(x))
        }
    );


# Create a histogram for the number of outlier genes per patient
five.outlier.gene.sum <- BoutrosLab.plotting.general::create.histogram(
    log10(outlier.patient.all.five.01.sum.gene + 1),
    breaks = 20,
    ylab.label = expression('Percent'),
    xlab.label = expression('Number of outlier genes per patient'),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.5,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    xat = c(0, log10(11), log10(101), log10(1001)),
    xaxis.lab = c(0, 10, 100, 1000),
    main = NULL
    );

# Display the histogram
five.outlier.gene.sum;

### OUTPUT ######################################################################
save.outlier.figure(
    five.outlier.gene.sum,
    c('Figure1d', '5_patient_per_outlier_gene_number', 'histogram'),
    width = 5.5,
    height = 5
    );

save.session.profile(file.path('output', 'Figure1d.txt'));
