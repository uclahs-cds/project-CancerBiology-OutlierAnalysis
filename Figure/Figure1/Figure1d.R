### HISTORY #####################################################################
# This script calculates the sum of outlier genes per patient and per gene, 
# and then generates a histogram that displays the distribution of the number of outlier genes 
# per patient on a logarithmic scale.
# Date: 2024-08-12

### DESCRIPTION #################################################################
# The script processes a matrix of outlier genes across multiple patients. It sums 
# the outlier genes per patient and per gene, and then plots a histogram to visualize 
# the distribution of outlier genes per patient.

### PREAMBLE ####################################################################
# Calculate the sum of outlier genes per patient and per gene
outlier.patient.all.five.01.sum <- apply(
    outlier.patient.all.five.01, 
    1, 
    function(x) { sum(na.omit(x)) }
    );

outlier.patient.all.five.01.sum.gene <- apply(
    outlier.patient.all.five.01, 
    2, 
    function(x) { sum(na.omit(x)) }
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

# Save the histogram as a PDF
pdf(
    file = generate.filename(
        '5_patient_per_outlier_gene_number', 
        'histogram', 
        'pdf'
        ), 
    width = 5.5, 
    height = 5
    );
five.outlier.gene.sum;
dev.off();

# Save the histogram as a PNG
png(
    file = generate.filename(
        '5_patient_per_outlier_gene_number', 
        'histogram', 
        'png'
        ), 
    width = 5.5, 
    height = 5, 
    unit = 'in', 
    res = 1200
    );
five.outlier.gene.sum;
dev.off();
