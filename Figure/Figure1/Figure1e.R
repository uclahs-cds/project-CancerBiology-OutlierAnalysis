### HISTORY #####################################################################
# This script generates a violin plot to visualize the number of outlier genes per patient
# across multiple datasets. The x-axis represents different datasets, and the y-axis
# shows the log2-transformed number of outlier genes.
# Date: 2024-08-12

### DESCRIPTION #################################################################
# The script processes data for the number of outlier genes per patient across
# multiple datasets. It then combines the data into a single data frame and creates
# a violin plot to compare the distributions of outlier genes per patient.

library(BoutrosLab.utilities);

# Source the helper library
args <- commandArgs();
source(file.path(
    dirname(dirname(normalizePath(sub('^--file=', '', args[grep('^--file=', args)])))),
    'common_functions.R'
    ));
# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-09-10_Figure1.rda'));

### PREAMBLE ####################################################################
# Create data frames for each dataset representing the number of outlier genes per patient
outlier.patient.tag.01.brca.gene.per.patient.sum <- apply(outlier.patient.tag.01.brca, 2, sum);
brca.outlier.patient.number.frame <- data.frame(
    sample = rep('TCGA-BRCA', length(patient.part.brca)),
    value = as.numeric(outlier.patient.tag.01.brca.gene.per.patient.sum)
    );


outlier.patient.tag.01.ispy.gene.per.patient.sum <- apply(outlier.patient.tag.01.ispy, 2, sum);
ispy.outlier.patient.number.frame <- data.frame(
    sample = rep('ISPY', length(patient.part.ispy)),
    value = as.numeric(outlier.patient.tag.01.ispy.gene.per.patient.sum)
    );


outlier.patient.tag.01.meta.gene.per.patient.sum <- apply(outlier.patient.tag.01.meta, 2, sum);
meta.outlier.patient.number.frame <- data.frame(
    sample = rep('METABRIC', length(patient.part.meta)),
    value = as.numeric(outlier.patient.tag.01.meta.gene.per.patient.sum)
    );


outlier.patient.tag.01.metador.gene.per.patient.sum <- apply(outlier.patient.tag.01.metador, 2, sum);
metador.outlier.patient.number.frame <- data.frame(
    sample = rep('METADOR', length(patient.part.metador)),
    value = as.numeric(outlier.patient.tag.01.metador.gene.per.patient.sum)
    );


outlier.patient.tag.01.icgc.gene.per.patient.sum <- apply(outlier.patient.tag.01.icgc, 2, sum);
icgc.outlier.patient.number.frame <- data.frame(
    sample = rep('ICGC BRCA-EU', length(patient.part.icgc)),
    value = as.numeric(outlier.patient.tag.01.icgc.gene.per.patient.sum)
    );

# Combine all datasets into a single data frame
five.outlier.patient.number.frame <- rbind(
    meta.outlier.patient.number.frame,
    brca.outlier.patient.number.frame,
    ispy.outlier.patient.number.frame,
    metador.outlier.patient.number.frame,
    icgc.outlier.patient.number.frame
    );

# Add an order column to distinguish between different datasets
five.outlier.patient.number.frame.order <- cbind(
    five.outlier.patient.number.frame,
    order = c(
        rep('a', nrow(meta.outlier.patient.number.frame)),
        rep('b', nrow(brca.outlier.patient.number.frame)),
        rep('c', nrow(ispy.outlier.patient.number.frame)),
        rep('d', nrow(metador.outlier.patient.number.frame)),
        rep('e', nrow(icgc.outlier.patient.number.frame))
        )
    );

### VIOLIN PLOT #################################################################

# Define colors for the plot
five.col <- c(
    grDevices::adjustcolor('deepskyblue4', alpha.f = 0.7),
    grDevices::adjustcolor('firebrick3', alpha.f = 0.7),
    grDevices::adjustcolor('gold2', alpha.f = 0.7),
    grDevices::adjustcolor('darkgreen', alpha.f = 0.7),
    grDevices::adjustcolor('mediumpurple', alpha.f = 0.7)
    );

# Create a violin plot to visualize the number of outlier genes per patient
five.outlier.patient.violin <- BoutrosLab.plotting.general::create.violinplot(
    formula = log2(value + 1) ~ order,
    data = five.outlier.patient.number.frame.order,
    main = expression('Number of outlier genes per patient'),
    main.cex = 1.4,
    xaxis.cex = 1.1,
    yaxis.cex = 1,
    xaxis.lab = c(
        'METABRIC\n n = 1991',
        'TCGA-BRCA\n n = 1085',
        'I-SPY2\n n = 988',
        'METADOR\n n = 528',
        'ICGC BRCA-EU\n n = 342'
        ),
    yaxis.lab = c(
        expression('2'^'0'),
        expression('2'^'2'),
        expression('2'^'4'),
        expression('2'^'6'),
        expression('2'^'8'),
        expression('2'^'10')
        ),
    yat = c(0, 2, 4, 6, 8, 10),
    ylimits = c(-1.2, 11),
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 3.5),
    xright.rectangle = c(2.5, 4.5),
    ybottom.rectangle = -50,
    ytop.rectangle = 1000,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    ylab.label = expression('Number of outlier genes'),
    xlab.label = NULL,
    xaxis.rot = 90,
    lwd = 1.2,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    col = five.col
    );

# Display the violin plot
five.outlier.patient.violin;

### OUTPUT ######################################################################
save.outlier.figure(
    five.outlier.patient.violin,
    c('Figure1e', '5_patient_number', 'violin'),
    width = 4.3,
    height = 5.3
    );

save.session.profile(file.path('output', 'Figure1e.txt'));
