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
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

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


outlier.patient.tag.01.matador.gene.per.patient.sum <- apply(outlier.patient.tag.01.matador, 2, sum);
matador.outlier.patient.number.frame <- data.frame(
    sample = rep('matador', length(patient.part.matador)),
    value = as.numeric(outlier.patient.tag.01.matador.gene.per.patient.sum)
    );


outlier.patient.tag.01.icgc.gene.per.patient.sum <- apply(outlier.patient.tag.01.icgc, 2, sum);
icgc.outlier.patient.number.frame <- data.frame(
    sample = rep('ICGC BRCA-EU', length(patient.part.icgc)),
    value = as.numeric(outlier.patient.tag.01.icgc.gene.per.patient.sum)
    );


outlier.patient.tag.01.cheng.gene.per.patient.sum <- apply(outlier.patient.tag.01.cheng, 2, sum);
cheng.outlier.patient.number.frame <- data.frame(
    sample = rep('Cheng', length(patient.part.cheng)),
    value = as.numeric(outlier.patient.tag.01.cheng.gene.per.patient.sum)
    );


outlier.patient.tag.01.kao.gene.per.patient.sum <- apply(outlier.patient.tag.01.kao, 2, sum);
kao.outlier.patient.number.frame <- data.frame(
    sample = rep('kao', length(patient.part.kao)),
    value = as.numeric(outlier.patient.tag.01.kao.gene.per.patient.sum)
    );


outlier.patient.tag.01.hatzis.gene.per.patient.sum <- apply(outlier.patient.tag.01.hatzis, 2, sum);
hatzis.outlier.patient.number.frame <- data.frame(
    sample = rep('hatzis', length(patient.part.hatzis)),
    value = as.numeric(outlier.patient.tag.01.hatzis.gene.per.patient.sum)
    );


outlier.patient.tag.01.sjostrom.gene.per.patient.sum <- apply(outlier.patient.tag.01.sjostrom, 2, sum);
sjostrom.outlier.patient.number.frame <- data.frame(
    sample = rep('sjostrom', length(patient.part.sjostrom)),
    value = as.numeric(outlier.patient.tag.01.sjostrom.gene.per.patient.sum)
    );




# Combine all datasets into a single data frame
nine.outlier.patient.number.frame <- rbind(
    meta.outlier.patient.number.frame,
    brca.outlier.patient.number.frame,
    ispy.outlier.patient.number.frame,
    sjostrom.outlier.patient.number.frame,
    cheng.outlier.patient.number.frame,
    matador.outlier.patient.number.frame,
    icgc.outlier.patient.number.frame,
    kao.outlier.patient.number.frame,
    hatzis.outlier.patient.number.frame
    );

# Add an order column to distinguish between different datasets
nine.outlier.patient.number.frame.order <- cbind(
    nine.outlier.patient.number.frame,
    order = c(
        rep('a', nrow(meta.outlier.patient.number.frame)),
        rep('b', nrow(brca.outlier.patient.number.frame)),
        rep('c', nrow(ispy.outlier.patient.number.frame)),
        rep('d', nrow(sjostrom.outlier.patient.number.frame)),
        rep('e', nrow(cheng.outlier.patient.number.frame)),
        rep('f', nrow(matador.outlier.patient.number.frame)),
        rep('g', nrow(icgc.outlier.patient.number.frame)),
        rep('h', nrow(kao.outlier.patient.number.frame)),
        rep('i', nrow(hatzis.outlier.patient.number.frame))
        )
    );



### VIOLIN PLOT #################################################################

# Define colors for the plot
nine.col <- c(
    grDevices::adjustcolor('deepskyblue4', alpha.f = 0.7),
    grDevices::adjustcolor('firebrick3', alpha.f = 0.7),
    grDevices::adjustcolor('gold2', alpha.f = 0.7),
    grDevices::adjustcolor('darkgreen', alpha.f = 0.7),
    grDevices::adjustcolor('mediumpurple3', alpha.f = 0.7),
    grDevices::adjustcolor(c('darkorange2'), alpha.f = 0.7), 
    grDevices::adjustcolor(c('darkblue'), alpha.f = 0.7), 
    grDevices::adjustcolor(c('lightpink2'), alpha.f = 0.8),
    grDevices::adjustcolor(c('yellowgreen'), alpha.f = 0.8)
    );

# Create a violin plot to visualize the number of outlier genes per patient
nine.outlier.patient.violin <- BoutrosLab.plotting.general::create.violinplot(
    formula = log2(value + 1) ~ order,
    data = nine.outlier.patient.number.frame.order,
    main = expression('Number of outlier genes per patient'),
    main.cex = 1.4,
    xaxis.cex = 1.1,
    yaxis.cex = 1,
    xaxis.lab = c(
        'METABRIC\n n = 1991',
        'TCGA-BRCA\n n = 1085',
        'I-SPY2\n n = 988',
        'Sjostrom\n n = 765', 
        'Cheng\n n = 638', 
        'matador\n n = 528',
        'ICGC BRCA-EU\n n = 342',
        'Kao\n n = 327', 
        'Hatzis\n n = 310'
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
    ylimits = c(-1.3, 11.5),
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 3.5, 5.5, 7.5),
    xright.rectangle = c(2.5, 4.5, 6.5, 8.5),
    ybottom.rectangle = -10,
    ytop.rectangle = 15,
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
    col = nine.col
    );

# Display the violin plot
nine.outlier.patient.violin;

### OUTPUT ######################################################################
save.outlier.figure(
    nine.outlier.patient.violin,
    c('Figure1e', '9_patient_number', 'violin'),
    width = 5.5,
    height = 5.3
    );

save.session.profile(file.path('output', 'Figure1e.txt'));
