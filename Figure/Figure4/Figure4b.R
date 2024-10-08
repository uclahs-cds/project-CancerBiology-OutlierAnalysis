### HISTORY #####################################################################
# This script analyzes the number of outlier genes per patient from CCLE and tissue
# datasets.
# Date: 2024-08-16

### DESCRIPTION #################################################################
# This script processes data from the Cancer Cell Line Encyclopedia (CCLE) and
# multiple tissue datasets (TCGA-BRCA, METABRIC, I-SPY2, MATADOR, ICGC BRCA-EU)
# to analyze and visualize the number of outlier genes per patient. It generates
# a density plot comparing the distribution of outlier genes in CCLE and tissue samples.

### PREAMBLE ####################################################################
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-10-08_Figure1_2_3_4_min_input.rda'));

# Outlier patient number
outlier.patient.tag.sum.05 <- apply(ccle.sample.outlier.status, 2, sum);
outlier.patient.sum <- apply(ccle.sample.outlier.status, 1, sum);
outlier.patient.sum <- data.frame(table(outlier.patient.sum));

colnames(outlier.patient.sum) <- c('outlier.num', 'number');

# 1. TCGA-BRCA
outlier.patient.tag.01.brca.gene.per.patient.sum <- apply(outlier.patient.tag.01.brca, 2, sum);

# 2. METABRIC
outlier.patient.tag.01.meta.gene.per.patient.sum <- apply(outlier.patient.tag.01.meta, 2, sum);

# 3. ISPY
outlier.patient.tag.01.ispy.gene.per.patient.sum <- apply(outlier.patient.tag.01.ispy, 2, sum);

# 4. METADOR
outlier.patient.tag.01.metador.gene.per.patient.sum <- apply(outlier.patient.tag.01.metador, 2, sum);

# 5. ICGC BRCA-EU
outlier.patient.tag.01.icgc.gene.per.patient.sum <- apply(outlier.patient.tag.01.icgc, 2, sum);

# Number of outlier genes per patient
outlier.patient.number.violin.tissue <- data.frame(
    out = c(
        outlier.patient.tag.01.brca.gene.per.patient.sum,
        outlier.patient.tag.01.meta.gene.per.patient.sum,
        outlier.patient.tag.01.ispy.gene.per.patient.sum,
        outlier.patient.tag.01.metador.gene.per.patient.sum,
        outlier.patient.tag.01.icgc.gene.per.patient.sum
        ),
    strip = 'tissue'
    );

outlier.patient.tag.sum.05.data <- data.frame(
    out = outlier.patient.tag.sum.05,
    strip = rep('outliers', length(outlier.patient.tag.sum.05))
    );
outlier.patient.number.violin.ccle.tissue <- rbind(outlier.patient.tag.sum.05.data, outlier.patient.number.violin.tissue);
outlier.patient.number.violin.ccle.tissue.log <- outlier.patient.number.violin.ccle.tissue;
outlier.patient.number.violin.ccle.tissue.log$out <- log2(outlier.patient.number.violin.ccle.tissue.log$out + 1);

data <- list(
    x = log2(outlier.patient.tag.sum.05.data$out + 1),
    y = log2(outlier.patient.number.violin.tissue$out + 1)
    );

ccle.col <- grDevices::adjustcolor(c('mediumpurple'), alpha.f = 0.7);

outlier.number.density <- create.densityplot(
    x = data,
    lty = c('solid', 'solid'),
    main = expression('Number of outlier genes per sample'),
    main.cex = 1.5,
    ylab.label = expression('Density'),
    xlab.label = expression('Number of outlier genes'),
    xaxis.lab = c(expression('2'^'0'), expression('2'^'2'), expression('2'^'4'), expression('2'^'6'), expression('2'^'8'), expression('2'^'10')),
    xat = c(0, 2, 4, 6, 8, 10),
    xlimits = c(-2, 10),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    width = 10,
    height = 10,
    lwd = 3,
    col = c(grDevices::adjustcolor(c('navy'), alpha.f = 1), grDevices::adjustcolor(c(ccle.col), alpha.f = 1)),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = c(grDevices::adjustcolor(c('navy'), alpha.f = 1), grDevices::adjustcolor(c(ccle.col), alpha.f = 1)),
                        pch = 21,
                        cex = 2,
                        fill = c(grDevices::adjustcolor(c('navy'), alpha.f = 1), grDevices::adjustcolor(c(ccle.col), alpha.f = 1))
                        ),
                    text = list(
                        lab = c('CCLE (n = 45)', 'All patients (n = 4592)')
                        ),
                    padding.text = c(0, 5, 0),
                    cex = 1
                    )
                ),
            x = 0.45,
            y = 0.97,
            draw = FALSE
            )
        )
    );


save.outlier.figure(
    outlier.number.density,
    c('Figure4b', 'CCLE', 'outlier', 'density'),
    width = 5,
    height = 5
    );

save.session.profile(file.path('output', 'Figure4b.txt'));
