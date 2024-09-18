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


genes <- c('IGF2', 'TMEM30A', 'NRAS', 'IGF2R', 'GAPDH', 'B2M'); # 여섯 개의 유전자

gene.data <- lapply(seq_along(genes), function(gene.index) {
    # Extract data for the specified gene across different datasets
    gene.alias <- genes[gene.index];

    # fpkm.tumor.symbol.filter.XXX: RNA abundance matrix with rows of genes and
    # columns of patients. Row names are Ensembl IDs.

    # fpkm.tumor.symbol.filter.XXX.symbol is the same, except it has one
    # additional `Symbol` column of gene aliases.

    metador.i <- fpkm.tumor.symbol.filter.metador.symbol[
        gene.alias == fpkm.tumor.symbol.filter.metador.symbol$Symbol,
        !('Symbol' == colnames(fpkm.tumor.symbol.filter.metador.symbol))
        ];

    meta.i <- na.omit(
        fpkm.tumor.symbol.filter.meta.symbol[
            gene.alias == fpkm.tumor.symbol.filter.meta.symbol$Symbol,
            !('Symbol' == colnames(fpkm.tumor.symbol.filter.meta.symbol))
            ]
        );

    brca.i <- fpkm.tumor.symbol.filter.brca[
        gene.alias == fpkm.tumor.symbol.filter.brca$Symbol,
        !('Symbol' == colnames(fpkm.tumor.symbol.filter.brca))
        ];

    # ispy an exception in that the row names are
    # the gene aliases.
    ispy.i <- fpkm.tumor.symbol.filter.ispy[
        gene.alias == rownames(fpkm.tumor.symbol.filter.ispy),
        ];

    # fpkm.tumor.symbol.filter.symbol.icgc is also an exception in that the row
    # names are (seemingly arbitrary) integers
    icgc.i <- fpkm.tumor.symbol.filter.symbol.icgc[
        gene.alias == fpkm.tumor.symbol.filter.symbol.icgc$Symbol,
        !('Symbol' == colnames(fpkm.tumor.symbol.filter.symbol.icgc))
        ];

    ### OUTLIER STATUS ##############################################################

    # outlier.patient.tag.01.XXX: Outlier status matrix of XXX dataset with
    # rows of genes and columns of patients. Values are 1 for outlier events
    # and 0 otherwise.
    outlier.status.metador <- as.numeric(outlier.patient.tag.01.metador[rownames(metador.i), ]);
    outlier.status.metador[is.na(outlier.status.metador)] <- 0;

    outlier.status.ispy <- as.numeric(outlier.patient.tag.01.ispy[rownames(ispy.i), ]);
    outlier.status.ispy[is.na(outlier.status.ispy)] <- 0;

    outlier.status.meta <- as.numeric(outlier.patient.tag.01.meta[rownames(meta.i), ]);
    outlier.status.meta[is.na(outlier.status.meta)] <- 0;

    outlier.status.brca <- as.numeric(outlier.patient.tag.01.brca[rownames(brca.i), ]);
    outlier.status.brca[is.na(outlier.status.brca)] <- 0;

    outlier.status.icgc <- as.numeric(outlier.patient.tag.01.icgc[rownames(icgc.i), ]);
    outlier.status.icgc[is.na(outlier.status.icgc)] <- 0;

    ### COMBINING OUTLIER STATUS ###################################################

    # Combine outlier status from all datasets
    outlier.status.all <- c(
        outlier.status.metador,
        outlier.status.ispy,
        outlier.status.meta,
        outlier.status.brca,
        outlier.status.icgc
        );

    ### CALCULATING MEANS AND SD ###################################################

    # Calculate the mean and standard deviation for non-outlier patients in each dataset
    outlier.status.brca.mean <- mean(as.numeric(brca.i)[outlier.status.brca == 0]);
    outlier.status.meta.mean <- mean(as.numeric(meta.i)[outlier.status.meta == 0]);
    outlier.status.ispy.mean <- mean(na.omit(as.numeric(ispy.i)[outlier.status.ispy == 0]));
    outlier.status.metador.mean <- mean(as.numeric(metador.i)[outlier.status.metador == 0]);
    outlier.status.icgc.mean <- mean(as.numeric(icgc.i)[outlier.status.icgc == 0]);

    outlier.status.brca.sd <- sd(as.numeric(brca.i)[outlier.status.brca == 0]);
    outlier.status.meta.sd <- sd(as.numeric(meta.i)[outlier.status.meta == 0]);
    outlier.status.ispy.sd <- sd(na.omit(as.numeric(ispy.i)[outlier.status.ispy == 0]));
    outlier.status.metador.sd <- sd(as.numeric(metador.i)[outlier.status.metador == 0]);
    outlier.status.icgc.sd <- sd(as.numeric(icgc.i)[outlier.status.icgc == 0]);

    ### CALCULATING Z-SCORES #######################################################

    # Calculate the z-scores for each dataset
    metador.i.z <- (metador.i - outlier.status.metador.mean) / outlier.status.metador.sd;
    ispy.i.z <- (ispy.i - outlier.status.ispy.mean) / outlier.status.ispy.sd;
    meta.i.z <- (meta.i - outlier.status.meta.mean) / outlier.status.meta.sd;
    brca.i.z <- (brca.i - outlier.status.brca.mean) / outlier.status.brca.sd;
    icgc.i.z <- (icgc.i - outlier.status.icgc.mean) / outlier.status.icgc.sd;

    # Combine all z-scores into a single vector
    five.i.z <- c(
        as.numeric(metador.i.z),
        as.numeric(ispy.i.z),
        as.numeric(meta.i.z),
        as.numeric(brca.i.z),
        as.numeric(icgc.i.z)
        );

    ### DATA FRAME PREPARATION #####################################################

    list(
        z.score = data.frame(
            sample = gene.alias,
            value = as.numeric(five.i.z),
            order = letters[gene.index]
            ),
        outlier.status = outlier.status.all
        );
    });



### DATA FRAME PREPARATION #####################################################

all.gene.z.scores.ordered <- do.call(rbind, lapply(gene.data, `[[`, 'z.score'));

all.gene.outlier.status <- unlist(lapply(gene.data, `[[`, 'outlier.status'))

# Define point colors, shapes, and sizes based on outlier status
patient.colors <- ifelse(all.gene.outlier.status == 1, 'red2', 'black');
patient.shapes <- ifelse(all.gene.outlier.status == 1, 23, 21);
patient.sizes <- ifelse(all.gene.outlier.status == 1, 0.9, 0.75);

### DATA ANALYSIS ###############################################################

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure1b')));

# Create the strip plot
stripplot.gene.z.scores <- BoutrosLab.plotting.general::create.stripplot(
    formula = as.numeric(value) ~ order,
    data = all.gene.z.scores.ordered,
    xaxis.cex = 1.1,
    yaxis.cex = 1,
    xaxis.lab = c('IGF2', 'TMEM30A', 'NRAS', 'IGF2R', 'GAPDH', 'B2M'),
    yat = seq(0, 200, 20),
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    cex = patient.sizes,
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 3.5, 5.5),
    xright.rectangle = c(2.5, 4.5, 6.5),
    ybottom.rectangle = -50,
    ytop.rectangle = 1000,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    ylab.label = expression('z-score'),
    xlab.label = NULL,
    col = patient.colors,
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
    pch = patient.shapes,
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
