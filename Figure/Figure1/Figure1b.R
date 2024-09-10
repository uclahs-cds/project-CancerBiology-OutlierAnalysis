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

### DATA PREPARATION ############################################################

genes <- c('IGF2', 'TMEM30A', 'NRAS', 'IGF2R', 'GAPDH', 'B2M'); # 여섯 개의 유전자

for (i in genes) {
    
    # Extract data for the specified gene across different datasets
    metador.i <- fpkm.tumor.symbol.filter.metador.symbol[
        fpkm.tumor.symbol.filter.metador.symbol$Symbol %in% i, 
        patient.part.metador
        ];
    
    ispy.i <- fpkm.tumor.symbol.filter.ispy[
        rownames(fpkm.tumor.symbol.filter.ispy) %in% i, 
        patient.part.ispy
        ];
    
    meta.i <- na.omit(
        fpkm.tumor.symbol.filter.meta.symbol[
            fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% i, 
            patient.part.meta
            ]
        );
        
    brca.i <- fpkm.tumor.symbol.filter.brca[
        fpkm.tumor.symbol.filter.brca$Symbol %in% i, 
        patient.part.brca
        ];
    
    icgc.i <- fpkm.tumor.symbol.filter.icgc[
        fpkm.tumor.symbol.filter.symbol.icgc$Symbol %in% i, 
        ];
    
    ### OUTLIER STATUS ##############################################################
    
    # Calculate the outlier status for each dataset
    outlier.status.brca <- outlier.patient.tag.01.brca[
        rownames(fpkm.tumor.symbol.filter.brca[
            fpkm.tumor.symbol.filter.brca$Symbol %in% i, 
            ]), 
        ];
        
    outlier.status.meta <- outlier.patient.tag.01.meta[
        rownames(fpkm.tumor.symbol.filter.meta.symbol[
            fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% i, 
            patient.part.meta
            ]), 
        ];
        
    outlier.status.ispy <- outlier.patient.tag.01.ispy[
        gsub(".*_", "", rownames(outlier.patient.tag.01.ispy)) %in% i, 
        ];
    
    outlier.status.metador <- outlier.patient.tag.01.metador[
        gsub(".*_", "", rownames(outlier.patient.tag.01.metador)) %in% i, 
        ];
    
    outlier.status.icgc <- outlier.patient.tag.01.icgc[
        fpkm.tumor.symbol.filter.symbol.icgc[
            rownames(outlier.patient.tag.01.icgc), 
            ]$Symbol %in% i, 
        ];
    
    # Convert outlier statuses to numeric and handle missing values
    outlier.status.metador <- as.numeric(outlier.status.metador);
    outlier.status.metador[is.na(outlier.status.metador)] <- 0;
    
    outlier.status.ispy <- as.numeric(outlier.status.ispy);
    outlier.status.ispy[is.na(outlier.status.ispy)] <- 0;
    
    outlier.status.meta <- as.numeric(outlier.status.meta);
    outlier.status.meta[is.na(outlier.status.meta)] <- 0;
    
    outlier.status.brca <- as.numeric(outlier.status.brca);
    outlier.status.brca[is.na(outlier.status.brca)] <- 0;
    
    outlier.status.icgc <- as.numeric(outlier.status.icgc);
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
    
    # Store the z-scores and outlier status in variables
    gene.z.score <- paste(i, '.z.score', sep = '');
    assign(gene.z.score, five.i.z);
    
    gene.outlier.status <- paste(i, '.outlier.status', sep = '');
    assign(gene.outlier.status, outlier.status.all);
    
    ### DATA FRAME PREPARATION #####################################################
    
    # Create data frames for each gene's z-scores
    assign(paste(i, ".i.frame.z", sep = ''), data.frame(
        sample = rep(i, length(get(gene.z.score))), 
        value = as.numeric(get(gene.z.score))
        )
    );
    
    }



### DATA FRAME PREPARATION #####################################################

# Create data frames for each gene's z-scores
# Create data frames for each gene's z-scores
IGF2.i.frame.z <- data.frame(
    sample = rep('IGF2', length(IGF2.z.score)), 
    value = as.numeric(IGF2.z.score)
    );

IGF2R.i.frame.z <- data.frame(
    sample = rep('IGF2R', length(IGF2R.z.score)), 
    value = as.numeric(IGF2R.z.score)
    );

TMEM30A.i.frame.z <- data.frame(
    sample = rep('TMEM30A', length(TMEM30A.z.score)), 
    value = as.numeric(TMEM30A.z.score)
    );

GAPDH.i.frame.z <- data.frame(
    sample = rep('GAPDH', length(GAPDH.z.score)), 
    value = as.numeric(GAPDH.z.score)
    );

B2M.i.frame.z <- data.frame(
    sample = rep('B2M', length(B2M.z.score)), 
    value = as.numeric(B2M.z.score)
    );

NRAS.i.frame.z <- data.frame(
    sample = rep('NRAS', length(NRAS.z.score)), 
    value = as.numeric(NRAS.z.score)
    );



### COMBINING DATA #############################################################

# Combine the data for all genes into a single data frame
all.gene.z.scores <- rbind(
    IGF2.i.frame.z,
    TMEM30A.i.frame.z,
    NRAS.i.frame.z,
    IGF2R.i.frame.z,
    GAPDH.i.frame.z,
    B2M.i.frame.z
    );

# Add an order column to distinguish between different genes
all.gene.z.scores.ordered <- cbind(
    all.gene.z.scores,
    order = c(
        rep('a', nrow(IGF2.i.frame.z)),
        rep('b', nrow(TMEM30A.i.frame.z)),
        rep('c', nrow(NRAS.i.frame.z)),
        rep('d', nrow(IGF2R.i.frame.z)),
        rep('e', nrow(GAPDH.i.frame.z)),
        rep('f', nrow(B2M.i.frame.z))
        )
    );

# Combine outlier status from all genes
all.gene.outlier.status <- c(
    IGF2.outlier.status,
    TMEM30A.outlier.status,
    NRAS.outlier.status,
    IGF2R.outlier.status,
    GAPDH.outlier.status,
    B2M.outlier.status
    );

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

### OUTPUT ######################################################################

# Save the plot as a PNG
output.directory <- get0('output.directory', ifnotfound = 'figures');

write.plot(
    trellis.object = stripplot.gene.z.scores,
    filename = file.path(output.directory, 'Figure_1_b.png'),
    width = 4.7,
    height = 5.3,
    size.units = 'in',
    resolution = 1200
);
