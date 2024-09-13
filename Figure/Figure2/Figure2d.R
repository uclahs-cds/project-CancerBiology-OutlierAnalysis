### ECDNA ANALYSIS #########################################################
# This script analyzes RNA abundance of genes on chromosome 6 in an outlier patient 
# with ecDNA. It calculates robust z-scores for the genes and visualizes the results 
# using a bar plot. 
# Date: 2024-08-13

### Load Required Libraries ######################################################
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));


brca.cnv.chr.new.gis.fpkm.order.match.chr6 <- brca.cnv.chr.new.gis.fpkm.order.match[brca.cnv.chr.new.gis.fpkm.order.match.chr$chromosome == 'chr6',];
brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location <- brca.cnv.chr.new.gis.fpkm.order.match.chr[brca.cnv.chr.new.gis.fpkm.order.match.chr$chromosome %in% 'chr6', ];

# Filter FPKM data for genes on chromosome 6
fpkm.chr6 <- fpkm.tumor.symbol.filter.brca[
    match(brca.cnv.chr.new.gis.fpkm.order.match.chr6$Hugo_Symbol, fpkm.tumor.symbol.filter.brca$Symbol), 
    patient.part.brca
    ];

# Calculate robust z-scores for each gene
fpkm.chr6.median.z <- NULL;
for (i in 1:nrow(fpkm.chr6)) {
    gene.median <- median(as.numeric(fpkm.chr6[i,]));
    gene.mad <- mad(as.numeric(fpkm.chr6[i,]));
    robust.z <- (as.numeric(fpkm.chr6[i,]) - gene.median) / gene.mad;
    fpkm.chr6.median.z <- rbind(fpkm.chr6.median.z, robust.z);
    }

# Set row and column names for the z-score matrix
rownames(fpkm.chr6.median.z) <- rownames(fpkm.chr6);
colnames(fpkm.chr6.median.z) <- colnames(fpkm.chr6);


ecdna.patient <- 'TCGA.A2.A3XX.01A';

# Prepare data for plotting
fpkm.chr6.median.z.bar <- data.frame(
    cbind(fpkm.chr6.median.z[, ecdna.patient, drop = FALSE], sample = c(1:nrow(fpkm.chr6.median.z)))
    );



# Set bar colors, highlighting outlier genes in red
bar.colours <- rep('black', nrow(fpkm.chr6.median.z.bar));
outlier.gene.indices <- which(
    brca.cnv.chr.new.gis.fpkm.order.match.chr6$Hugo_Symbol %in% 
        ecDNA.outlier.patient.location$gene_name[
            ecDNA.outlier.patient.location$sample_barcode %in% substr(ecdna.patient, 1, 15)
            ]
    );
bar.colours[outlier.gene.indices] <- 'red';

# Define ecDNA region start and end genes
chr6.ecdna.start <- as.numeric(ecDNA.outlier.patient.location[
    ecDNA.outlier.patient.location$sample_barcode %in% substr(ecdna.patient, 1, 15),]$start[1]
    );
chr6.ecdna.end <- as.numeric(ecDNA.outlier.patient.location[
    ecDNA.outlier.patient.location$sample_barcode %in% substr(ecdna.patient, 1, 15),]$end[1]
    );
chr6.ecdna.start.gene <- which(
    as.numeric(brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location$start) > chr6.ecdna.start
    )[1];
chr6.ecdna.end.gene <- which(
    as.numeric(brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location$end) < chr6.ecdna.end
)[length(which(as.numeric(brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location$end) < chr6.ecdna.end))];


### PLOTTING ####################################################################

ecdna.bar.chr6 <- create.barplot(
    formula = TCGA.A2.A3XX.01A ~ sample,
    data = fpkm.chr6.median.z.bar,
    main = expression('RNA abundance of genes on chr6 - outlier patient with ecDNA'),
    col = bar.colours,
    main.cex = 1.4,
    border.col = bar.colours,
    border.lwd = 1,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Robust z-score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0, 0),
    add.rectangle = TRUE,
    xleft.rectangle = c(chr6.ecdna.start.gene - 1),
    xright.rectangle = c(chr6.ecdna.end.gene + 1),
    ybottom.rectangle = -20,
    ytop.rectangle = 170,
    col.rectangle = "grey",
    alpha.rectangle = 0.5,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = 'red',
                        pch = 22,
                        cex = 1.5,
                        fill = c('red')
                        ),
                    text = list(
                        lab = c('Outlier genes on ecDNA')
                        ),
                    padding.text = 1,
                    cex = 1
                    )
                ),
            x = 0.65,
            y = 0.95
            )
        )
    );


save.outlier.figure(
    ecdna.bar.chr6,
    c('ecdna', 'chr6', 'barplot'),
    width = 8, 
    height = 3.5
    );
