### FGFR2_ZSCORE_ANALYSIS ########################################################
# This script calculates the z-scores for the expression of the FGFR2 gene across different 
# datasets. It identifies outlier and non-outlier samples and generates bar plots based on 
# the calculated z-scores.
# Date: 2024-08-13

### PREAMBLE ####################################################################
library(BoutrosLab.utilities)
library(BoutrosLab.plotting.general)


i <- 'FGFR2';

meta.i <- fpkm.tumor.symbol.filter.meta.symbol[
    fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% i, 1:(ncol(fpkm.tumor.symbol.filter.meta.symbol) - 1)
    ];
brca.i <- fpkm.tumor.symbol.filter.brca[
    fpkm.tumor.symbol.filter.brca$Symbol %in% i, 1:(ncol(fpkm.tumor.symbol.filter.brca) - 1)
    ];
icgc.i <- fpkm.tumor.symbol.filter.symbol.icgc[
    fpkm.tumor.symbol.filter.symbol.icgc$Symbol %in% i, 1:(ncol(fpkm.tumor.symbol.filter.symbol.icgc) - 1)
    ];


outlier.status.brca <- outlier.patient.tag.01.brca[
    rownames(fpkm.tumor.symbol.filter.brca[
        fpkm.tumor.symbol.filter.brca$Symbol %in% i,]), 
    ];
outlier.status.meta <- outlier.patient.tag.01.meta[
    rownames(fpkm.tumor.symbol.filter.meta.symbol[
        fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% i,]), 
    ];
outlier.status.icgc <- outlier.patient.tag.01.icgc[
    rownames(fpkm.tumor.symbol.filter.symbol.icgc[
    fpkm.tumor.symbol.filter.symbol.icgc$Symbol %in% i,]), 
    ];

# Convert outlier status to numeric and handle NAs
outlier.status.meta <- as.numeric(outlier.status.meta);
outlier.status.meta[is.na(outlier.status.meta)] <- 0;

outlier.status.brca <- as.numeric(outlier.status.brca);
outlier.status.brca[is.na(outlier.status.brca)] <- 0;

outlier.status.icgc <- as.numeric(outlier.status.icgc);
outlier.status.icgc[is.na(outlier.status.icgc)] <- 0;

outlier.status.all.cnv <- c(outlier.status.meta, outlier.status.brca, outlier.status.icgc);

### MEAN AND STANDARD DEVIATION CALCULATION 
outlier.status.brca.mean <- mean(as.numeric(brca.i)[outlier.status.brca == 0]);
outlier.status.meta.mean <- mean(as.numeric(meta.i)[outlier.status.meta == 0]);
outlier.status.icgc.mean <- mean(as.numeric(icgc.i)[outlier.status.icgc == 0]);

outlier.status.brca.sd <- sd(as.numeric(brca.i)[outlier.status.brca == 0]);
outlier.status.meta.sd <- sd(as.numeric(meta.i)[outlier.status.meta == 0]);
outlier.status.icgc.sd <- sd(as.numeric(icgc.i)[outlier.status.icgc == 0]);

### Z-SCORE CALCULATION 
meta.i.z <- (meta.i - outlier.status.meta.mean) / outlier.status.meta.sd;
brca.i.z <- (brca.i - outlier.status.brca.mean) / outlier.status.brca.sd;
icgc.i.z <- (icgc.i - outlier.status.icgc.mean) / outlier.status.icgc.sd;

three.i.z <- c(as.numeric(meta.i.z), as.numeric(brca.i.z), as.numeric(icgc.i.z));

i.fpkm.merge.data.order.out.three <- three.i.z[outlier.status.all.cnv %in% 1];
i.fpkm.merge.data.order.non.three <- three.i.z[outlier.status.all.cnv %in% 0];

### QUANTILE CALCULATION AND PLOTTING 
unequal.quan <- rev(seq(0, 0.9, 0.1));
i.fpkm.merge.data.order.group.three <- quantile(i.fpkm.merge.data.order.non.three, p = unequal.quan);
i.fpkm.merge.data.order.group.all.three <- c(mean(i.fpkm.merge.data.order.out.three), i.fpkm.merge.data.order.group.three);
i.fpkm.merge.data.order.group.all.three <- c(0, i.fpkm.merge.data.order.group.three);
i.fpkm.merge.data.order.group.all.df.three <- data.frame(
    gene = i.fpkm.merge.data.order.group.all.three,
    sample = LETTERS[1:length(i.fpkm.merge.data.order.group.all.three)]
    );


i.fpkm.merge.data.order.quan.unequal.three <- create.barplot(
    formula = gene ~ sample,
    main = as.expression(substitute(paste(var), list(var = i))),
    ylab.label = expression('z-score'),
    xlab.label = NULL,
    xaxis.lab = rev(c('Outlier samples', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100')),
    xaxis.rot = 90,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.2,
    main.cex = 1.4,
    col = c('darkred', rep('grey10', 10)),
    border.col = 'black',
    border.lwd = 0.25,
    xaxis.fontface = 1, 
    yaxis.fontface = 1, 
    data = i.fpkm.merge.data.order.group.all.df.three,
    sample.order = 'decreasing', 
    ylimits = c(-12, 140)
    );

# Add points to plot
your.points <- xyplot(i.fpkm.merge.data.order.out.three ~ c(0.9, 1, 1, 1.1, 1), 
                      pch = 23, col = "black", fill = 'red2', cex = 1.3);
i.fpkm.merge.data.order.quan.unequal.three.dot <- i.fpkm.merge.data.order.quan.unequal.three + as.layer(your.points);


# Save the bar plot as a PDF
pdf(
    file = generate.filename(
        i, 
        'barplot', 
        'pdf'
        ), 
    width = 4.5, 
    height = 5
    );
i.fpkm.merge.data.order.quan.unequal.three.dot;
dev.off();

# Save the bar plot as a PNG
png(
    file = generate.filename(
        i, 
        'barplot', 
        'png'
        ), 
    width = 4.5, 
    height = 5,
    unit = 'in', 
    res = 1200
    );
i.fpkm.merge.data.order.quan.unequal.three.dot;
dev.off();


