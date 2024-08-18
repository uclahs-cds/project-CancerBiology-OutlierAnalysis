### HISTORY #####################################################################
# This script performs meta-analysis for multiple datasets across chromosomes
# Date: 2024-08-13


library(metafor);

### DATA PREPARATION ############################################################

# 1. Chromosomal enrichment
chr.odd.se.5 <- list();

# Loop through each row of the chromosome odds ratio dataframe
for (i in 1:nrow(p.value.chr.brca.odd.sub.df)) {
    
    # Extract the odds ratios and standard errors for each dataset
    chr.odd <- c(
        ln.odd.brca[i],
        ln.odd.meta[i],
        ln.odd.ispy[i],
        ln.odd.metador[i],
        ln.odd.icgc[i]
        );
    chr.se <- c(
        se.odd.brca[i],
        se.odd.meta[i],
        se.odd.ispy[i],
        se.odd.metador[i],
        se.odd.icgc[i]
        );
    
    # Combine odds ratios and standard errors into a single data frame
    chr.all <- data.frame(cbind(chr.odd, chr.se));
    
    # Store the data frame in the list
    chr.odd.se.5[[i]] <- chr.all;
    
    }

metafor.chr.odd.ci.p.5 <- NULL;

# Loop through each chromosome's data for meta-analysis
for (i in 1:nrow(p.value.chr.brca.odd.sub.df)) {
    
    # Extract the sample data for the current chromosome
    chr.odd.se.sample <- chr.odd.se.5[[i]];
    
    # Remove infinite values from the sample data
    chr.odd.se.sample.inf <- chr.odd.se.sample[
        !is.infinite(chr.odd.se.sample$chr.odd) & !is.infinite(chr.odd.se.sample$chr.se)
        ];
    
    # Perform meta-analysis if there are more than one valid sample
    if (nrow(chr.odd.se.sample.inf) > 1) {
        
        metafor.chr <- rma.uni(
            yi = chr.odd, 
            sei = chr.se, 
            data = chr.odd.se.sample.inf, 
            method = 'DL'
            );
        
        metafor.chr.odd <- exp(metafor.chr$beta);
        metafor.chr.lower <- exp(metafor.chr$ci.lb);
        metafor.chr.upper <- exp(metafor.chr$ci.ub);
        metafor.chr.p <- metafor.chr$pval;
        
        metafor.all <- c(
            metafor.chr.odd, 
            metafor.chr.lower, 
            metafor.chr.upper, 
            metafor.chr.p
            );
        
        } 
    else {
        
        metafor.all <- c(NA, NA, NA, NA);
        
        }
    
    # Append the results to the matrix
    metafor.chr.odd.ci.p.5 <- rbind(metafor.chr.odd.ci.p.5, metafor.all);
    
    }

# Create a data frame to store meta-analysis results with labels
metafor.chr.odd.ci.p.data.5 <- data.frame(
    p.value = metafor.chr.odd.ci.p.5[1:24, 4],
    odd = metafor.chr.odd.ci.p.5[1:24, 1],
    ci.min = metafor.chr.odd.ci.p.5[1:24, 2], 
    ci.max = metafor.chr.odd.ci.p.5[1:24, 3]
    );

metafor.chr.odd.ci.p.data.label.5 <- data.frame(
    metafor.chr.odd.ci.p.data.5,
    labels = as.factor(paste('Chr', chr.name[1:24], sep = ''))
    );




# 2. Exon number/Gene length/GC content/RNA abundance
metafor.smd.all.matrix.5 <- data.frame(cbind(
    p.value = c(smd.metafor.exon.5$pval,
                smd.metafor.length.5$pval,
                smd.metafor.gc.5$pval,
                smd.metafor.rna.5$pval),
    odd = c(exp(smd.metafor.exon.5$beta),
            exp(smd.metafor.length.5$beta),
            exp(smd.metafor.gc.5$beta),
            exp(smd.metafor.rna.5$beta)),
    ci.min = c(exp(smd.metafor.exon.5$ci.lb),
               exp(smd.metafor.length.5$ci.lb),
               exp(smd.metafor.gc.5$ci.lb),
               exp(smd.metafor.rna.5$ci.lb)),
    ci.max = c(exp(smd.metafor.exon.5$ci.ub),
               exp(smd.metafor.length.5$ci.ub),
               exp(smd.metafor.gc.5$ci.ub),
               exp(smd.metafor.rna.5$ci.ub))
    ));

metafor.smd.all.matrix.label.5 <- data.frame(cbind(
    metafor.smd.all.matrix.5,
    labels = as.factor(c('Exon number', 'Gene length', 'GC content', 'RNA abundance'))
    ));

# FDR correction for the SMD matrix
metafor.smd.all.matrix.label.fdr.5 <- p.adjust(metafor.smd.all.matrix.label.5$p.value, method = 'BH');



### PLOTTING RESULTS ############################################################

# Perform FDR correction on the p-values
fdr.metafor.schr <- p.adjust(metafor.chr.odd.ci.p.data.label.rev.5$p.value, method = 'BH');
dot.colours <- vector(length=nrow(metafor.chr.odd.ci.p.data.label.rev.5));
dot.colours <- rep('grey70',nrow(metafor.chr.odd.ci.p.data.label.rev.5));
dot.colours[fdr.metafor.schr.5 < 0.01 & metafor.chr.odd.ci.p.data.label.rev.5$odd < 1] <- 'dodgerblue2';
dot.colours[fdr.metafor.schr.5 < 0.01 & metafor.chr.odd.ci.p.data.label.rev.5$odd > 1] <- 'red';

metafor.all.segplot.multi.5 <- BoutrosLab.plotting.general::create.segplot(
    formula = labels ~ log2(ci.min) + log2(ci.max),
    data = metafor.chr.odd.ci.p.data.label.rev.5,
    centers = log2(metafor.chr.odd.ci.p.data.label.rev.5$odd),
    main.cex = 0,
    ylab.label = expression('Chromosome'),
    yaxis.lab = metafor.chr.odd.ci.p.data.label.rev.5$labels,
    yaxis.fontface = 1,
    xlab.cex = 0,
    xlab.label = NULL,
    ylab.cex = 1.3,
    yaxis.cex = 1,
    xaxis.cex = 0,
    xlimits = c(-1.1, 1.3),
    xaxis.lab = c('0.5', '1', '2'),
    xat = seq(-1, 1, 1),
    xaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    segments.col = dot.colours,
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -1.7,
    xright.rectangle = 13,
    ybottom.rectangle = seq(1.5, 23.5, 2),
    ytop.rectangle = seq(2.5, 24.5, 2),
    # set rectangle colour
    col.rectangle = "grey",
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
    )
metafor.all.segplot.multi.5;


dot.colours <- vector(length=nrow(metafor.smd.all.matrix.label.5));
dot.colours <- rep('grey70',nrow(metafor.smd.all.matrix.label.5));
dot.colours[metafor.smd.all.matrix.label.fdr.5 < 0.01 & metafor.smd.all.matrix.label.5$odd < 1] <- 'dodgerblue2';
dot.colours[metafor.smd.all.matrix.label.fdr.5 < 0.01 & metafor.smd.all.matrix.label.5$odd > 1] <- 'red';
metafor.smd.all.matrix.label.segplot.multi.5 <- BoutrosLab.plotting.general::create.segplot(
    formula = labels ~ log2(ci.min) + log2(ci.max),
    data = metafor.smd.all.matrix.label.5,
    # add middle dots
    centers = log2(metafor.smd.all.matrix.label.5$odd),
    main.cex = 1.4,
    ylab.label = NULL,
    yaxis.lab = metafor.smd.all.matrix.label.5$labels,
    yaxis.fontface = 1,
    xlab.cex = 1.3,
    xlab.label = expression('Odds Ratio'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1,
    xlimits = c(-1.1, 1.3),
    xaxis.lab = c('0.5', '1', '2'),
    xat = seq(-1, 1, 1),
    xaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    segments.col = dot.colours,
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -1.7,
    xright.rectangle = 13,
    ybottom.rectangle = seq(1.5, 23.5, 2),
    ytop.rectangle = seq(2.5, 24.5, 2),
    # set rectangle colour
    col.rectangle = "grey",
    # set rectangle alpha (transparency)
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
)
metafor.smd.all.matrix.label.segplot.multi.5;



metafor.multi.chr.smd.5 <- create.multipanelplot(
    list(metafor.all.segplot.multi.5, metafor.smd.all.matrix.label.segplot.multi.5),
    resolution = 300,
    layout.height = 2,
    layout.width = 1,
    ylab.cex = 0,
    xlab.cex = 0,
    layout.skip = c(FALSE, FALSE),
    plot.objects.heights = c(5, 2),
    ylab.axis.padding = -7,
    xlab.axis.padding = 1,
    y.spacing = -3.5,
    bottom.padding = -0.5,
    top.padding = -0.5
    );

metafor.multi.chr.smd.5;


# Save the multi plot as a PDF
pdf(
    file = generate.filename(
        'metafor.multi.chr.smd.5', 
        'multipanel', 
        'pdf'
        ), 
    width = 5, 
    height = 6.5
    );
metafor.multi.chr.smd.5;
dev.off();

# Save the multi plot as a PNG
png(
    file = generate.filename(
        'metafor.multi.chr.smd.5', 
        'multipanel', 
        'png'
        ), 
    width = 5, 
    height = 6.5,
    unit = 'in', 
    res = 1200
    );
metafor.multi.chr.smd.5;
dev.off();
