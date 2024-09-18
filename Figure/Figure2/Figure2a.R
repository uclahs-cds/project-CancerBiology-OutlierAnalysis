### HISTORY ######################################################################
# This script performs CNA analysis for METABRIC/TCGA-BRCA/ICGC BRCA-EU data.
# Date: 2024-08-23

### DESCRIPTION ##################################################################
# This script conducts Copy Number Aberration (CNA) analysis across three breast 
# cancer datasets: METABRIC, TCGA-BRCA, and ICGC BRCA-EU. It processes CNA data 
# for outlier and non-outlier samples, performs meta-analysis, and creates 
# visualizations to represent the results.

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-08-23_Figure2a-d.rda'));

### 1. METABRIC
# Convert lists to numeric vectors and calculate medians for meta data
meta.outlier.cnv <- as.numeric(
    unlist(meta.outlier.sample.cnv.new.gis)
    ); # Outlier sample CNA values
meta.non.outlier.cnv <- as.numeric(
    unlist(meta.non.outlier.sample.cnv.new.gis)
    ); # Non-outlier sample CNA values

# Combine outlier and non-outlier values into a single data frame
meta.gis.value <- data.frame(new.gis.value = c(meta.non.outlier.cnv, meta.outlier.cnv));
meta.gis.box <- data.frame(
    cbind(
        meta.gis.value$new.gis.value,
        c(
            rep('non', length(meta.non.outlier.cnv)),
            rep('out', length(meta.outlier.cnv))
            )
        )
    );
colnames(meta.gis.box) <- c('new.gis.value', 'status');
meta.gis.box[, 1] <- as.numeric(meta.gis.box[, 1]);

# Create frequency tables for the CNV values
meta.gis.table.df <- data.frame(table(meta.gis.box));
meta.gis.table <- data.frame(table(meta.gis.box));

# Calculate the relative frequencies
meta.gis.table$Freq[meta.gis.table$status == 'non'] <- meta.gis.table$Freq[meta.gis.table$status == 'non'] / sum(meta.gis.table$Freq[meta.gis.table$status == 'non']);
meta.gis.table$Freq[meta.gis.table$status == 'out'] <- meta.gis.table$Freq[meta.gis.table$status == 'out'] / sum(meta.gis.table$Freq[meta.gis.table$status == 'out']);

# Add sample sizes and calculate standard errors
meta.gis.table.num <- cbind(meta.gis.table, sample = meta.gis.table.df$Freq);
meta.gis.table.num$se <- sqrt(meta.gis.table.num$Freq * (1 - meta.gis.table.num$Freq) / meta.gis.table.num$sample);



### 2. TCGA-BRCA
brca.outlier.cnv <- as.numeric(
    unlist(brca.outlier.sample.cnv.new.gis)
    );
brca.non.outlier.cnv <- as.numeric(
    unlist(brca.non.outlier.sample.cnv.new.gis)
    );

brca.gis.value <- data.frame(new.gis.value = c(brca.non.outlier.cnv, brca.outlier.cnv));
brca.gis.box <- data.frame(
    cbind(
        brca.gis.value$new.gis.value,
        c(
            rep('non', length(brca.non.outlier.cnv)),
            rep('out', length(brca.outlier.cnv))
            )
        )
    );
colnames(brca.gis.box) <- c('new.gis.value', 'status');
brca.gis.box[, 1] <- as.numeric(brca.gis.box[, 1]);

brca.gis.table.df <- data.frame(table(brca.gis.box));
brca.gis.table <- data.frame(table(brca.gis.box));

brca.gis.table$Freq[brca.gis.table$status == 'non'] <- brca.gis.table$Freq[brca.gis.table$status == 'non'] / sum(brca.gis.table$Freq[brca.gis.table$status == 'non']);
brca.gis.table$Freq[brca.gis.table$status == 'out'] <- brca.gis.table$Freq[brca.gis.table$status == 'out'] / sum(brca.gis.table$Freq[brca.gis.table$status == 'out']);

brca.gis.table.num <- cbind(brca.gis.table, sample = brca.gis.table.df$Freq);
brca.gis.table.num$se <- sqrt(brca.gis.table.num$Freq * (1 - brca.gis.table.num$Freq) / brca.gis.table.num$sample);



### 3. ICGC BRCA-EU
icgc.outlier.cnv <- as.numeric(
    unlist(icgc.outlier.sample.cnv.new.gis)
    );
icgc.non.outlier.cnv <- as.numeric(
    unlist(icgc.non.outlier.sample.cnv.new.gis)
    );

icgc.gis.value <- data.frame(new.gis.value = c(icgc.non.outlier.cnv, icgc.outlier.cnv));
icgc.gis.box <- data.frame(
    cbind(
        icgc.gis.value$new.gis.value,
        c(
            rep('non', length(icgc.non.outlier.cnv)),
            rep('out', length(icgc.outlier.cnv))
            )
        )
    );

colnames(icgc.gis.box) <- c('new.gis.value', 'status');
icgc.gis.box[, 1] <- as.numeric(icgc.gis.box[, 1]);
icgc.gis.table.df <- data.frame(table(icgc.gis.box));
icgc.gis.table <- data.frame(table(icgc.gis.box));

icgc.gis.table$Freq[icgc.gis.table$status == 'non'] <- icgc.gis.table$Freq[icgc.gis.table$status == 'non'] / sum(icgc.gis.table$Freq[icgc.gis.table$status == 'non']);
icgc.gis.table$Freq[icgc.gis.table$status == 'out'] <- icgc.gis.table$Freq[icgc.gis.table$status == 'out'] / sum(icgc.gis.table$Freq[icgc.gis.table$status == 'out']);

icgc.gis.table.num <- cbind(icgc.gis.table, sample = icgc.gis.table.df$Freq);
icgc.gis.table.num$se <- sqrt(icgc.gis.table.num$Freq * (1 - icgc.gis.table.num$Freq) / icgc.gis.table.num$sample);



### META ANALYSIS

library(metafor);

# Initialize vectors to store meta-analysis results
metafor.cnv.odd.all <- NULL;
metafor.cnv.lower.all <- NULL;
metafor.cnv.upper.all <- NULL;
metafor.cnv.p.all <- NULL;

# Perform meta-analysis for each of the five CNV categories
for (i in 1:5) {
    # Prepare the data for meta-analysis
    data <- data.frame(
        measure = c('OR', 'OR', 'OR'),
        ai = c(
            meta.gis.table.num$sample[meta.gis.table.num$status == 'out'][i],
            brca.gis.table.num$sample[brca.gis.table.num$status == 'out'][i],
            icgc.gis.table.num$sample[icgc.gis.table.num$status == 'out'][i]
            ),
        n1i = c(
            sum(meta.gis.table.num$sample[meta.gis.table.num$status == 'out']),
            sum(brca.gis.table.num$sample[brca.gis.table.num$status == 'out']),
            sum(icgc.gis.table.num$sample[icgc.gis.table.num$status == 'out'])
            ),
        ci = c(
            meta.gis.table.num$sample[meta.gis.table.num$status == 'non'][i],
            brca.gis.table.num$sample[brca.gis.table.num$status == 'non'][i],
            icgc.gis.table.num$sample[icgc.gis.table.num$status == 'non'][i]
            ),
        n2i = c(
            sum(meta.gis.table.num$sample[meta.gis.table.num$status == 'non']),
            sum(brca.gis.table.num$sample[brca.gis.table.num$status == 'non']),
            sum(icgc.gis.table.num$sample[icgc.gis.table.num$status == 'non'])
            )
        );

    # Calculate effect sizes and variances for meta-analysis
    es.data <- escalc(
        measure = 'OR',
        ai = ai,
        n1i = n1i,
        ci = ci,
        n2i = n2i,
        data = data
        );

    # Perform the meta-analysis using a random-effects model (DerSimonian-Laird method)
    meta.res <- rma.uni(yi, vi, data = es.data, method = 'DL');

    # Extract and store the results: odds ratio, confidence intervals, and p-value
    metafor.cnv.odd <- exp(meta.res$beta);
    metafor.cnv.lower <- exp(meta.res$ci.lb);
    metafor.cnv.upper <- exp(meta.res$ci.ub);
    metafor.cnv.p <- meta.res$pval;

    # Append the results to the corresponding vectors
    metafor.cnv.odd.all <- c(metafor.cnv.odd.all, metafor.cnv.odd);
    metafor.cnv.lower.all <- c(metafor.cnv.lower.all, metafor.cnv.lower);
    metafor.cnv.upper.all <- c(metafor.cnv.upper.all, metafor.cnv.upper);
    metafor.cnv.p.all <- c(metafor.cnv.p.all, metafor.cnv.p);
    }

# Prepare data for visualization (stacked bar plot)

# Adjust p-values for multiple testing using the Benjamini-Hochberg method (FDR)
metafor.cnv.p.all.fdr <- p.adjust(metafor.cnv.p.all, method = 'BH');



### PLOTTING RESULTS ############################################################
# - Use weighted mean
# Initialize an empty vector to store the weighted means
# Initialize vectors to store the weighted sums
brca.sum.vector <- sapply(unique(brca.gis.table.num$new.gis.value), function(x) {
    sum(brca.gis.table.num$sample[brca.gis.table.num$new.gis.value == x])
    })
brca.sum.vector <- rep(brca.sum.vector, each = 2);

meta.sum.vector <- sapply(unique(meta.gis.table.num$new.gis.value), function(x) {
    sum(meta.gis.table.num$sample[meta.gis.table.num$new.gis.value == x])
    })
meta.sum.vector <- rep(meta.sum.vector, each = 2);

icgc.sum.vector <- sapply(unique(icgc.gis.table.num$new.gis.value), function(x) {
    sum(icgc.gis.table.num$sample[icgc.gis.table.num$new.gis.value == x])
    })
icgc.sum.vector <- rep(icgc.sum.vector, each = 2);

# Initialize vector to store the weighted means
all.weighted.mean <- numeric(10);

# Compute weighted means for each status and GIS value combination
for (status in c('non', 'out')) {
    for (i in 1:5) {
        gis.value <- c(-2, -1, 0, 1, 2)[i];

        weights <- c(
            sum(brca.gis.table.num$sample[brca.gis.table.num$new.gis.value == gis.value & brca.gis.table.num$status == status]),
            sum(meta.gis.table.num$sample[meta.gis.table.num$new.gis.value == gis.value & meta.gis.table.num$status == status]),
            sum(icgc.gis.table.num$sample[icgc.gis.table.num$new.gis.value == gis.value & icgc.gis.table.num$status == status])
            );

        values <- c(
            brca.gis.table.num$Freq[brca.gis.table.num$new.gis.value == gis.value & brca.gis.table.num$status == status],
            meta.gis.table.num$Freq[meta.gis.table.num$new.gis.value == gis.value & meta.gis.table.num$status == status],
            icgc.gis.table.num$Freq[icgc.gis.table.num$new.gis.value == gis.value & icgc.gis.table.num$status == status]
            );

        # Compute the weighted average and store it
        index <- ifelse(status == 'non', i, i + 5);
        all.weighted.mean[index] <- weighted.mean(values, weights);
        }
    }



# 1. Create non-outlier fraction bar plot
metafor.cnv.odd.all.df.fraction.non <- data.frame(
    group = rev(c('a', 'b', 'c', 'd', 'e')),
    sample = rep('a', 5),
    odd = all.weighted.mean[1:5]
    );

cnv.col.ramp.5.seg <- c('#2166AC', '#90B2D5', 'grey60', '#D88B95', '#B2182B')
cnv.col.ramp.5.bar <- c('#2166AC', '#90B2D5', 'white', '#D88B95', '#B2182B')

non.fraction.bar <- BoutrosLab.plotting.general:::create.barplot(
    formula = group ~ odd,
    data = metafor.cnv.odd.all.df.fraction.non,
    main = NULL,
    plot.horizontal = TRUE,
    col = cnv.col.ramp.5.bar,
    stack = FALSE,
    main.cex = 0,
    border.col = 'black',
    border.lwd = 1,
    ylab.label = NULL,
    xlab.cex = 1.1,
    xlab.label = expression('Fraction'),
    xlimits = c(-0.07, 0.82),
    ylab.cex = 1.3,
    xaxis.cex = 1.1,
    yaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.rectangle = TRUE,
    xleft.rectangle = -5,
    xright.rectangle = 10,
    ybottom.rectangle = c(1.5, 3.5),
    ytop.rectangle = c(2.5, 4.5),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    );
non.fraction.bar;



# 2. outlier fraction
metafor.cnv.odd.all.df.fraction.out <- data.frame(
    group = rev(c('a', 'b', 'c', 'd', 'e')),
    sample = rep('a', 5),
    odd = all.weighted.mean[6:10]
    );

out.fraction.bar <- BoutrosLab.plotting.general:::create.barplot(
    formula = group ~ odd,
    data = metafor.cnv.odd.all.df.fraction.out,
    main = NULL,
    plot.horizontal = TRUE,
    col = cnv.col.ramp.5.bar,
    stack = FALSE,
    main.cex = 0,
    border.col = 'black',
    border.lwd = 1,
    # main.just = 'right',
    ylab.label = NULL,
    xlab.cex = 1.1,
    xlab.label = expression('Fraction'),
    xlimits = c(-0.07, 0.82),
    ylab.cex = 1.3,
    xaxis.cex = 1.1,
    yaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    # xaxis.rot = 45,
    add.rectangle = TRUE,
    xleft.rectangle = -5,
    xright.rectangle = 10,
    ybottom.rectangle = c(1.5, 3.5),
    ytop.rectangle = c(2.5, 4.5),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    );
out.fraction.bar;



# 3. segment plot
metafor.cnv.odd.ci.p.data <- data.frame(cbind(
    p.value = metafor.cnv.p.all,
    odd = metafor.cnv.odd.all,
    ci.min = metafor.cnv.lower.all,
    ci.max = metafor.cnv.upper.all,
    fdr = metafor.cnv.p.all.fdr
    ));

metafor.cnv.odd.ci.p.data$label <- rev(c('a', 'b', 'c', 'd', 'e'));
metafor.cnv.odd.ci.p.data$label <- as.factor(metafor.cnv.odd.ci.p.data$label);

# Change the order from lowest to higest odds from meta analysis
metafor.cnv.odd.ci.p.data.rev <- metafor.cnv.odd.ci.p.data[rev(1:5), ];

metafor.all.segplot <- BoutrosLab.plotting.general::create.segplot(
    formula = label ~ log2(ci.min) + log2(ci.max),
    data = metafor.cnv.odd.ci.p.data.rev,
    # add middle dots
    centers = log2(metafor.cnv.odd.ci.p.data.rev$odd),
    main.cex = 0,
    yaxis.fontface = 1,
    xlab.cex = 1.3,
    xlab.label = expression('Odds Ratio'),
    symbol.cex = 1.2,
    ylab.cex = 0,
    yaxis.cex = 0,
    xaxis.cex = 1,
    xaxis.lab = c(0, 0.25, 0.5, 1, 2, 4, 8, 16),
    xaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    segments.col = rev(cnv.col.ramp.5.seg),
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -3,
    xright.rectangle = 13,
    ybottom.rectangle = seq(1.5, 23.5, 2),
    ytop.rectangle = seq(2.5, 24.5, 2),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
    );
metafor.all.segplot;


# 4. FDR
metafor.cnv.fdr <- data.frame(
    group = rev(c('a', 'b', 'c', 'd', 'e')),
    sample = rep('a', 5),
    odd = metafor.cnv.odd.ci.p.data$fdr
    );

fdr.bar <- BoutrosLab.plotting.general:::create.barplot(
    formula = group ~ -log10(odd),
    data = metafor.cnv.fdr,
    main = NULL,
    plot.horizontal = TRUE,
    col = 'grey20',
    stack = FALSE,
    main.cex = 0,
    border.col = 'black',
    border.lwd = 1,
    ylab.label = NULL,
    xlab.cex = 1.1,
    xlab.label = expression('-log'[10] * '(FDR)'),
    xlimits = c(-1, 13),
    ylab.cex = 1.3,
    xaxis.cex = 1.1,
    yaxis.cex = 0,
    abline.v = 2,
    abline.lty = 3,
    abline.lwd = 1.5,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    # xaxis.rot = 45,
    add.rectangle = TRUE,
    xleft.rectangle = -5,
    xright.rectangle = 15,
    ybottom.rectangle = c(1.5, 3.5),
    ytop.rectangle = c(2.5, 4.5),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    );
fdr.bar;




multi.gene <- create.multipanelplot(
    list(non.fraction.bar, out.fraction.bar, metafor.all.segplot, fdr.bar),
    main.cex = 0,
    main.y = 0.5,
    resolution = 300,
    layout.height = 1,
    layout.width = 4,
    layout.skip = c(FALSE, FALSE, FALSE, FALSE),
    plot.objects.widths = c(1, 1, 1, 1),
    ylab.axis.padding = -5,
    xlab.cex = 0,
    xlab.axis.padding = 1,
    x.spacing = 1,
    right.legend.padding = 0
    );

multi.gene;

save.outlier.figure(
    multi.gene,
    c('Figure2a', 'CNA', 'multipanel'),
    width = 10.4,
    height = 4.5
    );

save.session.profile(file.path('output', 'Figure2a.txt'));
