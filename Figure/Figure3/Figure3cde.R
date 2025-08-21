### HISTORY ######################################################################
# This script performs pooling analysis to examine the XEG enrichment in
# multiple breast cancer subtypes.
# Date: 2025-08-14

### DESCRIPTION ##################################################################
# This script analyzes XEG enrichment across different breast cancer subtypes
# using data from all nine datasets. It performs
# Fisher's exact tests, calculates odds ratios. The script generates visualizations
# of the results, including volcano plot for each subtype, contingency table of 
# example genes, and heatmap.

### PREAMBLE #####################################################################
# Load required libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
library(ConsensusClusterPlus);
library(cluster);
library(factoextra);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################


# XEG enrichment in each subtype
#   - volcano plot
# In all 9 dataset
outlier.patient.all.nine.01.sum.patient <- apply(outlier.patient.all.nine.01, 1, function(x) { sum(na.omit(x)); });
frequent.gene.nine.10 <- outlier.patient.all.nine.01.sum.patient[outlier.patient.all.nine.01.sum.patient > 9];
frequent.gene.nine.10 <- names(outlier.patient.all.nine.01.sum.patient);
colnames(subtype.5.total.outlier.num.meta) <- colnames(subtype.total.outlier.num.1.brca);

all.subtype.info <- rbind(
    subtype.total.outlier.num.1.brca,
    subtype.5.total.outlier.num.meta,
    subtype.total.outlier.num.1.ispy,
    subtype.total.outlier.num.1.cheng,
    subtype.total.outlier.num.1.sjostrom,
    subtype.total.outlier.num.1.matador,
    subtype.total.outlier.num.1.icgc,
    subtype.total.outlier.num.1.kao,
    subtype.total.outlier.num.1.hatzis
    );

os.group.combine.gene.subtype.enrich.list <- list();
os.group.combine.gene.subtype.enrich.basal <- NULL;
os.group.combine.gene.subtype.enrich.her2 <- NULL;
os.group.combine.gene.subtype.enrich.luma <- NULL;
os.group.combine.gene.subtype.enrich.lumb <- NULL;
os.group.combine.gene.subtype.enrich.normal <- NULL;

frequent.gene.nine.10 <- names(outlier.patient.all.nine.01.sum.patient);

for (i in frequent.gene.nine.10) {
    outlier.patient.all.nine.01.gene <- outlier.patient.all.nine.01[, match(rownames(all.subtype.info), colnames(outlier.patient.all.nine.01))][i, ];
    
    # Create outlier flag vector
    os.group.combine.gene <- all.subtype.info;
    os.group.combine.gene$outlier_flag <- as.numeric(unlist(outlier.patient.all.nine.01.gene));
    os.group.combine.gene <- na.omit(os.group.combine.gene);
    os.group.combine.gene$outlier_flag <- ifelse(os.group.combine.gene$outlier_flag == 1, TRUE, FALSE);
    
    # Define unique subtype levels
    pam50_levels <- c(1, 2, 3, 4, 5);
    
    # Initialize results dataframe
    results <- data.frame(
        subtype = pam50_levels,
        OR = NA,
        CI_low = NA,
        CI_high = NA,
        p_value = NA,
        subtype.out = NA,
        stringsAsFactors = FALSE
        );
    
    # Iterate for each subtype
    for (k in seq_along(pam50_levels)) {
        subtype_i <- pam50_levels[k];
        
        # Build 2x2 contingency table
        a <- sum(os.group.combine.gene$subtype == subtype_i & os.group.combine.gene$outlier_flag == TRUE);
        b <- sum(os.group.combine.gene$subtype != subtype_i & os.group.combine.gene$outlier_flag == TRUE);
        c <- sum(os.group.combine.gene$subtype == subtype_i & os.group.combine.gene$outlier_flag == FALSE);
        d <- sum(os.group.combine.gene$subtype != subtype_i & os.group.combine.gene$outlier_flag == FALSE);
        
        contingency <- matrix(c(a + 1, b + 1, c + 1, d + 1), nrow = 2, byrow = TRUE);
        
        # Fisher's exact test
        fisher_res <- fisher.test(contingency);
        
        # Save results
        results$OR[k] <- ifelse(is.finite(fisher_res$estimate), fisher_res$estimate, Inf);
        results$CI_low[k] <- fisher_res$conf.int[1];
        results$CI_high[k] <- fisher_res$conf.int[2];
        results$p_value[k] <- fisher_res$p.value;
        results$subtype.out[k] <- a;
        };
    
    os.group.combine.gene.subtype.enrich.list[i] <- results;
    
    os.group.combine.gene.subtype.enrich.basal <- rbind(os.group.combine.gene.subtype.enrich.basal, results[1, ]);
    os.group.combine.gene.subtype.enrich.her2 <- rbind(os.group.combine.gene.subtype.enrich.her2, results[2, ]);
    os.group.combine.gene.subtype.enrich.luma <- rbind(os.group.combine.gene.subtype.enrich.luma, results[3, ]);
    os.group.combine.gene.subtype.enrich.lumb <- rbind(os.group.combine.gene.subtype.enrich.lumb, results[4, ]);
    os.group.combine.gene.subtype.enrich.normal <- rbind(os.group.combine.gene.subtype.enrich.normal, results[5, ]);
    };


os.group.combine.gene.subtype.enrich.basal$fdr <- p.adjust(os.group.combine.gene.subtype.enrich.basal$p_value, method = 'BH');
os.group.combine.gene.subtype.enrich.her2$fdr <- p.adjust(os.group.combine.gene.subtype.enrich.her2$p_value, method = 'BH');
os.group.combine.gene.subtype.enrich.luma$fdr <- p.adjust(os.group.combine.gene.subtype.enrich.luma$p_value, method = 'BH');
os.group.combine.gene.subtype.enrich.lumb$fdr <- p.adjust(os.group.combine.gene.subtype.enrich.lumb$p_value, method = 'BH');
os.group.combine.gene.subtype.enrich.normal$fdr <- p.adjust(os.group.combine.gene.subtype.enrich.normal$p_value, method = 'BH');

rownames(os.group.combine.gene.subtype.enrich.basal) <- frequent.gene.nine.10;
rownames(os.group.combine.gene.subtype.enrich.her2) <- frequent.gene.nine.10;
rownames(os.group.combine.gene.subtype.enrich.luma) <- frequent.gene.nine.10;
rownames(os.group.combine.gene.subtype.enrich.lumb) <- frequent.gene.nine.10;
rownames(os.group.combine.gene.subtype.enrich.normal) <- frequent.gene.nine.10;

os.group.combine.gene.subtype.enrich.basal$out.num <- outlier.patient.all.nine.01.sum.patient;
os.group.combine.gene.subtype.enrich.her2$out.num <- outlier.patient.all.nine.01.sum.patient;
os.group.combine.gene.subtype.enrich.luma$out.num <- outlier.patient.all.nine.01.sum.patient;
os.group.combine.gene.subtype.enrich.lumb$out.num <- outlier.patient.all.nine.01.sum.patient;
os.group.combine.gene.subtype.enrich.normal$out.num <- outlier.patient.all.nine.01.sum.patient;


# Cut-off
# - > 0.1% total
# - > 0.5% at least one subtype
subtype.each.out.num <- data.frame(
    basal = os.group.combine.gene.subtype.enrich.basal$subtype.out,
    her2 = os.group.combine.gene.subtype.enrich.her2$subtype.out,
    luma = os.group.combine.gene.subtype.enrich.luma$subtype.out,
    lumb = os.group.combine.gene.subtype.enrich.lumb$subtype.out,
    normal = os.group.combine.gene.subtype.enrich.normal$subtype.out
    );

subtype.each.out.num$basal.status <- ifelse(subtype.each.out.num$basal > sum(na.omit(all.subtype.info$subtype == 1)) * 0.005, TRUE, FALSE);
subtype.each.out.num$her2.status <- ifelse(subtype.each.out.num$her2 > sum(na.omit(all.subtype.info$subtype == 2)) * 0.005, TRUE, FALSE);
subtype.each.out.num$luma.status <- ifelse(subtype.each.out.num$luma > sum(na.omit(all.subtype.info$subtype == 3)) * 0.005, TRUE, FALSE);
subtype.each.out.num$lumb.status <- ifelse(subtype.each.out.num$lumb > sum(na.omit(all.subtype.info$subtype == 4)) * 0.005, TRUE, FALSE);
subtype.each.out.num$normal.status <- ifelse(subtype.each.out.num$normal > sum(na.omit(all.subtype.info$subtype == 5)) * 0.005, TRUE, FALSE);

subtype.each.out.num$all.status <- 
    subtype.each.out.num$basal.status +
    subtype.each.out.num$her2.status +
    subtype.each.out.num$luma.status +
    subtype.each.out.num$lumb.status +
    subtype.each.out.num$normal.status;

subtype.each.out.num$sum.patient <- outlier.patient.all.nine.01.sum.patient;

subtype.each.out.num.selected <- subtype.each.out.num[
    subtype.each.out.num$all.status > 0 &
    subtype.each.out.num$sum.patient > ncol(outlier.patient.all.nine.01) * 0.001, 
    ];

# Use different cut-off based on per-subtype selection
os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01 <- os.group.combine.gene.subtype.enrich.basal[
    subtype.each.out.num$all.status > 0 &
    subtype.each.out.num$sum.patient > ncol(outlier.patient.all.nine.01) * 0.001, 
    ];
os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$cut.fdr <- p.adjust(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$p_value, method = 'BH');

os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01 <- os.group.combine.gene.subtype.enrich.her2[
    subtype.each.out.num$all.status > 0 &
    subtype.each.out.num$sum.patient > ncol(outlier.patient.all.nine.01) * 0.001, 
    ];
os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$cut.fdr <- p.adjust(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$p_value, method = 'BH');

os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01 <- os.group.combine.gene.subtype.enrich.luma[
    subtype.each.out.num$all.status > 0 &
    subtype.each.out.num$sum.patient > ncol(outlier.patient.all.nine.01) * 0.001, 
    ];
os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$cut.fdr <- p.adjust(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$p_value, method = 'BH');

os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01 <- os.group.combine.gene.subtype.enrich.lumb[
    subtype.each.out.num$all.status > 0 &
    subtype.each.out.num$sum.patient > ncol(outlier.patient.all.nine.01) * 0.001, 
    ];
os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$cut.fdr <- p.adjust(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$p_value, method = 'BH');

os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01 <- os.group.combine.gene.subtype.enrich.normal[
    subtype.each.out.num$all.status > 0 &
    subtype.each.out.num$sum.patient > ncol(outlier.patient.all.nine.01) * 0.001, 
    ];
os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$cut.fdr <- p.adjust(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$p_value, method = 'BH');


# 1. Basal
dot.colours <- vector(length = nrow(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01));
dot.colours <- rep('grey80', nrow(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01));
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$OR) < 1 & as.numeric(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$cut.fdr) < 0.05] <- 'dodgerblue3';
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$OR) > 1 & as.numeric(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$cut.fdr) < 0.05] <- 'red2';

# Optional labeling block (kept commented)
# interesting.points <- os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$OR > 1 & os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$cut.fdr < 0.0001;
# text.x <- log2(na.omit(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$OR[interesting.points]));
# text.y <- -log10(na.omit(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$cut.fdr[interesting.points]));
# text.labels <- na.omit(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$labels[interesting.points]);

basal.mega.enrich.scatterplot <- create.scatterplot(
    formula = -log10(as.numeric(cut.fdr)) ~ as.numeric(log2(OR)),
    data = os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01,
    col = dot.colours,
    alpha = .8,
    ylimits = c(-1.3, 24.7),
    xlimits = c(-5.2, 7.1),
    yat = seq(0, 20, 5),
    yaxis.lab = expression(10^0, 10^-5, 10^-10, 10^-15, 10^-20),
    xat = seq(-4, 6, 2),
    xaxis.lab = expression(2^-4, 2^-2, 2^0, 2^2, 2^4, 2^6),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = TRUE,
    grid.colour = 'grey85',
    cex = 1,
    pch = 21,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1,
    main.cex = 1.5,
    left.padding = 0,
    add.text = TRUE,
    # text.x = text.x,
    # text.y = text.y,
    # text.labels = text.labels,
    # text.fontface = 1,
    main = NULL,
    xlab.label = expression('Odds Ratio'),
    ylab.label = expression('FDR'),
    abline.h = -log10(0.05),
    abline.v = c(log2(1)),
    abline.col = c('grey20'),
    abline.lwd = 1,
    abline.lty = 2
    );
basal.mega.enrich.scatterplot;

save.outlier.figure(
    basal.mega.enrich.scatterplot,
    c('Figure3c', 'basal', 'scatter'),
    width = 3.3,
    height = 5.2
    );


# 2. Her2
dot.colours <- vector(length = nrow(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01));
dot.colours <- rep('grey80', nrow(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01));
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$OR) < 1 & as.numeric(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$cut.fdr) < 0.05] <- 'dodgerblue3';
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$OR) > 1 & as.numeric(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$cut.fdr) < 0.05] <- 'red2';

# interesting.points <- os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$OR > 1 & os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$cut.fdr < 0.0001;
# text.x <- log2(na.omit(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$OR[interesting.points]));
# text.y <- -log10(na.omit(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$cut.fdr[interesting.points]));
# text.labels <- na.omit(os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$labels[interesting.points]);

her2.mega.enrich.scatterplot <- create.scatterplot(
    formula = -log10(as.numeric(cut.fdr)) ~ as.numeric(log2(OR)),
    data = os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01,
    col = dot.colours,
    alpha = .8,
    ylimits = c(-1.3, 24.7),
    xlimits = c(-5.2, 7.1),
    yat = seq(0, 20, 5),
    yaxis.lab = expression(10^0, 10^-5, 10^-10, 10^-15, 10^-20),
    xat = seq(-4, 6, 2),
    xaxis.lab = expression(2^-4, 2^-2, 2^0, 2^2, 2^4, 2^6),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = TRUE,
    grid.colour = 'grey85',
    cex = 1,
    pch = 21,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1,
    main.cex = 1.5,
    left.padding = 0,
    add.text = TRUE,
    # text.x = text.x,
    # text.y = text.y,
    # text.labels = text.labels,
    # text.fontface = 1,
    main = NULL,
    xlab.label = expression('Odds Ratio'),
    ylab.label = expression('FDR'),
    abline.h = -log10(0.05),
    abline.v = c(log2(1)),
    abline.col = c('grey20'),
    abline.lwd = 1,
    abline.lty = 2
    );
her2.mega.enrich.scatterplot;

save.outlier.figure(
    her2.mega.enrich.scatterplot,
    c('Figure3c', 'her2', 'scatter'),
    width = 3.3,
    height = 5.2
    );

# 3. Luma
dot.colours <- vector(length = nrow(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01));
dot.colours <- rep('grey80', nrow(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01));
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$OR) < 1 & as.numeric(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$cut.fdr) < 0.05] <- 'dodgerblue3';
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$OR) > 1 & as.numeric(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$cut.fdr) < 0.05] <- 'red2';

# interesting.points <- os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$OR > 1 & os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$cut.fdr < 0.0001;
# text.x <- log2(na.omit(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$OR[interesting.points]));
# text.y <- -log10(na.omit(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$cut.fdr[interesting.points]));
# text.labels <- na.omit(os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$labels[interesting.points]);

luma.mega.enrich.scatterplot <- create.scatterplot(
    formula = -log10(as.numeric(cut.fdr)) ~ as.numeric(log2(OR)),
    data = os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01,
    col = dot.colours,
    alpha = .8,
    ylimits = c(-1.3, 24.7),
    xlimits = c(-5.2, 7.1),
    yat = seq(0, 20, 5),
    yaxis.lab = expression(10^0, 10^-5, 10^-10, 10^-15, 10^-20),
    xat = seq(-4, 6, 2),
    xaxis.lab = expression(2^-4, 2^-2, 2^0, 2^2, 2^4, 2^6),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = TRUE,
    grid.colour = 'grey85',
    cex = 1,
    pch = 21,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1,
    main.cex = 1.5,
    left.padding = 0,
    add.text = TRUE,
    # text.x = text.x,
    # text.y = text.y,
    # text.labels = text.labels,
    # text.fontface = 1,
    main = NULL,
    xlab.label = expression('Odds Ratio'),
    ylab.label = expression('FDR'),
    abline.h = -log10(0.05),
    abline.v = c(log2(1)),
    abline.col = c('grey20'),
    abline.lwd = 1,
    abline.lty = 2
    );
luma.mega.enrich.scatterplot;

save.outlier.figure(
    luma.mega.enrich.scatterplot,
    c('Figure3c', 'luma', 'scatter'),
    width = 3.3,
    height = 5.2
    );

# 4. Lumb
dot.colours <- vector(length = nrow(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01));
dot.colours <- rep('grey80', nrow(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01));
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$OR) < 1 & as.numeric(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$cut.fdr) < 0.05] <- 'dodgerblue3';
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$OR) > 1 & as.numeric(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$cut.fdr) < 0.05] <- 'red2';

# interesting.points <- os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$OR > 1 & os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$cut.fdr < 0.0001;
# text.x <- log2(na.omit(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$OR[interesting.points]));
# text.y <- -log10(na.omit(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$cut.fdr[interesting.points]));
# text.labels <- na.omit(os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$labels[interesting.points]);

lumb.mega.enrich.scatterplot <- create.scatterplot(
    formula = -log10(as.numeric(cut.fdr)) ~ as.numeric(log2(OR)),
    data = os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01,
    col = dot.colours,
    alpha = .8,
    ylimits = c(-1.3, 24.7),
    xlimits = c(-5.2, 7.1),
    yat = seq(0, 20, 5),
    yaxis.lab = expression(10^0, 10^-5, 10^-10, 10^-15, 10^-20),
    xat = seq(-4, 6, 2),
    xaxis.lab = expression(2^-4, 2^-2, 2^0, 2^2, 2^4, 2^6),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = TRUE,
    grid.colour = 'grey85',
    cex = 1,
    pch = 21,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1,
    main.cex = 1.5,
    left.padding = 0,
    add.text = TRUE,
    # text.x = text.x,
    # text.y = text.y,
    # text.labels = text.labels,
    # text.fontface = 1,
    main = NULL,
    xlab.label = expression('Odds Ratio'),
    ylab.label = expression('FDR'),
    abline.h = -log10(0.05),
    abline.v = c(log2(1)),
    abline.col = c('grey20'),
    abline.lwd = 1,
    abline.lty = 2
    );
lumb.mega.enrich.scatterplot;

save.outlier.figure(
    lumb.mega.enrich.scatterplot,
    c('Figure3c', 'lumb', 'scatter'),
    width = 3.3,
    height = 5.2
    );

# 5. Normal
dot.colours <- vector(length = nrow(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01));
dot.colours <- rep('grey80', nrow(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01));
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$OR) < 1 & as.numeric(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$cut.fdr) < 0.05] <- 'dodgerblue3';
dot.colours[as.numeric(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$OR) > 1 & as.numeric(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$cut.fdr) < 0.05] <- 'red2';

# interesting.points <- os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$OR > 1 & os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$cut.fdr < 0.0001;
# text.x <- log2(na.omit(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$OR[interesting.points]));
# text.y <- -log10(na.omit(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$cut.fdr[interesting.points]));
# text.labels <- na.omit(os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$labels[interesting.points]);

normal.mega.enrich.scatterplot <- create.scatterplot(
    formula = -log10(as.numeric(cut.fdr)) ~ as.numeric(log2(OR)),
    data = os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01,
    col = dot.colours,
    alpha = .8,
    ylimits = c(-1.3, 24.7),
    xlimits = c(-5.2, 7.1),
    yat = seq(0, 20, 5),
    yaxis.lab = expression(10^0, 10^-5, 10^-10, 10^-15, 10^-20),
    xat = seq(-4, 6, 2),
    xaxis.lab = expression(2^-4, 2^-2, 2^0, 2^2, 2^4, 2^6),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = TRUE,
    grid.colour = 'grey85',
    cex = 1,
    pch = 21,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1,
    main.cex = 1.5,
    left.padding = 0,
    add.text = TRUE,
    # text.x = text.x,
    # text.y = text.y,
    # text.labels = text.labels,
    # text.fontface = 1,
    main = NULL,
    xlab.label = expression('Odds Ratio'),
    ylab.label = expression('FDR'),
    abline.h = -log10(0.05),
    abline.v = c(log2(1)),
    abline.col = c('grey20'),
    abline.lwd = 1,
    abline.lty = 2
    );
normal.mega.enrich.scatterplot;

save.outlier.figure(
    normal.mega.enrich.scatterplot,
    c('Figure3c', 'normal', 'scatter'),
    width = 3.3,
    height = 5.2
    );




# Make contingency table of example gene
# Combinded contingency table
i <- 'FGFR2';

outlier.patient.all.nine.01.gene <- outlier.patient.all.nine.01[, match(rownames(all.subtype.info), colnames(outlier.patient.all.nine.01))][i, ];

# Create outlier flag vector
os.group.combine.gene <- all.subtype.info;
os.group.combine.gene$outlier_flag <- as.numeric(unlist(outlier.patient.all.nine.01.gene));
os.group.combine.gene <- na.omit(os.group.combine.gene);
os.group.combine.gene$outlier_flag <- ifelse(os.group.combine.gene$outlier_flag == 1, TRUE, FALSE);
os.group.combine.gene.table <- table(os.group.combine.gene[, c(1, 3)]);
os.group.combine.gene.table <- os.group.combine.gene.table[c(1, 2, 4, 5, 3), ];

text.size <- 1;
main.heatmap2 <- create.heatmap(
    x = os.group.combine.gene.table / sum(os.group.combine.gene.table),
    clustering = 'none',
    colour.scheme = c('white', 'dodgerblue'),
    colour.alpha = 1,
    at = seq(0, 1, 0.1),
    cell.text = data.frame(os.group.combine.gene.table)$Freq,
    text.cex = 1.2,
    text.fontface = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    grid.row = TRUE,
    grid.col = TRUE,
    ylab.label = expression("Subtype"),
    xlab.label = expression("Number of outliers per patient"),
    yaxis.lab = rownames(contingency),
    xaxis.tck = 0,
    yaxis.tck = 0,
    # xaxis.lab = colnames(plot.data.all),
    xaxis.cex = 1.1,
    yaxis.cex = 1.1,
    ylab.cex = 1.5,
    xlab.cex = 1.5,
    col.pos = rep(1:ncol(os.group.combine.gene.table), each = nrow(os.group.combine.gene.table)),
    row.pos = rep(nrow(os.group.combine.gene.table):1, times = ncol(os.group.combine.gene.table)),
    print.colour.key = FALSE,
    same.as.matrix = TRUE,
    yaxis.rot = 90,
    # xaxis.rot = 0,
    use.legacy.settings = FALSE,
    x.alternating = 2,
    xaxis.rot.top = 0,
    # xaxis.lab.top = c(colnames(outlier.num.subtype.5.patient.table.meta.5.only)[1:3], "more than 3"),
    xlab.top.cex = 1.3,
    xlab.top.y = 2.4
    );
# Row total heatmap
row.total <- create.heatmap(
    x = data.frame(rowSums(os.group.combine.gene.table), rowSums(os.group.combine.gene.table), rowSums(os.group.combine.gene.table)) / sum(os.group.combine.gene.table),
    clustering = 'none',
    colour.scheme = c('white', 'dodgerblue'),
    colour.alpha = 1,
    at = seq(0, 1, 0.1),
    cell.text = rev(rowSums(os.group.combine.gene.table)),
    text.cex = text.size,
    text.fontface = 1,
    xaxis.fontface = 1,
    xaxis.lab = c('    ', 'Total', rep('', 4)),
    xaxis.rot = 0,
    yaxis.lab = rep('', 8),
    xaxis.cex = 1,
    grid.row = TRUE,
    grid.col = FALSE,
    col.pos = rep(2, nrow(os.group.combine.gene.table)),
    row.pos = 1:nrow(os.group.combine.gene.table),
    print.colour.key = FALSE,
    same.as.matrix = TRUE,
    xaxis.tck = c(0, 0),
    yaxis.tck = c(0, 0),
    use.legacy.settings = FALSE
    );
# Column total heatmap
col.total <- create.heatmap(
    x = data.frame(colSums(os.group.combine.gene.table)) / sum(os.group.combine.gene.table),
    clustering = 'none',
    colour.scheme = c('white', 'dodgerblue'),
    colour.alpha = 1,
    at = seq(0, 1, 0.1),
    cell.text = colSums(os.group.combine.gene.table),
    text.cex = text.size,
    text.fontface = 1,
    xaxis.fontface = 1,
    xaxis.lab = rep('', 8),
    yaxis.lab = expression('Total'),
    yaxis.cex = 1.1,
    yat = 1.5,
    grid.row = FALSE,
    grid.col = TRUE,
    col.pos = 1:ncol(os.group.combine.gene.table),
    row.pos = rep(1.5, ncol(os.group.combine.gene.table)),
    print.colour.key = FALSE,
    yaxis.tck = c(0, 0),
    xaxis.tck = c(0, 0),
    xaxis.top.tck = c(0, 0),
    same.as.matrix = TRUE,
    use.legacy.settings = FALSE
    );

# Define unique subtype levels
pam50_levels <- c(1, 2, 4, 5, 3);

# Initialize results dataframe
results <- data.frame(
    subtype = pam50_levels,
    OR = NA,
    CI_low = NA,
    CI_high = NA,
    p_value = NA,
    stringsAsFactors = FALSE
    );

# Iterate for each subtype
for (k in seq_along(pam50_levels)) {
    subtype_i <- pam50_levels[k];

    # Build 2x2 contingency table
    a <- sum(os.group.combine.gene$subtype == subtype_i & os.group.combine.gene$outlier_flag == TRUE);
    b <- sum(os.group.combine.gene$subtype != subtype_i & os.group.combine.gene$outlier_flag == TRUE);
    c <- sum(os.group.combine.gene$subtype == subtype_i & os.group.combine.gene$outlier_flag == FALSE);
    d <- sum(os.group.combine.gene$subtype != subtype_i & os.group.combine.gene$outlier_flag == FALSE);

    contingency <- matrix(c(a + 1, b + 1, c + 1, d + 1), nrow = 2, byrow = TRUE);

    # Fisher's exact test
    fisher_res <- fisher.test(contingency);

    # Save results
    results$OR[k] <- ifelse(is.finite(fisher_res$estimate), fisher_res$estimate, Inf);
    results$CI_low[k] <- fisher_res$conf.int[1];
    results$CI_high[k] <- fisher_res$conf.int[2];
    results$p_value[k] <- fisher_res$p.value;
    };

# Column for the odds ratio (colored by FDR scale placeholder)
row.total.4 <- create.heatmap(
    x = data.frame(cbind(
        log2(c(results$OR)),
        log2(c(results$OR)),
        log2(c(results$OR))
        )),
    clustering = 'none',
    colour.scheme = c('#107090', 'white', '#b2402b'),
    colour.alpha = 1,
    at = seq(-2, 2, 0.1),
    cell.text = c(rev(round(c(results$OR), digits = 2))),
    text.cex = text.size,
    text.fontface = 1,
    xaxis.fontface = 1,
    xaxis.lab = c('    ', 'Odd Ratio', rep('', 4)),
    # xlab.top.label = 'FDR',
    # xlab.top.cex = 1,
    xaxis.rot = 0,
    yaxis.lab = rep('', 8),
    xaxis.cex = 1,
    grid.row = TRUE,
    grid.col = FALSE,
    col.pos = rep(2, nrow(os.group.combine.gene.table)),
    row.pos = 1:nrow(os.group.combine.gene.table),
    print.colour.key = FALSE,
    same.as.matrix = TRUE,
    xaxis.tck = c(0, 0),
    yaxis.tck = c(0, 0),
    use.legacy.settings = FALSE
    );

new.plot <- create.multipanelplot(
    list(main.heatmap2, row.total, row.total.4, col.total),
    main = expression("Cancer subtype and the number of outliers"),
    main.cex = 1.7,
    main.y = 1.6,
    resolution = 300,
    layout.height = 2,
    layout.width = 3,
    layout.skip = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    plot.objects.heights = c(10, 2),
    plot.objects.widths = c(10, 2.5, 2.5),
    x.spacing = -4,
    y.spacing = -10,
    bottom.padding = 0,
    bottom.legend.padding = -0.5,
    top.padding = 2,
    right.padding = 0
    );
new.plot;

save.outlier.figure(
    normal.mega.enrich.scatterplot,
    c('Figure3d', i, 'contingency'),
    width = 5,
    height = 5.2
    );





# Clustered heatmap
os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01 <- data.frame(
    basal = os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01$OR,
    her2 = os.group.combine.gene.subtype.enrich.her2.cut01.eachsub01$OR,
    luma = os.group.combine.gene.subtype.enrich.luma.cut01.eachsub01$OR,
    lumb = os.group.combine.gene.subtype.enrich.lumb.cut01.eachsub01$OR,
    normal = os.group.combine.gene.subtype.enrich.normal.cut01.eachsub01$OR
    );
rownames(os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01) <- rownames(os.group.combine.gene.subtype.enrich.basal.cut01.eachsub01);

distance.matrix.t <- dist(t(os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01), method = 'euclidean');
fit.t <- hclust(distance.matrix.t, method = 'ward.D2');
os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01.order <- os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01[, fit.t$order];

distance.matrix <- dist(os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01.order, method = 'euclidean');
fit <- hclust(distance.matrix, method = 'ward.D2');
os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01.order.two <- os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01.order[fit$order, ];

# Consensus clustering
# Prepare data - use log2-transformed matrix
data_matrix <- log2(t(os.group.combine.gene.subtype.enrich.allsub.cut01.eachsub01.order.two));

# Run Consensus Clustering (k from 2 to maxK)
consensus_result <- ConsensusClusterPlus(
    d = data_matrix,
    maxK = 6,
    reps = 1000,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = 'hc',
    distance = 'pearson',
    seed = 12345,
    plot = 'png',
    writeTable = TRUE,
    title = 'consensus_clustering_outlier_genes'
    );

optimal_k <- 4;

# Extract final cluster assignment for chosen k
final_clusters <- consensus_result[[optimal_k]]$consensusClass;

# Build cluster info and order samples by cluster
cluster_info <- data.frame(
    Sample = names(final_clusters),
    Cluster = final_clusters
    );
cluster_order <- order(final_clusters);
data_matrix_ordered <- data_matrix[, cluster_order];
clusters_ordered <- final_clusters[cluster_order];

# Create color bar for clusters
cluster_colors <- c('#FF6B6B', '#4ECDC4', '#FFEAA7', '#96CEB4');
names(cluster_colors) <- 1:4;
row_colors <- cluster_colors[as.character(clusters_ordered)];

cluster.heat <- create.heatmap(
    data.frame(t(clusters_ordered)),
    clustering.method = 'none',
    colour.scheme = cluster_colors,
    total.colours = 5,
    row.colour = 'black',
    col.colour = 'black',
    grid.row = TRUE,
    grid.col = TRUE,
    print.colour.key = FALSE,
    yaxis.tck = 0,
    xaxis.tck = 0
    );

# Main heatmap with ordered samples using BoutrosLab
main.heat <- BoutrosLab.plotting.general:::create.heatmap(
    x = data_matrix_ordered,
    clustering.method = 'none',
    cluster.dimensions = 'none',
    xaxis.lab = rownames(data_matrix_ordered),
    xaxis.fontface = 1,
    yaxis.lab = NULL,
    yaxis.cex = 0.2,
    xaxis.cex = 2,
    main = 'Outlier status - Consensus Clusters (k=5)',
    main.cex = 1.5,
    grid.col = TRUE,
    print.colour.key = FALSE,
    ylab.label = 'The overlap outlier gene set (clustered)',
    ylab.cex = 1.3,
    xlab.label = 'CCLE',
    force.grid.row = TRUE,
    force.grid.col = TRUE,
    grid.colour = 'white',
    xlab.cex = 1.3,
    yaxis.tck = 0,
    xaxis.tck = 0,
    axes.lwd = 0.8,
    colour.scheme = c('#107090', 'white', '#b2402b'),
    colour.centering.value = 0,
    at = seq(-3, 3, 0.001),
    colourkey.cex = 1,
    # row.colour = row_colors,
    resolution = 1000
    );

legend.sample.grob <- BoutrosLab.plotting.general:::legend.grob(
    list(
        legend = list(
            title = expression(underline('OR')),
            continuous = TRUE,
            colours = c('#107090', 'white', '#b2402b'),
            total.colours = 100,
            labels = expression(2^-3, 2^3),
            cex = 0.9,
            height = 3
            )
        ),
    title.just = 'left',
    title.fontface = 'plain'
    );

clustering.heat.all <- BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(main.heat, cluster.heat),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(2, 1),
    panel.widths = c(1, 0.1),
    ylab.padding = 2,
    xlab.to.xaxis.padding = -1.5,
    x.spacing = 0.3,
    main.cex = 0,
    xaxis.cex = 0,
    # xaxis.lab = c(rownames(data_matrix_ordered)),
    xaxis.lab = NULL,
    xaxis.fontface = 1,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.rot = 90,
    xaxis.tck = 0,
    legend = list(right = list(fun = legend.sample.grob)),
    print.new.legend = TRUE,
    resolution = 500
    );
clustering.heat.all;

save.outlier.figure(
    clustering.heat.all,
    c('Figure3c', 'XEG_cluster', 'heatmap'),
    width = 5,
    height = 5.2
    );


save.session.profile(file.path('output', 'Figure3cde.txt'));
