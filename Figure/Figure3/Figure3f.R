### HISTORY ######################################################################
# This script performs analysis of subtype enrichment across multiple datasets
# including TCGA-BRCA, METABRIC, I-SPY2, MATADOR, and ICGC BRCA-EU.
# Date: 2024-08-15

### DESCRIPTION ##################################################################
# This script analyzes breast cancer subtype enrichment across multiple datasets.
# It processes and standardizes subtype data, conducts Fisher's exact tests,
# calculates odds ratios, performs meta-analysis, and generates visualizations
# including heatmaps and forest plots to compare subtype enrichment across datasets.

### PREAMBLE #####################################################################
# Load required libraries
library(metafor);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-09-11_Figure3e-i.rda'));


# 1. TCGA-BRCA
brca.clinic.order.data <- data.frame(brca.clinic.order$Subtype);
brca.clinic.order.data[is.na(brca.clinic.order.data$brca.clinic.order.Subtype), ] <- 6;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_Basal', ] <- 1;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_Her2', ] <- 2;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_LumA', ] <- 3;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_LumB', ] <- 4;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_Normal', ] <- 5;
brca.clinic.order.data.num <- data.frame(as.numeric(brca.clinic.order.data$brca.clinic.order.Subtype));
rownames(brca.clinic.order.data.num) <- colnames(outlier.patient.tag.01.brca);

outlier.patient.tag.01.brca.sum <- apply(outlier.patient.tag.01.brca, 2, sum);
subtype.total.outlier.num.brca <- data.frame(cbind(
    subtype = brca.clinic.order.data.num,
    outlier = outlier.patient.tag.01.brca.sum
    ));
colnames(subtype.total.outlier.num.brca) <- c('subtype', 'outlier');
subtype.total.outlier.num.1.brca <- subtype.total.outlier.num.brca;
subtype.total.outlier.num.1.brca$outlier[subtype.total.outlier.num.1.brca$outlier > 0] <- 1;
outlier.subtype.brca.status <- data.frame(table(subtype.total.outlier.num.1.brca));
subtype.brca.status <- data.frame(table(brca.clinic.order$Subtype));


# 2. METABRIC
meta.clinic.5.order.data <- data.frame(meta.clinic.5.order.combine$pam50);
meta.clinic.5.order.data[is.na(meta.clinic.5.order.data$meta.clinic.5.order.combine.pam50), ] <- 6;
meta.clinic.5.order.data[meta.clinic.5.order.data$meta.clinic.5.order.combine.pam50 == 'Basal', ] <- 1;
meta.clinic.5.order.data[meta.clinic.5.order.data$meta.clinic.5.order.combine.pam50 == 'Her2', ] <- 2;
meta.clinic.5.order.data[meta.clinic.5.order.data$meta.clinic.5.order.combine.pam50 == 'LumA', ] <- 3;
meta.clinic.5.order.data[meta.clinic.5.order.data$meta.clinic.5.order.combine.pam50 == 'LumB', ] <- 4;
meta.clinic.5.order.data[meta.clinic.5.order.data$meta.clinic.5.order.combine.pam50 == 'Normal', ] <- 5;
meta.clinic.5.order.data.num <- data.frame(as.numeric(meta.clinic.5.order.data$meta.clinic.5.order.combine.pam50));
rownames(meta.clinic.5.order.data.num) <- colnames(outlier.patient.tag.01.meta);

outlier.patient.tag.01.meta.sum <- apply(outlier.patient.tag.01.meta, 2, sum);
subtype.5.total.outlier.num.meta <- data.frame(cbind(
    subtype.5 = meta.clinic.5.order.data.num,
    outlier = outlier.patient.tag.01.meta.sum
    ));
colnames(subtype.5.total.outlier.num.meta) <- c('subtype.5', 'outlier');
subtype.5.total.outlier.num.1.meta <- subtype.5.total.outlier.num.meta;
subtype.5.total.outlier.num.1.meta$outlier[subtype.5.total.outlier.num.1.meta$outlier > 0] <- 1;
outlier.subtype.5.meta.status <- data.frame(table(subtype.5.total.outlier.num.1.meta));
subtype.5.meta.status <- data.frame(table(meta.clinic.5.order$subtype));
subtype.5.meta.status <- subtype.5.meta.status[match(c('Basal', 'Her2', 'LumA', 'LumB', 'Normal'), subtype.5.meta.status$Var1), ];


# 3. ICGC BRCA-EU
icgc.clinic.subtype.order <- icgc.clinic.order[match(colnames(outlier.patient.tag.01.icgc), icgc.clinic.order$sample), ];

icgc.clinic.subtype.order.data <- data.frame(as.character(icgc.clinic.subtype.order$subtype));
icgc.clinic.subtype.order.data[is.na(icgc.clinic.subtype.order.data$as.character.icgc.clinic.subtype.order.subtype.), ] <- 6;
icgc.clinic.subtype.order.data[icgc.clinic.subtype.order.data$as.character.icgc.clinic.subtype.order.subtype. == 'Basal', ] <- 1;
icgc.clinic.subtype.order.data[icgc.clinic.subtype.order.data$as.character.icgc.clinic.subtype.order.subtype. == 'Her2', ] <- 2;
icgc.clinic.subtype.order.data[icgc.clinic.subtype.order.data$as.character.icgc.clinic.subtype.order.subtype. == 'LumA', ] <- 3;
icgc.clinic.subtype.order.data[icgc.clinic.subtype.order.data$as.character.icgc.clinic.subtype.order.subtype. == 'LumB', ] <- 4;
icgc.clinic.subtype.order.data[icgc.clinic.subtype.order.data$as.character.icgc.clinic.subtype.order.subtype. == 'Normal', ] <- 5;

subtype.total.outlier.num.1.icgc <- subtype.total.outlier.num.icgc;
subtype.total.outlier.num.1.icgc$outlier[subtype.total.outlier.num.1.icgc$outlier > 0] <- 1;
outlier.subtype.icgc.status <- data.frame(table(subtype.total.outlier.num.1.icgc));
total.subtype <- icgc.clinic.subtype.order.data[, 1];
total.subtype <- total.subtype[total.subtype != '']
subtype.icgc.status <- data.frame(table(total.subtype));


# 4. I-SPY2
# ispy.clinic <-  read.csv(file = '/Users/jee/Documents/1.Project/ISPY/clinical/1-s2.0-S1535610822002161-mmc3.csv', header = TRUE, stringsAsFactors = F, sep = ',');
# rownames(ispy.clinic) <- paste('X', ispy.clinic$Patient.Identifier, sep = '');
ispy.clinic.order <- ispy.clinic[colnames(outlier.patient.tag.01.ispy), ];

ispy.clinic.order.data <- data.frame(ispy.clinic.order$PAM50.Subtype);
ispy.clinic.order.data[is.na(ispy.clinic.order.data$ispy.clinic.order.PAM50.Subtype), ] <- 6;
ispy.clinic.order.data[ispy.clinic.order.data$ispy.clinic.order.PAM50.Subtype == 'Basal', ] <- 1;
ispy.clinic.order.data[ispy.clinic.order.data$ispy.clinic.order.PAM50.Subtype == 'Her2', ] <- 2;
ispy.clinic.order.data[ispy.clinic.order.data$ispy.clinic.order.PAM50.Subtype == 'LumA', ] <- 3;
ispy.clinic.order.data[ispy.clinic.order.data$ispy.clinic.order.PAM50.Subtype == 'LumB', ] <- 4;
ispy.clinic.order.data[ispy.clinic.order.data$ispy.clinic.order.PAM50.Subtype == 'Normal', ] <- 5;
ispy.clinic.order.data.num <- data.frame(as.numeric(ispy.clinic.order.data$ispy.clinic.order.PAM50.Subtype));
rownames(ispy.clinic.order.data.num) <- colnames(outlier.patient.tag.01.ispy);

outlier.patient.tag.01.sum.ispy <- apply(outlier.patient.tag.01.ispy, 2, sum);
subtype.total.outlier.num.ispy <- data.frame(cbind(
    subtype = ispy.clinic.order.data.num,
    outlier = outlier.patient.tag.01.sum.ispy
    ));
colnames(subtype.total.outlier.num.ispy) <- c('subtype', 'outlier');
subtype.total.outlier.num.1.ispy <- subtype.total.outlier.num.ispy;
subtype.total.outlier.num.1.ispy$outlier[subtype.total.outlier.num.1.ispy$outlier > 0] <- 1;
outlier.subtype.ispy.status <- data.frame(table(subtype.total.outlier.num.1.ispy));
subtype.ispy.status <- data.frame(table(ispy.clinic.order$PAM50.Subtype));


# 5. MATADOR

metador.clinic.order.data <- data.frame(metador.clinic.order$subtype);
metador.clinic.order.data[is.na(metador.clinic.order.data$metador.clinic.order.PAM50.Subtype),] <- 6;
metador.clinic.order.data[metador.clinic.order.data$metador.clinic.order.PAM50.Subtype == 'Basal',] <- 1;
metador.clinic.order.data[metador.clinic.order.data$metador.clinic.order.PAM50.Subtype == 'Her2',] <- 2;
metador.clinic.order.data[metador.clinic.order.data$metador.clinic.order.PAM50.Subtype == 'LumA',] <- 3;
metador.clinic.order.data[metador.clinic.order.data$metador.clinic.order.PAM50.Subtype == 'LumB',] <- 4;
metador.clinic.order.data[metador.clinic.order.data$metador.clinic.order.PAM50.Subtype == 'Normal',] <- 5;
metador.clinic.order.data.num <- data.frame(as.numeric(metador.clinic.order.data$metador.clinic.order.subtype));
rownames(metador.clinic.order.data.num) <- colnames(outlier.patient.tag.01.matador);

outlier.patient.tag.01.matador.sum <- apply(outlier.patient.tag.01.matador, 2, sum);
subtype.total.outlier.num.matador <- data.frame(cbind(subtype = metador.clinic.order.data.num,
                                              outlier = outlier.patient.tag.01.matador.sum));
colnames(subtype.total.outlier.num.matador) <- c("subtype", "outlier");

subtype.total.outlier.num.1.matador <- subtype.total.outlier.num.matador;
subtype.total.outlier.num.1.matador$outlier[subtype.total.outlier.num.1.matador$outlier > 0] <- 1;
outlier.subtype.matador.status <- data.frame(table(subtype.total.outlier.num.1.matador));
subtype.matador.status <- data.frame(table(subtype.total.outlier.num.matador$subtype));



# Function to perform Fisher's exact test and odds ratio calculation
perform.fisher.test <- function(total.patients, total.outliers, subtype.freq, outlier.freq) {
    fisher.test.result <- fisher.test(
        matrix(
            c(
                outlier.freq,
                total.outliers - outlier.freq,
                subtype.freq - outlier.freq,
                total.patients - total.outliers - subtype.freq + outlier.freq
                ),
            nrow = 2
            ),
        alternative = 'two.sided'
        )
    list(p.value = fisher.test.result$p.value, odd.ratio = fisher.test.result$estimate, ci = fisher.test.result$conf.int)
    }

# Function to perform subtype analysis for each dataset
perform.subtype.analysis <- function(subtype.data, subtype.freq, outlier.freq) {
    total.patients <- nrow(subtype.data); # Total number of patients
    total.outliers <- sum(subtype.data$outlier > 0); # Total number of patients with outliers

    p.value.list <- NULL;
    odd.ratio.list <- NULL;
    ci.list <- NULL;

    for (i in 1:5) {
        results <- perform.fisher.test(total.patients, total.outliers, subtype.freq[i], outlier.freq[i]);
        p.value.list <- c(p.value.list, results$p.value);
        odd.ratio.list <- c(odd.ratio.list, results$odd.ratio);
        ci.list <- rbind(ci.list, results$ci);
        }

    p.value.fdr <- p.adjust(p.value.list, method = 'BH');

    return(list(p.value = p.value.list, odd.ratio = odd.ratio.list, ci = ci.list, fdr = p.value.fdr));
    }

# Function to perform meta-analysis
perform.meta.analysis <- function(ln.odd.list, se.odd.list) {
    chr.odd.se <- list();

    for (i in 1:length(ln.odd.list[[1]])) {
        chr.odd <- sapply(ln.odd.list, function(x) x[i]);
        chr.se <- sapply(se.odd.list, function(x) x[i]);
        chr.all <- data.frame(cbind(chr.odd, chr.se));
        chr.odd.se[[i]] <- chr.all;
        }

    metafor.chr.odd.ci.p <- NULL;

    for (i in 1:length(ln.odd.list[[1]])) {
        chr.odd.se.sample <- chr.odd.se[[i]];
        chr.odd.se.sample.inf <- chr.odd.se.sample[!is.infinite(chr.odd.se.sample$chr.odd) & !is.infinite(chr.odd.se.sample$chr.se), ];

        if (nrow(chr.odd.se.sample.inf) > 1) {
            metafor.chr <- rma.uni(yi = chr.odd, sei = chr.se, data = chr.odd.se.sample.inf, method = 'DL');
            metafor.all <- c(exp(metafor.chr$beta), exp(metafor.chr$ci.lb), exp(metafor.chr$ci.ub), metafor.chr$pval);
            } else {
            metafor.all <- c(NA, NA, NA, NA);
            }

        metafor.chr.odd.ci.p <- rbind(metafor.chr.odd.ci.p, metafor.all);
        }

    return(metafor.chr.odd.ci.p);
    }

# Perform Fisher's test for each dataset and subtype
# TCGA-BRCA analysis
brca.results <- perform.subtype.analysis(subtype.total.outlier.num.brca, subtype.brca.status$Freq, outlier.subtype.brca.status[outlier.subtype.brca.status$outlier == 1, ]$Freq);

# METABRIC analysis
meta.results <- perform.subtype.analysis(subtype.5.total.outlier.num.meta, subtype.5.meta.status$Freq, outlier.subtype.5.meta.status[outlier.subtype.5.meta.status$outlier == 1, ]$Freq);

# I-SPY2 analysis
ispy.results <- perform.subtype.analysis(subtype.total.outlier.num.ispy, subtype.ispy.status$Freq, outlier.subtype.ispy.status[outlier.subtype.ispy.status$outlier == 1, ]$Freq);

# MATADOR analysis
matador.results <- perform.subtype.analysis(subtype.total.outlier.num.matador, subtype.matador.status$Freq, outlier.subtype.matador.status[outlier.subtype.matador.status$outlier == 1, ]$Freq);

# ICGC analysis
icgc.results <- perform.subtype.analysis(subtype.total.outlier.num.icgc, subtype.icgc.status$Freq, outlier.subtype.icgc.status[outlier.subtype.icgc.status$outlier == 1, ]$Freq);


all.odd.subtype <- cbind(
    brca.results$odd.ratio,
    meta.results$odd.ratio,
    ispy.results$odd.ratio,
    matador.results$odd.ratio,
    icgc.results$odd.ratio
    );
all.odd.subtype.table <- as.table(all.odd.subtype);
rownames(all.odd.subtype.table) <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal');
all.odd.subtype.table <- all.odd.subtype.table[order(all.odd.subtype.table[, 1], decreasing = TRUE), ]

# Create heatmap
odd.heat <- create.heatmap(
    x = log2(all.odd.subtype.table),
    clustering = 'none',
    colour.scheme = c('#107090', 'white', '#b2402b'),
    colour.alpha = 1,
    at = seq(-2, 2, 0.1),
    cell.text = round(data.frame(all.odd.subtype.table)$Freq, digits = 1),
    text.cex = 1.2,
    text.fontface = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    grid.row = TRUE,
    grid.col = TRUE,
    ylab.label = expression('Subtype'),
    yaxis.lab = rownames(all.odd.subtype.table),
    xaxis.lab.top = c('TCGA-BRCA', 'METABRIC', 'I-SPY2', 'MATADOR', 'ICGC BRCA-EU'),
    xaxis.tck = 0,
    yaxis.tck = 0,
    xaxis.cex = 1.1,
    yaxis.cex = 1.1,
    ylab.cex = 1.5,
    xlab.cex = 1.5,
    col.pos = rep(1:ncol(all.odd.subtype.table), each = nrow(all.odd.subtype.table)),
    row.pos = rep(nrow(all.odd.subtype.table):1, times = ncol(all.odd.subtype.table)),
    print.colour.key = FALSE,
    same.as.matrix = TRUE,
    yaxis.rot = 90,
    use.legacy.settings = FALSE,
    x.alternating = 2,
    xaxis.rot.top = 0,
    xlab.top.cex = 1.3,
    xlab.top.y = 2.4
    );


ln.odd.list <- list(
    log(brca.results$odd.ratio),
    log(meta.results$odd.ratio),
    log(ispy.results$odd.ratio),
    log(matador.results$odd.ratio),
    log(icgc.results$odd.ratio)
    );

se.odd.list <- list(
    (log(brca.results$ci[, 2]) - log(brca.results$ci[, 1])) / 3.92,
    (log(meta.results$ci[, 2]) - log(meta.results$ci[, 1])) / 3.92,
    (log(ispy.results$ci[, 2]) - log(ispy.results$ci[, 1])) / 3.92,
    (log(matador.results$ci[, 2]) - log(matador.results$ci[, 1])) / 3.92,
    (log(icgc.results$ci[, 2]) - log(icgc.results$ci[, 1])) / 3.92
    );

# Perform meta-analysis
metafor.chr.odd.ci.p <- perform.meta.analysis(ln.odd.list, se.odd.list);


metafor.chr.odd.ci.p.data <- data.frame(cbind(
    p.value = metafor.chr.odd.ci.p[, 4],
    odd = metafor.chr.odd.ci.p[, 1],
    ci.min = metafor.chr.odd.ci.p[, 2],
    ci.max = metafor.chr.odd.ci.p[, 3],
    fdr = p.adjust(metafor.chr.odd.ci.p[, 4], method = 'BH')
    ));

metafor.chr.odd.ci.p.data$labels <- as.factor(c('Basal', 'Her2', 'LuminalA', 'LuminalB', 'Normal'));

# Change the order from lowest to highest odds from meta analysis
metafor.chr.odd.ci.p.data.label.rev <- metafor.chr.odd.ci.p.data[rev(1:5), ];
metafor.chr.odd.ci.p.data.label.rev <- metafor.chr.odd.ci.p.data.label.rev[order(metafor.chr.odd.ci.p.data.label.rev$odd), ];


dot.colours <- vector(length = nrow(metafor.chr.odd.ci.p.data.label.rev));
dot.colours <- rep('grey70', nrow(metafor.chr.odd.ci.p.data.label.rev));
dot.colours[metafor.chr.odd.ci.p.data.label.rev$fdr < 0.05 & metafor.chr.odd.ci.p.data.label.rev$odd < 1] <- '#107090';
dot.colours[metafor.chr.odd.ci.p.data.label.rev$fdr < 0.05 & metafor.chr.odd.ci.p.data.label.rev$odd > 1] <- 'red3';

# segment plot
metafor.all.segplot <- BoutrosLab.plotting.general::create.segplot(
    formula = labels ~ log2(ci.min) + log2(ci.max),
    data = metafor.chr.odd.ci.p.data.label.rev,
    centers = log2(metafor.chr.odd.ci.p.data.label.rev$odd),
    segments.col = dot.colours,
    main.cex = 0,
    yaxis.fontface = 1,
    xlab.cex = 1.3,
    xlab.label = expression('Odds Ratio'),
    ylab.cex = 0,
    yaxis.cex = 0,
    xaxis.cex = 1,
    xaxis.lab = c(0, 0.25, 0.5, 1, 2, 4),
    xaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -3,
    xright.rectangle = 13,
    ybottom.rectangle = seq(1.5, 23.5, 2),
    ytop.rectangle = seq(2.5, 24.5, 2),
    # set rectangle colour
    col.rectangle = 'grey',
    # set rectangle alpha (transparency)
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
    );
metafor.all.segplot;

# Combine heatmap and segment plots
multi.gene <- create.multipanelplot(
    list(odd.heat, metafor.all.segplot),
    main.cex = 0,
    resolution = 300,
    layout.height = 1,
    layout.width = 2,
    layout.skip = c(FALSE, FALSE),
    plot.objects.widths = c(3, 2),
    ylab.axis.padding = -10,
    x.spacing = 6,
    right.legend.padding = 0
    );

multi.gene;

save.outlier.figure(
    multi.gene,
    c('Figure3f', 'subtype', 'multi'),
    width = 7,
    height = 5
    );

save.session.profile(file.path('output', 'Figure3f.txt'));
