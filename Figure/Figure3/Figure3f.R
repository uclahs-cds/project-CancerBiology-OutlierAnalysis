### HISTORY #####################################################################
# This script performs analysis of subtype enrichment across multiple datasets
# including TCGA-BRCA, METABRIC, I-SPY2, MATADOR, and ICGC BRCA-EU. 
# Date: 2024-08-15


# Load required libraries
library(metafor);
library(BoutrosLab.plotting.general);

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
            nrow=2), 
        alternative="two.sided")$p.value;
    return(list(p.value = fisher.test.result$p.value, odd.ratio = fisher.test.result$estimate, ci = fisher.test.result$conf.int))
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
            } 
        else {
            metafor.all <- c(NA, NA, NA, NA);
            }
        
        metafor.chr.odd.ci.p <- rbind(metafor.chr.odd.ci.p, metafor.all);
        }
    
    return(metafor.chr.odd.ci.p);
    }

# Perform Fisher's test for each dataset and subtype
# TCGA-BRCA analysis
brca.results <- perform.subtype.analysis(subtype.total.outlier.num.brca, subtype.brca.status$Freq, outlier.subtype.brca.status$Freq);

# METABRIC analysis
meta.results <- perform.subtype.analysis(subtype.5.total.outlier.num.meta, subtype.5.meta.status$Freq, outlier.subtype.5.meta.status$Freq);

# I-SPY2 analysis
ispy.results <- perform.subtype.analysis(subtype.total.outlier.num.ispy, subtype.ispy.status$Freq, outlier.subtype.ispy.status$Freq);

# MATADOR analysis
matador.results <- perform.subtype.analysis(subtype.total.outlier.num.matador, subtype.matador.status$Freq, outlier.subtype.matador.status$Freq);

# ICGC analysis
icgc.results <- perform.subtype.analysis(subtype.total.outlier.num.icgc, subtype.icgc.status$Freq, outlier.subtype.icgc.status$Freq);


all.odd.subtype <- cbind(
    brca.results$odd.ratio,
    meta.results$odd.ratio,
    ispy.results$odd.ratio,
    matador.results$odd.ratio,
    icgc.results$odd.ratio
    );
all.odd.subtype.table <- as.table(all.odd.subtype);
rownames(all.odd.subtype.table) <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal');
all.odd.subtype.table <- all.odd.subtype.table[order(all.odd.subtype.table[,1], decreasing = TRUE),]

# Adjusted p-values (FDR) for each dataset
all.fdr.subtype <- cbind(
    brca.results$fdr,
    meta.results$fdr,
    ispy.results$fdr,
    matador.results$fdr,
    icgc.results$fdr
    );

# Create heatmap
odd.heat <- create.heatmap(
    x = log2(all.odd.subtype.table),
    clustering = 'none',
    colour.scheme = c("#107090", "white","#b2402b"),
    colour.alpha = 1,
    at = seq(-2, 2, 0.1),
    cell.text = round(data.frame(all.odd.subtype.table)$Freq, digits = 1),
    text.cex = 1.2,
    text.fontface = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    grid.row = TRUE,
    grid.col = TRUE,
    ylab.label = expression("Subtype"),
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
    ln(brca.results$odd.ratio),
    ln(meta.results$odd.ratio),
    ln(ispy.results$odd.ratio),
    ln(matador.results$odd.ratio),
    ln(icgc.results$odd.ratio)
    );

se.odd.list <- list(
    (ln(brca.results$ci[,2])-ln(brca.results$ci[,1])) / 3.92,
    (ln(meta.results$ci[,2])-ln(meta.results$ci[,1])) / 3.92,
    (ln(ispy.results$ci[,2])-ln(ispy.results$ci[,1])) / 3.92,
    (ln(matador.results$ci[,2])-ln(matador.results$ci[,1])) / 3.92,
    (ln(icgc.results$ci[,2])-ln(icgc.results$ci[,1])) / 3.92
    );

# Perform meta-analysis
metafor.chr.odd.ci.p <- perform.meta.analysis(ln.odd.list, se.odd.list);


metafor.chr.odd.ci.p.data <- data.frame(cbind(
    p.value = metafor.chr.odd.ci.p[,4],
    odd = metafor.chr.odd.ci.p[,1],
    ci.min = metafor.chr.odd.ci.p[,2], 
    ci.max = metafor.chr.odd.ci.p[,3],
    fdr = p.adjust(metafor.chr.odd.ci.p[,4], method = 'BH')
    ));

metafor.chr.odd.ci.p.data$labels <- as.factor(c('Basal', 'Her2', 'LuminalA', 'LuminalB', 'Normal'));

# Change the order from lowest to highest odds from meta analysis
metafor.chr.odd.ci.p.data.label.rev <- metafor.chr.odd.ci.p.data[rev(1:5),];
metafor.chr.odd.ci.p.data.label.rev <- metafor.chr.odd.ci.p.data.label.rev[order(metafor.chr.odd.ci.p.data.label.rev$odd),];


dot.colours <- vector(length=nrow(metafor.chr.odd.ci.p.data.label.rev));
dot.colours <- rep('grey70',nrow(metafor.chr.odd.ci.p.data.label.rev));
dot.colours[metafor.chr.odd.ci.p.data.label.rev$fdr < 0.05 & metafor.chr.odd.ci.p.data.label.rev$odd < 1] <- "#107090";
dot.colours[metafor.chr.odd.ci.p.data.label.rev$fdr < 0.05 & metafor.chr.odd.ci.p.data.label.rev$odd > 1] <- "red3";

# segment plot
metafor.all.segplot <- BoutrosLab.plotting.general::create.segplot(
    formula = labels ~ log2(ci.min) + log2(ci.max),
    data = metafor.chr.odd.ci.p.data.label.rev,
    centers = log2(metafor.chr.odd.ci.p.data.label.rev$odd),
    segments.col = dot.colours,
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -3,
    xright.rectangle = 13,
    col.rectangle = "grey",
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


# Save the multi plot as a PDF
pdf(
    file = generate.filename(
        'subtype', 
        'multi', 
        'pdf'
        ), 
    width = 7, 
    height = 5
    );
multi.gene;
dev.off();

# Save the multi plot as a PNG
png(
    file = generate.filename(
        'subtype', 
        'multi', 
        'png'
        ), 
    width = 7, 
    height = 5,
    unit = 'in', 
    res = 1200
    );
multi.gene;
dev.off();



