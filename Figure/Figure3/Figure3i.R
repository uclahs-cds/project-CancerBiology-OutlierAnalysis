### HISTORY #####################################################################
# This script performs a meta-meta analysis to combine hazard ratios 
# across different breast cancer subtypes (e.g., Basal, Her2, LuminalA, LuminalB, Normal) 
# Date: 2024-08-15



# Load required libraries
library(metafor);


combine.surv.data.subtype <- combine.surv.data[2:6, ];
ln.hr.combine.surv.data.subtype <- ln(combine.surv.data.subtype$mean);
se.hr.combine.surv.data.subtype <- (ln(combine.surv.data.subtype$upper) - ln(combine.surv.data.subtype$lower)) / 3.92;

# Perform meta-meta analysis
res.meta.meta.all.combine.cox <- rma.uni(yi = ln.hr.combine.surv.data.subtype, sei = se.hr.combine.surv.data.subtype, method = "DL");
res.meta.meta.all.combine <- summary(res.meta.meta.all.combine.cox);


all.data.combine.meta <- data.frame(
    mean = exp(as.numeric(res.meta.meta.all.combine$beta)),
    lower = exp(as.numeric(res.meta.meta.all.combine$ci.lb)),
    upper = exp(as.numeric(res.meta.meta.all.combine$ci.ub)),
    Features = 'Outlier patients',
    number = paste(as.character(sum(os.group.combine[,'patient'] == 2)), '/', as.character(length(os.group.combine[,'patient'])), sep = ''),
    HR = as.character(round(exp(as.numeric(res.meta.meta.all.combine$beta)), digits = 2)),
    Pvalue = as.character(round(as.numeric(res.meta.meta.all.combine$pval), digits = 4)),
    event = paste(as.character(sum(os.group.combine[os.group.combine$patient == 2, 'status'] == 1)), '/', as.character(sum(os.group.combine[, 'status'] == 1)), sep = ''),
    assumption = as.character(round(cox.rf.group.assumption.sub$table[1, 3], digits = 3))
    );


combine.surv.data.meta <- rbind(
    all = all.data.combine.meta,
    basal = Basal.data.combine,
    her2 = Her2.data.combine,
    luma = LumA.data.combine,
    lumb = LumB.data.combine,
    normal = Normal.data.combine
    );

combine.surv.data.meta$Features <- c('All patients', 'Basal', 'Her2', 'LuminalA', 'LuminalB', 'Normal');


combine.surv.data.meta.rev <- combine.surv.data.meta[c(6, 5, 4, 3, 2, 1), ];

dot.colours <- rep('grey70', nrow(combine.surv.data.meta.rev));
dot.colours[as.numeric(combine.surv.data.meta.rev$Pvalue) < 0.05 & combine.surv.data.meta.rev$mean < 1] <- 'dodgerblue2';
dot.colours[as.numeric(combine.surv.data.meta.rev$Pvalue) < 0.05 & combine.surv.data.meta.rev$mean > 1] <- 'red';

# Generate segment plot
merge.surv.seg.log <- BoutrosLab.plotting.general::create.segplot(
    formula = as.factor(Features) ~ log2(lower) + log2(upper),
    data = combine.surv.data.meta.rev,
    centers = log2(combine.surv.data.meta.rev$mean),
    main = NULL,
    ylab.label = expression('Subtype'),
    yaxis.lab = combine.surv.data.rev$Features,
    yaxis.fontface = 1,
    xlab.label = expression('Hazard Ratio'),
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    yaxis.cex = 1,
    xaxis.cex = 1,
    xlimits = c(-0.8, 1.8),
    xaxis.lab = c('0.75', '1', '1.5', '2', '2.5'),
    xaxis.fontface = 1,
    xat = c(log2(0.75), log2(1), log2(1.5), log2(2), log2(2.5)),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    segments.col = dot.colours,
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -1,
    xright.rectangle = 4,
    ybottom.rectangle = seq(1.5, 23.5, 2),
    ytop.rectangle = seq(2.5, 24.5, 2),
    col.rectangle = "grey",
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
    );


# Save the segment plot as a PDF
pdf(
    file = generate.filename(
        'survival_subtype', 
        'segment', 
        'pdf'
        ), 
    width = 5, 
    height = 3.8
    );
merge.surv.seg.log;
dev.off();

# Save the segment plot as a PNG
png(
    file = generate.filename(
        'survival_subtype', 
        'segment', 
        'png'
        ), 
    width = 5, 
    height = 3.8,
    unit = 'in', 
    res = 1200
    );
merge.surv.seg.log;
dev.off();

