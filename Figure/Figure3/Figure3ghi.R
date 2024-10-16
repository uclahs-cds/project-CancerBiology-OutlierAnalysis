### HISTORY #####################################################################
# This script processes survival data from the TCGA-BRCA and METABRIC datasets
# to analyze overall survival by combining data from both sources. The analysis
# includes Kaplan-Meier survival plots among different patient groups.
# Date: 2024-08-15

### DESCRIPTION #################################################################
# The script loads and processes survival data from TCGA-BRCA and METABRIC datasets,
# combines them, and performs a Kaplan-Meier survival analysis. It creates a plot
# comparing survival outcomes between outlier and non-outlier patient groups.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.survival);
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

### 1. TCGA-BRCA
os.data.brca <- data.frame(cbind(
    status = substr(brca.clinic.order$Overall.Survival.Status, 1, 1),
    os = brca.clinic.order$Overall.Survival..Months.,
    sum = apply(outlier.patient.tag.01.brca, 2, sum),
    pam50 = sub('^BRCA_', '', brca.clinic.order$Subtype),
    age = brca.clinic.order$Diagnosis.Age
    ));

os.data.meta <- data.frame(cbind(
    status = substr(meta.clinic.5.order.combine$Overall.Survival.Status, 1, 1),
    os = meta.clinic.5.order.combine$Overall.Survival..Months.,
    sum = apply(outlier.patient.tag.01.meta, 2, sum),
    pam50 = meta.clinic.5.order.combine$pam50,
    age = meta.clinic.5.order.combine$Age.at.Diagnosis
    ));

os.group.combine <- rbind(
    os.data.brca,
    os.data.meta
    );

# Ensure 'status', 'os', 'sum', and 'age' are numeric
os.group.combine$status <- as.numeric(os.group.combine$status)
os.group.combine$os <- as.numeric(os.group.combine$os)
os.group.combine$sum <- as.numeric(os.group.combine$sum)
os.group.combine$age <- as.numeric(os.group.combine$age)

# Filter out any NC pam50 values, then convert that to a factor
os.group.combine <- os.group.combine[!(os.group.combine$pam50 %in% 'NC'), ];
os.group.combine$pam50 <- relevel(as.factor(os.group.combine$pam50), ref = 'LumA');

# Group risk groups based on total sums
os.group.combine$patient <- ifelse(os.group.combine$sum == 0, 1, 2);

os.group.combine <- na.omit(os.group.combine);

message('First plot')
message(paste(str(os.group.combine)))

### 4. Kaplan-Meier Survival Analysis
km.os.group.combine <- create.km.plot(
    survival.object = Surv(os.group.combine$os, os.group.combine$status),
    main = as.expression(substitute(paste('Kaplan-Meier estimate (Combined datasets)'))),
    show.risktable = TRUE,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    xlab.label = expression('Overall survival (Months)'),
    ylab.label = expression('Estimated proportion'),
    xlimits = c(0, 280),
    xat = seq(0, 240, 80),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.5,
    key.stats.cex = 1.1,
    patient.groups = as.factor(os.group.combine$patient),
    risktable.fontsize = 11.5,
    show.key.groups = TRUE,
    risk.label.pos = -70,
    ylab.axis.padding = 2,
    risk.label.fontface = 1,
    left.padding = 5.5,
    key.groups.labels = rev(c('Outlier patients', 'Non-outlier patients')),
    key.groups.cex = 1,
    line.colours = rev(c('red3', 'dodgerblue3'))
    );
km.os.group.combine;

save.outlier.figure(
    km.os.group.combine,
    c('Figure3ghi', 'os', 'merge', 'km'),
    width = 7.5,
    height = 7
    );

i <- 'Basal'

os.group.basal <- os.group.combine[os.group.combine$pam50 %in% 'Basal', ];

km.os.group.combine <- create.km.plot(
    survival.object = Surv(os.group.basal$os, os.group.basal$status),
    main = as.expression(substitute(paste('Kaplan-Meier estimate (Combined datasets) - ', var), list(var = 'Basal'))),
    show.risktable = TRUE,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    xlab.label = expression('Overall survival (Months)'),
    ylab.label = expression('Estimated proportion'),
    xlimits = c(0, 280),
    xat = seq(0, 240, 80),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.5,
    key.stats.cex = 1.1,
    patient.groups = as.factor(os.group.basal$patient),
    risktable.fontsize = 11.5,
    show.key.groups = TRUE,
    risk.label.pos = -70,
    ylab.axis.padding = 2,
    risk.label.fontface = 1,
    left.padding = 5.5,
    key.groups.labels = c('Non-outlier patients', 'Outlier patients'),
    key.groups.cex = 1,
    line.colours = rev(c('red3', 'dodgerblue3'))
    );
km.os.group.combine;

save.outlier.figure(
    km.os.group.combine,
    c('Figure3ghi', i, 'km'),
    width = 7.5,
    height = 7
    );

# This section combines TCGA-BRCA and METABRIC datasets, performs Cox proportional
# hazards regression for each breast cancer subtype, and then conducts a meta-analysis
# to combine the hazard ratios. It generates a segment plot to visualize the results.

os.model.subtype <- coxph(Surv(os, status == TRUE) ~ patient, data = os.group.combine);

cox.rf.group.assumption.sub <- cox.zph(os.model.subtype);

summary.cox.sub <- summary(os.model.subtype);

all.data.combine <- data.frame(
    mean = as.numeric(summary.cox.sub$conf.int[1]),
    lower = as.numeric(summary.cox.sub$conf.int[3]),
    upper = as.numeric(summary.cox.sub$conf.int[4]),
    Features = 'Outlier patients',
    number = paste(as.character(sum(os.group.combine[, 'patient'] == 2)), '/', as.character(length(os.group.combine[, 'patient'])), sep = ''),
    HR = as.character(round(as.numeric(summary.cox.sub$conf.int[1]), digits = 2)),
    Pvalue = as.character(round(as.numeric(summary.cox.sub$coefficients[5]), digits = 4)),
    event = paste(as.character(sum(os.group.combine[os.group.combine$patient == 2, 'status'] == 1)), '/', as.character(sum(os.group.combine[, 'status'] == 1)), sep = ''),
    assumption = as.character(round(cox.rf.group.assumption.sub$table[1, 3], digits = 3))
    );

input.subtype <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal')

for (i in input.subtype) {
    os.model.subtype <- coxph(Surv(os, status == TRUE) ~ patient, data = os.group.combine[os.group.combine$pam50 == i, ]);
    cox.rf.group.assumption.sub <- cox.zph(os.model.subtype);
    summary.cox.sub <- summary(os.model.subtype);

    data <- data.frame(
        mean = as.numeric(summary.cox.sub$conf.int[1]),
        lower = as.numeric(summary.cox.sub$conf.int[3]),
        upper = as.numeric(summary.cox.sub$conf.int[4]),
        Features = 'Outlier patients',
        number = paste(as.character(sum(os.group.combine[os.group.combine$pam50 == i, 'patient'] == 2)), '/', as.character(length(os.group.combine[os.group.combine$pam50 == i, 'patient'])), sep = ''),
        HR = as.character(round(as.numeric(summary.cox.sub$conf.int[1]), digits = 2)),
        Pvalue = as.character(round(as.numeric(summary.cox.sub$coefficients[5]), digits = 3)),
        event = paste(as.character(sum(os.group.combine[os.group.combine$pam50 == i & os.group.combine$patient == 2, 'status'] == 1)), '/', as.character(sum(os.group.combine[os.group.combine$pam50 == i, 'status'] == 1)), sep = ''),
        assumption = as.character(round(cox.rf.group.assumption.sub$table[1, 3], digits = 3))
        );

    data.name <- paste(i, '.data.combine', sep = '');
    assign(data.name, data);
    }


combine.surv.data <- rbind(
    all = all.data.combine,
    basal = Basal.data.combine,
    her2 = Her2.data.combine,
    luma = LumA.data.combine,
    lumb = LumB.data.combine,
    normal = Normal.data.combine
    );
combine.surv.data$Features <- c('All patients', 'Basal', 'Her2', 'LuminalA', 'LuminalB', 'Normal');

combine.surv.data.subtype <- combine.surv.data[2:6, ];
ln.hr.combine.surv.data.subtype <- log(combine.surv.data.subtype$mean);
se.hr.combine.surv.data.subtype <- (log(combine.surv.data.subtype$upper) - log(combine.surv.data.subtype$lower)) / 3.92;

combine.surv.data.meta <- rbind(
    all = all.data.combine,
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
    yaxis.lab = combine.surv.data.meta.rev$Features,
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
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
    );


save.outlier.figure(
    merge.surv.seg.log,
    c('Figure3ghi', 'survival', 'subtype', 'segment'),
    width = 5,
    height = 3.8
    );


save.session.profile(file.path('output', 'Figure3ghi.txt'));
