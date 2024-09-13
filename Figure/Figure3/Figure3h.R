### HISTORY #####################################################################
# This script generates Kaplan-Meier survival curves for each PAM50 subtype 
# in combined datasets (TCGA-BRCA and METABRIC). 
# The analysis corresponds to Figure 3g in the associated study.
# Date: 2024-08-15

library(BoutrosLab.plotting.survival);

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));


### 3. Combine TCGA-BRCA and METABRIC Datasets
os.group.combine <- data.frame(rbind(
    os.group.brca,
    os.group.meta
    ));


os.group.combine$pam50 <- factor(os.group.combine$pam50);
os.group.combine$pam50 <- relevel(os.group.combine$pam50, ref = "LumA");
os.group.combine <- os.group.combine[!(os.group.combine$pam50 %in% "NC"),];
os.group.combine <- na.omit(os.group.combine);
os.group.combine$pam50 <- factor(os.group.combine$pam50);



i <- 'Basal'


km.os.group.combine <- create.km.plot(
    survival.object = Surv(os.group.combine[os.group.combine$pam50 == i,]$os, as.numeric(os.group.combine[os.group.combine$pam50 == i,]$status)),
    main = as.expression(substitute(paste('Kaplan-Meier estimate (Combined datasets) - ',var), list(var = i))),
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
    patient.groups = as.factor(os.group.combine[os.group.combine$pam50 == i,]$patient),
    risktable.fontsize = 11.5,
    show.key.groups = TRUE,
    risk.label.pos = -70,
    ylab.axis.padding = 2,
    risk.label.fontface = 1,
    left.padding = 5.5,
    key.groups.labels = c("Non-outlier patients","Outlier patients"),
    key.groups.cex = 1,
    line.colours = c('red3', 'dodgerblue3')
    );
km.os.group.combine;

save.outlier.figure(
    km.os.group.combine,
    c(i, 'km'),
    width = 7.5,
    height = 7
    );
