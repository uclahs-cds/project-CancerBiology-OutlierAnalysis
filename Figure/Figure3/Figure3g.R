### HISTORY #####################################################################
# This script processes survival data from the TCGA-BRCA and METABRIC datasets
# to analyze overall survival by combining data from both sources. The analysis
# includes Kaplan-Meier survival plots among different patient groups.
# Date: 2024-08-15

library(BoutrosLab.plotting.survival)

### 1. TCGA-BRCA 
brca.outlier.patient.tag.01.t.p.order.sum <- apply(brca.outlier.patient.tag.01.t.p.order, 2, sum);
os.data.brca <- data.frame(cbind(
    status = substr(brca.clinic.order$Overall.Survival.Status, 1, 1),
    os = brca.clinic.order$Overall.Survival..Months.,
    patient = brca.outlier.patient.tag.01.t.p.order.sum, 
    pam50 = brca.clinic.order$Subtype, 
    age = brca.clinic.order$Diagnosis.Age
    ));

# PAM50 subtypes
os.data.brca$pam50[os.data.brca$pam50 == 'BRCA_Basal'] <- 'Basal';
os.data.brca$pam50[os.data.brca$pam50 == 'BRCA_Her2'] <- 'Her2';
os.data.brca$pam50[os.data.brca$pam50 == 'BRCA_LumA'] <- 'LumA';
os.data.brca$pam50[os.data.brca$pam50 == 'BRCA_LumB'] <- 'LumB';
os.data.brca$pam50[os.data.brca$pam50 == 'BRCA_Normal'] <- 'Normal';

os.data.brca[,1] <- as.numeric(os.data.brca[,1]);
os.data.brca[,2] <- as.numeric(os.data.brca[,2]);
os.data.brca[,3] <- as.numeric(os.data.brca[,3]);
os.data.brca[,5] <- as.numeric(os.data.brca[,5]);
os.data.brca <- na.omit(os.data.brca);

# Group patients
os.group.brca <- os.data.brca;
for (i in 1:nrow(os.group.brca)) {
    if (os.group.brca[i,3] == 0) {
        os.group.brca[i,3] <- 1;
        } 
    else {
        os.group.brca[i,3] <- 2;
        }
    }


os.group.brca$pam50 <- factor(os.group.brca$pam50);
os.group.brca$pam50 <- relevel(os.group.brca$pam50, ref = "LumA");
os.group.brca <- os.group.brca[!(os.group.brca$pam50 %in% "NC"),];
os.group.brca <- na.omit(os.group.brca);



### 2. METABRIC Data Preparation
outlier.patient.tag.01.meta.sum <- apply(outlier.patient.tag.01.meta, 2, sum)
os.data.meta <- data.frame(cbind(
    status = substr(meta.clinic.5.order.combine$Overall.Survival.Status, 1, 1),
    os = meta.clinic.5.order.combine$Overall.Survival..Months.,
    patient = outlier.patient.tag.01.meta.sum, 
    pam50 = meta.clinic.5.order.combine$pam50, 
    age = meta.clinic.5.order.combine$Age.at.Diagnosis
));


os.data.meta[,1] <- as.numeric(os.data.meta[,1]);
os.data.meta[,2] <- as.numeric(os.data.meta[,2]);
os.data.meta[,3] <- as.numeric(os.data.meta[,3]);
os.data.meta[,5] <- as.numeric(os.data.meta[,5]);
os.data.meta <- na.omit(os.data.meta);

# Group patients
os.group.meta <- os.data.meta;
for (i in 1:nrow(os.group.meta)) {
    if (os.group.meta[i,3] == 0) {
        os.group.meta[i,3] <- 1;
        } 
    else {
        os.group.meta[i,3] <- 2;
        }
    }


os.group.meta$pam50 <- factor(os.group.meta$pam50);
os.group.meta$pam50 <- relevel(os.group.meta$pam50, ref = "LumA");
os.group.meta <- os.group.meta[!(os.group.meta$pam50 %in% "NC"),];
os.group.meta <- na.omit(os.group.meta);




### 3. Combine TCGA-BRCA and METABRIC Datasets
os.group.combine <- data.frame(rbind(
    os.group.brca,
    os.group.meta
    ));


os.group.combine$pam50 <- factor(os.group.combine$pam50);
os.group.combine$pam50 <- relevel(os.group.combine$pam50, ref = "LumA");
os.group.combine <- os.group.combine[!(os.group.combine$pam50 %in% "NC"),];
os.group.combine <- na.omit(os.group.combine);

### 4. Kaplan-Meier Survival Analysis
km.os.group.combine <- create.km.plot(
    survival.object = Surv(os.group.combine$os, as.numeric(os.group.combine$status)),
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
    key.groups.labels = c("Outlier patients", "Non-outlier patients"),
    key.groups.cex = 1,
    line.colours = rev(c('red3', 'dodgerblue3'))
    );
km.os.group.combine;



# Save the km plot as a PDF
pdf(
    file = generate.filename(
        'os_merge', 
        'km', 
        'pdf'
        ), 
    width = 7.5, 
    height = 7
    );
print(km.os.group.combine);
dev.off();

# Save the km plot as a PNG
png(
    file = generate.filename(
        'os_merge', 
        'km', 
        'png'
        ), 
    width = 7.5, 
    height = 7,
    unit = 'in', 
    res = 1200
    );
print(km.os.group.combine);
dev.off();
