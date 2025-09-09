### HISTORY ######################################################################
# This script performs meta-analysis to compare the tumour stages and subtypes across
# multiple breast cancer subtypes.
# Date: 2024-08-15

### DESCRIPTION ##################################################################
# This script analyzes tumor stage data across different breast cancer subtypes
# using data from all nine datasets. It performs
# Fisher's exact tests, calculates odds ratios, and conducts meta-analysis
# to compare tumor stages across subtypes. The script generates visualizations
# of the results, including dotmaps for each subtype and for all patients.

### PREAMBLE #####################################################################
# Load required libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
library(metafor);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

### 1. TCGA-BRCA
convert.stage <- function(stage) {
    result <- rep(NA, length(stage));
    result[stage %in% c('STAGE I', 'STAGE IA', 'STAGE IB')] <- 1;
    result[stage %in% c('STAGE II', 'STAGE IIA', 'STAGE IIB')] <- 2;
    result[stage %in% c('STAGE III', 'STAGE IIIA', 'STAGE IIIB', 'STAGE IIIC')] <- 3;
    result[stage == 'STAGE IV'] <- 4;
    return(result);
    }

perform.fisher.test.brca <- function(data, subtype = NULL) {
    if (!is.null(subtype)) {
        data <- data[data$pam50 == subtype, ]
        }

    stages <- stages <- c(1, 2, 3);
    p.values <- numeric();
    odd.ratios <- numeric();
    ci.intervals <- matrix(nrow = 0, ncol = 2);

    for (i in stages[-1]) {
        table.1 <- sum(data$patient == 0 & data$stage.convert == stages[1]);
        table.2 <- sum(data$patient == 0 & data$stage.convert == i);
        table.3 <- sum(data$patient > 0 & data$stage.convert == stages[1]);
        table.4 <- sum(data$patient > 0 & data$stage.convert == i);

        fisher.result <- fisher.test(matrix(c(table.1 + 1, table.2 + 1, table.3 + 1, table.4 + 1), nrow = 2), alternative = 'two.sided');

        p.values <- c(p.values, fisher.result$p.value);
        odd.ratios <- c(odd.ratios, fisher.result$estimate);
        ci.intervals <- rbind(ci.intervals, fisher.result$conf.int);
        }

    p.values.fdr <- p.adjust(p.values, method = 'BH');

    return(list(p.values = p.values, odd.ratios = odd.ratios, ci.intervals = ci.intervals, p.values.fdr = p.values.fdr));
    }


outlier.patient.tag.01.brca.sum <- apply(outlier.patient.tag.01.brca, 2, sum);

os.data.stage.pseudo.brca <- data.frame(
    status = substr(brca.clinic.order$Overall.Survival.Status, 1, 1),
    os = brca.clinic.order$Overall.Survival..Months.,
    patient = outlier.patient.tag.01.brca.sum,
    pam50 = brca.clinic.order$Subtype,
    age = brca.clinic.order$Diagnosis.Age,
    stage = brca.clinic.order$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code
    );

os.data.stage.pseudo.brca$pam50 <- as.character(os.data.stage.pseudo.brca$pam50);
os.data.stage.pseudo.brca$pam50[os.data.stage.pseudo.brca$pam50 == 'BRCA_Basal'] <- 'Basal';
os.data.stage.pseudo.brca$pam50[os.data.stage.pseudo.brca$pam50 == 'BRCA_Her2'] <- 'Her2';
os.data.stage.pseudo.brca$pam50[os.data.stage.pseudo.brca$pam50 == 'BRCA_LumA'] <- 'LumA';
os.data.stage.pseudo.brca$pam50[os.data.stage.pseudo.brca$pam50 == 'BRCA_LumB'] <- 'LumB';
os.data.stage.pseudo.brca$pam50[os.data.stage.pseudo.brca$pam50 == 'BRCA_Normal'] <- 'Normal';

os.data.stage.pseudo.brca$stage.convert <- convert.stage(os.data.stage.pseudo.brca$stage);

os.data.stage.pseudo.brca$status <- as.numeric(os.data.stage.pseudo.brca$status);
os.data.stage.pseudo.brca$os <- as.numeric(os.data.stage.pseudo.brca$os);
os.data.stage.pseudo.brca$patient <- as.numeric(os.data.stage.pseudo.brca$patient);
os.data.stage.pseudo.brca$age <- as.numeric(os.data.stage.pseudo.brca$age);
os.data.stage.pseudo.brca <- na.omit(os.data.stage.pseudo.brca);

# All BRCA
all.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca);
p.value.stage.pseudo.odd.sub.t1.brca <- all.results$odd.ratios;
p.value.stage.pseudo.ci.sub.t1.brca <- all.results$ci.intervals;

# Basal
basal.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, 'Basal');
p.value.stage.pseudo.basal.odd.sub.t1.brca <- basal.results$odd.ratios;
p.value.stage.pseudo.basal.ci.sub.t1.brca <- basal.results$ci.intervals;

# Her2
her2.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, 'Her2');
p.value.stage.pseudo.her2.odd.sub.t1.brca <- her2.results$odd.ratios;
p.value.stage.pseudo.her2.ci.sub.t1.brca <- her2.results$ci.intervals;

# LumA
luma.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, 'LumA');
p.value.stage.pseudo.luma.odd.sub.t1.brca <- luma.results$odd.ratios;
p.value.stage.pseudo.luma.ci.sub.t1.brca <- luma.results$ci.intervals;

# LumB
lumb.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, 'LumB');
p.value.stage.pseudo.lumb.odd.sub.t1.brca <- lumb.results$odd.ratios;
p.value.stage.pseudo.lumb.ci.sub.t1.brca <- lumb.results$ci.intervals;

# Normal
normal.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, 'Normal');
p.value.stage.pseudo.normal.odd.sub.t1.brca <- normal.results$odd.ratios;
p.value.stage.pseudo.normal.ci.sub.t1.brca <- normal.results$ci.intervals;




### 2. METABRIC
outlier.patient.tag.01.meta.sum <- apply(outlier.patient.tag.01.meta, 2, sum)


perform.fisher.test.meta <- function(data, subtype = NULL) {
    if (!is.null(subtype)) {
        data <- data[data$pam50 == subtype, ];
        }

    stages <- c(1, 2, 3);
    p.values <- numeric();
    odd.ratios <- numeric();
    ci.intervals <- matrix(nrow = 0, ncol = 2);

    for (i in stages[-1]) {
        table.1 <- sum(data$patient == 0 & data$stage == stages[1]);
        table.2 <- sum(data$patient == 0 & data$stage == i);
        table.3 <- sum(data$patient > 0 & data$stage == stages[1]);
        table.4 <- sum(data$patient > 0 & data$stage == i);

        fisher.result <- fisher.test(matrix(c(table.1 + 1, table.2 + 1, table.3 + 1, table.4 + 1), nrow = 2), alternative = 'two.sided')

        p.values <- c(p.values, fisher.result$p.value);
        odd.ratios <- c(odd.ratios, fisher.result$estimate);
        ci.intervals <- rbind(ci.intervals, fisher.result$conf.int);
        }

    p.values.fdr <- p.adjust(p.values, method = 'BH');

    return(list(p.values = p.values, odd.ratios = odd.ratios, ci.intervals = ci.intervals, p.values.fdr = p.values.fdr));
    }

# Data preparation
os.data.stage.pseudo.meta <- data.frame(
    status = substr(meta.clinic.5.order.combine$Overall.Survival.Status, 1, 1),
    os = meta.clinic.5.order.combine$Overall.Survival..Months.,
    patient = outlier.patient.tag.01.meta.sum,
    pam50 = meta.clinic.5.order.combine$pam50,
    age = meta.clinic.5.order.combine$Age.at.Diagnosis,
    stage = meta.clinic.5.order.combine$stage
    );

os.data.stage.pseudo.meta[, 1] <- as.numeric(os.data.stage.pseudo.meta[, 1]);
os.data.stage.pseudo.meta[, 2] <- as.numeric(os.data.stage.pseudo.meta[, 2]);
os.data.stage.pseudo.meta[, 3] <- as.numeric(os.data.stage.pseudo.meta[, 3]);
os.data.stage.pseudo.meta[, 5] <- as.numeric(os.data.stage.pseudo.meta[, 5]);

os.data.stage.pseudo.meta.stage <- data.frame(table(os.data.stage.pseudo.meta$stage));
os.data.stage.pseudo.meta.stage <- os.data.stage.pseudo.meta.stage[2:5, ];


# Metabric - All subtypes
meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta);
p.value.stage.pseudo.odd.sub.t1.meta <- meta.results$odd.ratios;
p.value.stage.pseudo.ci.sub.t1.meta <- meta.results$ci.intervals;

# Basal subtype
basal.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, 'Basal');
p.value.stage.pseudo.basal.odd.sub.t1.meta <- basal.meta.results$odd.ratios;
p.value.stage.pseudo.basal.ci.sub.t1.meta <- basal.meta.results$ci.intervals;

# Her2 subtype
her2.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, 'Her2');
p.value.stage.pseudo.her2.odd.sub.t1.meta <- her2.meta.results$odd.ratios;
p.value.stage.pseudo.her2.ci.sub.t1.meta <- her2.meta.results$ci.intervals;

# Luminal A subtype
luma.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, 'LumA');
p.value.stage.pseudo.luma.odd.sub.t1.meta <- luma.meta.results$odd.ratios;
p.value.stage.pseudo.luma.ci.sub.t1.meta <- luma.meta.results$ci.intervals;

# Luminal B subtype
lumb.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, 'LumB');
p.value.stage.pseudo.lumb.odd.sub.t1.meta <- lumb.meta.results$odd.ratios;
p.value.stage.pseudo.lumb.ci.sub.t1.meta <- lumb.meta.results$ci.intervals;

# Normal subtype
normal.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, 'Normal');
p.value.stage.pseudo.normal.odd.sub.t1.meta <- normal.meta.results$odd.ratios;
p.value.stage.pseudo.normal.ci.sub.t1.meta <- normal.meta.results$ci.intervals;





# ### 3. ICGC BRCA-EU
icgc.sample.num <- substr(icgc.clinic.order$sample_name, 3, nchar(icgc.clinic.order$sample_name));

# Prepare ICGC sample numbers
icgc.sample.num.nra <- gsub('PR(\\d+)a(\\.RNA|\\.2)?', '\\1', colnames(outlier.patient.tag.01.icgc));
icgc.sample.num.nra <- gsub('PR(\\d+)b(\\.RNA|\\.2)?', '\\1', icgc.sample.num.nra);
icgc.sample.num.nra <- gsub('PR(\\d+)c(\\.RNA|\\.2)?', '\\1', icgc.sample.num.nra);

icgc.clinic.all.order <- icgc.clinic.order[match(icgc.sample.num.nra, icgc.sample.num), ];
rownames(icgc.clinic.all.order) <- colnames(outlier.patient.tag.01.icgc);

# Stage conversion
icgc.stage.convert.new.x <- rep(NA, nrow(icgc.clinic.all.order));
icgc.clinic.all.order.x <- icgc.clinic.all.order;
icgc.clinic.all.order.x$N_stage[icgc.clinic.all.order.x$N_stage %in% 'NX'] <- 'N0';
icgc.clinic.all.order.x$M_stage[icgc.clinic.all.order.x$M_stage %in% c('MX', 'Mx')] <- 'M0';

for (i in 1:nrow(icgc.clinic.all.order.x)) {
    T <- icgc.clinic.all.order.x$T_stage[i];
    N <- icgc.clinic.all.order.x$N_stage[i];
    M <- icgc.clinic.all.order.x$M_stage[i];

    if (is.na(T) | is.na(N) | is.na(M)) {
        next;
        }

    stage <- NA;
    if (T == 'Tis' && N == 'N0' && M == 'M0') {
        stage <- 0;
        } else if ((T == 'T1' && N == 'N0' && M == 'M0') || (T %in% c('T0', 'T1') && N == 'N1mi' && M == 'M0')) {
        stage <- 1;
        } else if ((T == 'T0' && N == 'N1' && M == 'M0') || (T == 'T1' && N == 'N1' && M == 'M0') || (T == 'T2' && N == 'N0' && M == 'M0')) {
        stage <- 2;
        } else if ((T == 'T2' && N == 'N1' && M == 'M0') || (T == 'T3' && N == 'N0' && M == 'M0')) {
        stage <- 2;
        } else if ((T %in% c('T0', 'T1', 'T2', 'T3') && N == 'N2' && M == 'M0') || (T == 'T3' && N == 'N1' && M == 'M0')) {
        stage <- 3;
        } else if (T == 'T4' && (N %in% c('N0', 'N1', 'N2')) && M == 'M0') {
        stage <- 3
        } else if ((T %in% c('T0', 'T1', 'T2', 'T3', 'T4') && N == 'N3' && M == 'M0')) {
        stage <- 3
        } else if (M %in% c('M1', 'M2')) {
        stage <- 4
        } else {
        stage <- NA # Set stage as NA if none of the conditions match
        }

    icgc.stage.convert.new.x[i] <- stage
    }


perform.fisher.test.icgc <- function(data, subtype = NULL) {
    if (!is.null(subtype)) {
        data <- data[data$subtype == subtype, ];
        }

    stages <- c(1, 2, 3);
    p.values <- numeric();
    odd.ratios <- numeric();
    ci.intervals <- matrix(nrow = 0, ncol = 2);

    for (i in stages[-1]) {
        table.1 <- sum(data$outlier == 0 & data$stage == stages[1]);
        table.2 <- sum(data$outlier == 0 & data$stage == i);
        table.3 <- sum(data$outlier > 0 & data$stage == stages[1]);
        table.4 <- sum(data$outlier > 0 & data$stage == i);

        fisher.result <- fisher.test(matrix(c(table.1 + 1, table.2 + 1, table.3 + 1, table.4 + 1), nrow = 2), alternative = 'two.sided');

        p.values <- c(p.values, fisher.result$p.value);
        odd.ratios <- c(odd.ratios, fisher.result$estimate);
        ci.intervals <- rbind(ci.intervals, fisher.result$conf.int);
        }

    p.values.fdr <- p.adjust(p.values, method = 'BH');

    return(list(p.values = p.values, odd.ratios = odd.ratios, ci.intervals = ci.intervals, p.values.fdr = p.values.fdr));
    }


icgc.clinic.order.data <- data.frame(as.character(icgc.clinic.order$subtype));
icgc.clinic.order.data[is.na(icgc.clinic.order.data$as.character.icgc.clinic.order.subtype.), ] <- 6;
icgc.clinic.order.data[icgc.clinic.order.data$as.character.icgc.clinic.order.subtype. == 'Basal', ] <- 1;
icgc.clinic.order.data[icgc.clinic.order.data$as.character.icgc.clinic.order.subtype. == 'Her2', ] <- 2;
icgc.clinic.order.data[icgc.clinic.order.data$as.character.icgc.clinic.order.subtype. == 'LumA', ] <- 3;
icgc.clinic.order.data[icgc.clinic.order.data$as.character.icgc.clinic.order.subtype. == 'LumB', ] <- 4;
icgc.clinic.order.data[icgc.clinic.order.data$as.character.icgc.clinic.order.subtype. == 'Normal', ] <- 5;
icgc.clinic.order.data.num <- data.frame(as.numeric(icgc.clinic.order.data$as.character.icgc.clinic.order.subtype.));

outlier.patient.tag.01.icgc.sum <- apply(outlier.patient.tag.01.icgc, 2, sum);
subtype.total.outlier.num.icgc <- data.frame(cbind(
    subtype = icgc.clinic.order.data.num,
    outlier = outlier.patient.tag.01.icgc.sum
    ));
colnames(subtype.total.outlier.num.icgc) <- c('subtype', 'outlier');

# Data preparation
os.data.stage.pseudo.icgc <- data.frame(cbind(subtype.total.outlier.num.icgc,
    stage = icgc.stage.convert.new.x
    ));
os.data.stage.pseudo.icgc <- na.omit(os.data.stage.pseudo.icgc);

# All ICGC
all.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc);
p.value.stage.pseudo.odd.sub.t1.icgc <- all.results.icgc$odd.ratios;
p.value.stage.pseudo.ci.sub.t1.icgc <- all.results.icgc$ci.intervals;

# Basal
basal.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 1);
p.value.stage.pseudo.basal.odd.sub.t1.icgc <- basal.results.icgc$odd.ratios;
p.value.stage.pseudo.basal.ci.sub.t1.icgc <- basal.results.icgc$ci.intervals;

# Lum A
luma.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 3);
p.value.stage.pseudo.luma.odd.sub.t1.icgc <- luma.results.icgc$odd.ratios;
p.value.stage.pseudo.luma.ci.sub.t1.icgc <- luma.results.icgc$ci.intervals;

# Lum B
lumb.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 4);
p.value.stage.pseudo.lumb.odd.sub.t1.icgc <- lumb.results.icgc$odd.ratios;
p.value.stage.pseudo.lumb.ci.sub.t1.icgc <- lumb.results.icgc$ci.intervals;

# Her2
her2.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 2);
p.value.stage.pseudo.her2.odd.sub.t1.icgc <- her2.results.icgc$odd.ratios;
p.value.stage.pseudo.her2.ci.sub.t1.icgc <- her2.results.icgc$ci.intervals;

# Normal
normal.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 5);
p.value.stage.pseudo.normal.odd.sub.t1.icgc <- normal.results.icgc$odd.ratios;
p.value.stage.pseudo.normal.ci.sub.t1.icgc <- normal.results.icgc$ci.intervals;




### 4. cheng
outlier.patient.tag.01.cheng.sum <- apply(outlier.patient.tag.01.cheng, 2, sum)


perform.fisher.test.cheng <- function(data, subtype = NULL) {
    if (!is.null(subtype)) {
        data <- data[data$pam50 == subtype, ];
        }

    stages <- c(1, 2, 3);
    p.values <- numeric();
    odd.ratios <- numeric();
    ci.intervals <- matrix(nrow = 0, ncol = 2);

    for (i in stages[-1]) {
        table.1 <- sum(data$patient == 0 & data$stage == stages[1]);
        table.2 <- sum(data$patient == 0 & data$stage == i);
        table.3 <- sum(data$patient > 0 & data$stage == stages[1]);
        table.4 <- sum(data$patient > 0 & data$stage == i);

        fisher.result <- fisher.test(matrix(c(table.1 + 1, table.2 + 1, table.3 + 1, table.4 + 1), nrow = 2), alternative = 'two.sided')

        p.values <- c(p.values, fisher.result$p.value);
        odd.ratios <- c(odd.ratios, fisher.result$estimate);
        ci.intervals <- rbind(ci.intervals, fisher.result$conf.int);
        }

    p.values.fdr <- p.adjust(p.values, method = 'BH');

    return(list(p.values = p.values, odd.ratios = odd.ratios, ci.intervals = ci.intervals, p.values.fdr = p.values.fdr));
    }

# Data preparation
os.data.stage.pseudo.cheng <- data.frame(
    patient = outlier.patient.tag.01.cheng.sum,
    pam50 = patient.cheng$subtype.genefu,
    stage = patient.cheng$stage
    );

os.data.stage.pseudo.cheng[, 1] <- as.numeric(os.data.stage.pseudo.cheng[, 1]);
# os.data.stage.pseudo.cheng[, 2] <- as.numeric(os.data.stage.pseudo.cheng[, 2]);
os.data.stage.pseudo.cheng[, 3] <- as.numeric(os.data.stage.pseudo.cheng[, 3]);

os.data.stage.pseudo.cheng.stage <- data.frame(table(os.data.stage.pseudo.cheng$stage));


# chengbric - All subtypes
cheng.results <- perform.fisher.test.cheng(os.data.stage.pseudo.cheng);
p.value.stage.pseudo.odd.sub.t1.cheng <- cheng.results$odd.ratios;
p.value.stage.pseudo.ci.sub.t1.cheng <- cheng.results$ci.intervals;

# Basal subtype
basal.cheng.results <- perform.fisher.test.cheng(os.data.stage.pseudo.cheng, 'Basal');
p.value.stage.pseudo.basal.odd.sub.t1.cheng <- basal.cheng.results$odd.ratios;
p.value.stage.pseudo.basal.ci.sub.t1.cheng <- basal.cheng.results$ci.intervals;

# Her2 subtype
her2.cheng.results <- perform.fisher.test.cheng(os.data.stage.pseudo.cheng, 'Her2');
p.value.stage.pseudo.her2.odd.sub.t1.cheng <- her2.cheng.results$odd.ratios;
p.value.stage.pseudo.her2.ci.sub.t1.cheng <- her2.cheng.results$ci.intervals;

# Luminal A subtype
luma.cheng.results <- perform.fisher.test.cheng(os.data.stage.pseudo.cheng, 'LumA');
p.value.stage.pseudo.luma.odd.sub.t1.cheng <- luma.cheng.results$odd.ratios;
p.value.stage.pseudo.luma.ci.sub.t1.cheng <- luma.cheng.results$ci.intervals;

# Luminal B subtype
lumb.cheng.results <- perform.fisher.test.cheng(os.data.stage.pseudo.cheng, 'LumB');
p.value.stage.pseudo.lumb.odd.sub.t1.cheng <- lumb.cheng.results$odd.ratios;
p.value.stage.pseudo.lumb.ci.sub.t1.cheng <- lumb.cheng.results$ci.intervals;

# Normal subtype
normal.cheng.results <- perform.fisher.test.cheng(os.data.stage.pseudo.cheng, 'Normal');
p.value.stage.pseudo.normal.odd.sub.t1.cheng <- normal.cheng.results$odd.ratios;
p.value.stage.pseudo.normal.ci.sub.t1.cheng <- normal.cheng.results$ci.intervals;


### 5. hatzis
outlier.patient.tag.01.hatzis.sum <- apply(outlier.patient.tag.01.hatzis, 2, sum)


perform.fisher.test.hatzis <- function(data, subtype = NULL) {
    if (!is.null(subtype)) {
        data <- data[data$pam50 == subtype, ];
        }

    data <- na.omit(data);
    stages <- c(1, 2, 3);
    p.values <- numeric();
    odd.ratios <- numeric();
    ci.intervals <- matrix(nrow = 0, ncol = 2);

    for (i in stages[-1]) {
        table.1 <- sum(data$patient == 0 & data$stage == stages[1]);
        table.2 <- sum(data$patient == 0 & data$stage == i);
        table.3 <- sum(data$patient > 0 & data$stage == stages[1]);
        table.4 <- sum(data$patient > 0 & data$stage == i);
        
        table.1[is.na(table.1)] <- 0;
        table.2[is.na(table.2)] <- 0;
        table.3[is.na(table.3)] <- 0;
        table.4[is.na(table.4)] <- 0;

        fisher.result <- fisher.test(matrix(c(table.1 + 1, table.2 + 1, table.3 + 1, table.4 + 1), nrow = 2), alternative = 'two.sided')

        p.values <- c(p.values, fisher.result$p.value);
        odd.ratios <- c(odd.ratios, fisher.result$estimate);
        ci.intervals <- rbind(ci.intervals, fisher.result$conf.int);
        }

    p.values.fdr <- p.adjust(p.values, method = 'BH');

    return(list(p.values = p.values, odd.ratios = odd.ratios, ci.intervals = ci.intervals, p.values.fdr = p.values.fdr));
    }

# Data preparation
os.data.stage.pseudo.hatzis <- data.frame(
    patient = outlier.patient.tag.01.hatzis.sum,
    pam50 = patient.hatzis$subtype,
    stage = substr(patient.hatzis$stage, 1, 1)
    );

os.data.stage.pseudo.hatzis[, 1] <- as.numeric(os.data.stage.pseudo.hatzis[, 1]);
# os.data.stage.pseudo.hatzis[, 2] <- as.numeric(os.data.stage.pseudo.hatzis[, 2]);
os.data.stage.pseudo.hatzis[, 3] <- as.numeric(os.data.stage.pseudo.hatzis[, 3]);

os.data.stage.pseudo.hatzis.stage <- data.frame(table(os.data.stage.pseudo.hatzis$stage));


# hatzisbric - All subtypes
hatzis.results <- perform.fisher.test.hatzis(os.data.stage.pseudo.hatzis);
p.value.stage.pseudo.odd.sub.t1.hatzis <- hatzis.results$odd.ratios;
p.value.stage.pseudo.ci.sub.t1.hatzis <- hatzis.results$ci.intervals;

# Basal subtype
basal.hatzis.results <- perform.fisher.test.hatzis(os.data.stage.pseudo.hatzis, 'Basal');
p.value.stage.pseudo.basal.odd.sub.t1.hatzis <- basal.hatzis.results$odd.ratios;
p.value.stage.pseudo.basal.ci.sub.t1.hatzis <- basal.hatzis.results$ci.intervals;

# Her2 subtype
her2.hatzis.results <- perform.fisher.test.hatzis(os.data.stage.pseudo.hatzis, 'Her2');
p.value.stage.pseudo.her2.odd.sub.t1.hatzis <- her2.hatzis.results$odd.ratios;
p.value.stage.pseudo.her2.ci.sub.t1.hatzis <- her2.hatzis.results$ci.intervals;

# Luminal A subtype
luma.hatzis.results <- perform.fisher.test.hatzis(os.data.stage.pseudo.hatzis, 'LumA');
p.value.stage.pseudo.luma.odd.sub.t1.hatzis <- luma.hatzis.results$odd.ratios;
p.value.stage.pseudo.luma.ci.sub.t1.hatzis <- luma.hatzis.results$ci.intervals;

# Luminal B subtype
lumb.hatzis.results <- perform.fisher.test.hatzis(os.data.stage.pseudo.hatzis, 'LumB');
p.value.stage.pseudo.lumb.odd.sub.t1.hatzis <- lumb.hatzis.results$odd.ratios;
p.value.stage.pseudo.lumb.ci.sub.t1.hatzis <- lumb.hatzis.results$ci.intervals;

# Normal subtype
normal.hatzis.results <- perform.fisher.test.hatzis(os.data.stage.pseudo.hatzis, 'Normal');
p.value.stage.pseudo.normal.odd.sub.t1.hatzis <- normal.hatzis.results$odd.ratios;
p.value.stage.pseudo.normal.ci.sub.t1.hatzis <- normal.hatzis.results$ci.intervals;



data.frame(table(os.data.stage.pseudo.hatzis[os.data.stage.pseudo.hatzis$pam50 == 'Basal',][,c(3)]));
data.frame(table(os.data.stage.pseudo.hatzis[os.data.stage.pseudo.hatzis$pam50 == 'Her2',][,c(3)]));
data.frame(table(os.data.stage.pseudo.hatzis[os.data.stage.pseudo.hatzis$pam50 == 'LumA',][,c(3)]));
data.frame(table(os.data.stage.pseudo.hatzis[os.data.stage.pseudo.hatzis$pam50 == 'LumB',][,c(3)]));
data.frame(table(os.data.stage.pseudo.hatzis[os.data.stage.pseudo.hatzis$pam50 == 'Normal',][,c(3)]));



data.frame(table(os.data.stage.pseudo.cheng[os.data.stage.pseudo.cheng$pam50 == 'Basal',][,c(3)]));
data.frame(table(os.data.stage.pseudo.cheng[os.data.stage.pseudo.cheng$pam50 == 'Her2',][,c(3)]));
data.frame(table(os.data.stage.pseudo.cheng[os.data.stage.pseudo.cheng$pam50 == 'LumA',][,c(3)]));
data.frame(table(os.data.stage.pseudo.cheng[os.data.stage.pseudo.cheng$pam50 == 'LumB',][,c(3)]));
data.frame(table(os.data.stage.pseudo.cheng[os.data.stage.pseudo.cheng$pam50 == 'Normal',][,c(3)]));




data.frame(table(os.data.stage.pseudo.brca[os.data.stage.pseudo.brca$pam50 == 'Basal',][,c(6)]));
data.frame(table(os.data.stage.pseudo.brca[os.data.stage.pseudo.brca$pam50 == 'Her2',][,c(6)]));
data.frame(table(os.data.stage.pseudo.brca[os.data.stage.pseudo.brca$pam50 == 'LumA',][,c(6)]));
data.frame(table(os.data.stage.pseudo.brca[os.data.stage.pseudo.brca$pam50 == 'LumB',][,c(6)]));
data.frame(table(os.data.stage.pseudo.brca[os.data.stage.pseudo.brca$pam50 == 'Normal',][,c(6)]));


data.frame(table(os.data.stage.pseudo.meta[os.data.stage.pseudo.meta$pam50 == 'Basal',][,c(6)]));
data.frame(table(os.data.stage.pseudo.meta[os.data.stage.pseudo.meta$pam50 == 'Her2',][,c(6)]));
data.frame(table(os.data.stage.pseudo.meta[os.data.stage.pseudo.meta$pam50 == 'LumA',][,c(6)]));
data.frame(table(os.data.stage.pseudo.meta[os.data.stage.pseudo.meta$pam50 == 'LumB',][,c(6)]));
data.frame(table(os.data.stage.pseudo.meta[os.data.stage.pseudo.meta$pam50 == 'Normal',][,c(6)]));


data.frame(table(os.data.stage.pseudo.icgc[os.data.stage.pseudo.icgc$subtype == 1,][,c(3)]));
data.frame(table(os.data.stage.pseudo.icgc[os.data.stage.pseudo.icgc$subtype == 2,][,c(3)]));
data.frame(table(os.data.stage.pseudo.icgc[os.data.stage.pseudo.icgc$subtype == 3,][,c(3)]));
data.frame(table(os.data.stage.pseudo.icgc[os.data.stage.pseudo.icgc$subtype == 4,][,c(3)]));
data.frame(table(os.data.stage.pseudo.icgc[os.data.stage.pseudo.icgc$subtype == 5,][,c(3)]));




### Meta-analysis
# Calculate log odds and standard error
calculate.ln.odd.se <- function(odd, ci) {
    ln.odd <- log(odd);
    se.odd <- (log(ci[, 2]) - log(ci[, 1])) / 3.92;
    return(list(ln.odd = ln.odd, se.odd = se.odd));
    }

# Perform meta-analysis using metafor
perform.meta.analysis <- function(ln.odd, se.odd) {
    stage.pseudo.odd.se.t1 <- list();

    for (i in 1:length(ln.odd[[1]])) {
        chr.odd <- sapply(ln.odd, function(x) x[i]);
        chr.se <- sapply(se.odd, function(x) x[i]);
        chr.all <- data.frame(chr.odd, chr.se);
        stage.pseudo.odd.se.t1[[i]] <- chr.all;
        }

    metafor.stage.pseudo.odd.ci.p.t1 <- NULL;

    for (i in 1:length(stage.pseudo.odd.se.t1)) {
        chr.odd.se.sample <- stage.pseudo.odd.se.t1[[i]];
        chr.odd.se.sample.inf <- chr.odd.se.sample[!is.infinite(chr.odd.se.sample$chr.odd) & !is.infinite(chr.odd.se.sample$chr.se), ];

        if (nrow(chr.odd.se.sample.inf) > 1) {
            metafor.chr <- rma.uni(yi = chr.odd, sei = chr.se, data = chr.odd.se.sample.inf, method = 'DL');
            metafor.chr.odd <- exp(metafor.chr$beta);
            metafor.chr.lower <- exp(metafor.chr$ci.lb);
            metafor.chr.upper <- exp(metafor.chr$ci.ub);
            metafor.chr.p <- metafor.chr$pval;
            metafor.all <- c(metafor.chr.odd, metafor.chr.lower, metafor.chr.upper, metafor.chr.p);
            } else {
            metafor.all <- c(NA, NA, NA, NA);
            }

        metafor.stage.pseudo.odd.ci.p.t1 <- rbind(metafor.stage.pseudo.odd.ci.p.t1, metafor.all);
        }

    metafor.stage.pseudo.odd.ci.p.data.t1 <- data.frame(
        p.value = metafor.stage.pseudo.odd.ci.p.t1[, 4],
        odd = metafor.stage.pseudo.odd.ci.p.t1[, 1],
        ci.min = metafor.stage.pseudo.odd.ci.p.t1[, 2],
        ci.max = metafor.stage.pseudo.odd.ci.p.t1[, 3]
        );

    metafor.stage.pseudo.odd.ci.p.data.fdr.t1 <- p.adjust(metafor.stage.pseudo.odd.ci.p.data.t1$p.value, method = 'BH');

    return(list(
        data = metafor.stage.pseudo.odd.ci.p.data.t1,
        fdr = metafor.stage.pseudo.odd.ci.p.data.fdr.t1
        ));
    }


# 1. Basal
ln.odd.brca.basal <- calculate.ln.odd.se(p.value.stage.pseudo.basal.odd.sub.t1.brca, p.value.stage.pseudo.basal.ci.sub.t1.brca);
ln.odd.meta.basal <- calculate.ln.odd.se(p.value.stage.pseudo.basal.odd.sub.t1.meta, p.value.stage.pseudo.basal.ci.sub.t1.meta);
ln.odd.icgc.basal <- calculate.ln.odd.se(p.value.stage.pseudo.basal.odd.sub.t1.icgc, p.value.stage.pseudo.basal.ci.sub.t1.icgc);
ln.odd.cheng.basal <- calculate.ln.odd.se(p.value.stage.pseudo.basal.odd.sub.t1.cheng, p.value.stage.pseudo.basal.ci.sub.t1.cheng);
ln.odd.hatzis.basal <- calculate.ln.odd.se(p.value.stage.pseudo.basal.odd.sub.t1.hatzis, p.value.stage.pseudo.basal.ci.sub.t1.hatzis);

basal.results <- perform.meta.analysis(
    list(ln.odd.brca.basal$ln.odd, ln.odd.meta.basal$ln.odd, ln.odd.icgc.basal$ln.odd, ln.odd.cheng.basal$ln.odd, ln.odd.hatzis.basal$ln.odd),
    list(ln.odd.brca.basal$se.odd, ln.odd.meta.basal$se.odd, ln.odd.icgc.basal$se.odd, ln.odd.cheng.basal$se.odd, ln.odd.hatzis.basal$se.odd)
    );

# 2. Her2 (excluding hatzis)
ln.odd.brca.her2 <- calculate.ln.odd.se(p.value.stage.pseudo.her2.odd.sub.t1.brca, p.value.stage.pseudo.her2.ci.sub.t1.brca);
ln.odd.meta.her2 <- calculate.ln.odd.se(p.value.stage.pseudo.her2.odd.sub.t1.meta, p.value.stage.pseudo.her2.ci.sub.t1.meta);
ln.odd.cheng.her2 <- calculate.ln.odd.se(p.value.stage.pseudo.her2.odd.sub.t1.cheng, p.value.stage.pseudo.her2.ci.sub.t1.cheng);
ln.odd.icgc.her2 <- calculate.ln.odd.se(p.value.stage.pseudo.her2.odd.sub.t1.icgc, p.value.stage.pseudo.her2.ci.sub.t1.icgc);

her2.results <- perform.meta.analysis(
    list(ln.odd.brca.her2$ln.odd, ln.odd.meta.her2$ln.odd, ln.odd.cheng.her2$ln.odd, ln.odd.icgc.her2$ln.odd),
    list(ln.odd.brca.her2$se.odd, ln.odd.meta.her2$se.odd, ln.odd.cheng.her2$se.odd, ln.odd.icgc.her2$se.odd)
    );

# 3. LumA 
ln.odd.brca.luma <- calculate.ln.odd.se(p.value.stage.pseudo.luma.odd.sub.t1.brca, p.value.stage.pseudo.luma.ci.sub.t1.brca);
ln.odd.meta.luma <- calculate.ln.odd.se(p.value.stage.pseudo.luma.odd.sub.t1.meta, p.value.stage.pseudo.luma.ci.sub.t1.meta);
ln.odd.icgc.luma <- calculate.ln.odd.se(p.value.stage.pseudo.luma.odd.sub.t1.icgc, p.value.stage.pseudo.luma.ci.sub.t1.icgc);
ln.odd.cheng.luma <- calculate.ln.odd.se(p.value.stage.pseudo.luma.odd.sub.t1.cheng, p.value.stage.pseudo.luma.ci.sub.t1.cheng);
ln.odd.hatzis.luma <- calculate.ln.odd.se(p.value.stage.pseudo.luma.odd.sub.t1.hatzis, p.value.stage.pseudo.luma.ci.sub.t1.hatzis);

luma.results <- perform.meta.analysis(
    list(ln.odd.brca.luma$ln.odd, ln.odd.meta.luma$ln.odd, ln.odd.icgc.luma$ln.odd, ln.odd.cheng.luma$ln.odd, ln.odd.hatzis.luma$ln.odd),
    list(ln.odd.brca.luma$se.odd, ln.odd.meta.luma$se.odd, ln.odd.icgc.luma$se.odd, ln.odd.cheng.luma$se.odd, ln.odd.hatzis.luma$se.odd)
    );

# 4. LumB (excluding hatzis)
ln.odd.brca.lumb <- calculate.ln.odd.se(p.value.stage.pseudo.lumb.odd.sub.t1.brca, p.value.stage.pseudo.lumb.ci.sub.t1.brca);
ln.odd.meta.lumb <- calculate.ln.odd.se(p.value.stage.pseudo.lumb.odd.sub.t1.meta, p.value.stage.pseudo.lumb.ci.sub.t1.meta);
ln.odd.icgc.lumb <- calculate.ln.odd.se(p.value.stage.pseudo.lumb.odd.sub.t1.icgc, p.value.stage.pseudo.lumb.ci.sub.t1.icgc);
ln.odd.cheng.lumb <- calculate.ln.odd.se(p.value.stage.pseudo.lumb.odd.sub.t1.cheng, p.value.stage.pseudo.lumb.ci.sub.t1.cheng);

lumb.results <- perform.meta.analysis(
    list(ln.odd.brca.lumb$ln.odd, ln.odd.meta.lumb$ln.odd, ln.odd.icgc.lumb$ln.odd, ln.odd.cheng.lumb$ln.odd),
    list(ln.odd.brca.lumb$se.odd, ln.odd.meta.lumb$se.odd, ln.odd.icgc.lumb$se.odd, ln.odd.cheng.lumb$se.odd)
    );

# 5. Normal (excluding ICGC)
ln.odd.brca.normal <- calculate.ln.odd.se(p.value.stage.pseudo.normal.odd.sub.t1.brca, p.value.stage.pseudo.normal.ci.sub.t1.brca);
ln.odd.meta.normal <- calculate.ln.odd.se(p.value.stage.pseudo.normal.odd.sub.t1.meta, p.value.stage.pseudo.normal.ci.sub.t1.meta);
ln.odd.cheng.normal <- calculate.ln.odd.se(p.value.stage.pseudo.normal.odd.sub.t1.cheng, p.value.stage.pseudo.normal.ci.sub.t1.cheng);
ln.odd.hatzis.normal <- calculate.ln.odd.se(p.value.stage.pseudo.normal.odd.sub.t1.hatzis, p.value.stage.pseudo.normal.ci.sub.t1.hatzis);

normal.results <- perform.meta.analysis(
    list(ln.odd.brca.normal$ln.odd, ln.odd.meta.normal$ln.odd, ln.odd.cheng.normal$ln.odd, ln.odd.hatzis.normal$ln.odd),
    list(ln.odd.brca.normal$se.odd, ln.odd.meta.normal$se.odd, ln.odd.cheng.normal$se.odd, ln.odd.hatzis.normal$se.odd)
    );


# 6. All datasets
ln.odd.brca <- calculate.ln.odd.se(p.value.stage.pseudo.odd.sub.t1.brca, p.value.stage.pseudo.ci.sub.t1.brca);
ln.odd.meta <- calculate.ln.odd.se(p.value.stage.pseudo.odd.sub.t1.meta, p.value.stage.pseudo.ci.sub.t1.meta);
ln.odd.icgc <- calculate.ln.odd.se(p.value.stage.pseudo.odd.sub.t1.icgc, p.value.stage.pseudo.ci.sub.t1.icgc);
ln.odd.cheng <- calculate.ln.odd.se(p.value.stage.pseudo.odd.sub.t1.cheng, p.value.stage.pseudo.ci.sub.t1.cheng);
ln.odd.hatzis <- calculate.ln.odd.se(p.value.stage.pseudo.odd.sub.t1.hatzis, p.value.stage.pseudo.ci.sub.t1.hatzis);

all.results <- perform.meta.analysis(
    list(ln.odd.brca$ln.odd, ln.odd.meta$ln.odd, ln.odd.icgc$ln.odd, ln.odd.cheng$ln.odd, ln.odd.hatzis$ln.odd),
    list(ln.odd.brca$se.odd, ln.odd.meta$se.odd, ln.odd.icgc$se.odd, ln.odd.cheng$se.odd, ln.odd.hatzis$se.odd)
    );
# Combine the results for heatmap generation
metafor.stage.pseudo.odd.ci.p.data.each.subtype.odd.t1.no3 <- data.frame(
    basal = basal.results$data$odd[1:2],
    her2 = her2.results$data$odd[1:2],
    luma = luma.results$data$odd[1:2],
    lumb = lumb.results$data$odd[1:2],
    normal = normal.results$data$odd[1:2]
    );


metafor.stage.pseudo.odd.ci.p.data.each.subtype.p.t1.no3 <- data.frame(
    basal = basal.results$data$p.value[1:2],
    her2 = her2.results$data$p.value[1:2],
    luma = luma.results$data$p.value[1:2],
    lumb = lumb.results$data$p.value[1:2],
    normal = normal.results$data$p.value[1:2]
    );

metafor.stage.pseudo.odd.ci.p.data.each.subtype.fdr.t1.no3 <- data.frame(
    basal = basal.results$fdr[1:2],
    her2 = her2.results$fdr[1:2],
    luma = luma.results$fdr[1:2],
    lumb = lumb.results$fdr[1:2],
    normal = normal.results$fdr[1:2]
    );

metafor.stage.pseudo.odd.ci.p.data.t1 <- data.frame(
    odd = all.results$data$odd[1:2],
    p.value = all.results$data$p.value[1:2]
    );


# Plotting the results
background.cutoff <- 2;

colourkey.labels.at <- seq(0, background.cutoff, by = 2);
colourkey.labels <- sapply(
    X = colourkey.labels.at,
    FUN = function(x) {
        if (x == 0) {
            return(expression('10'^'0'));
            } else if (x != background.cutoff) {
            return(as.expression(bquote('10'^-.(as.character(x)))));
            } else {
            return(as.expression(bquote('<10'^-.(as.character(x)))));
            }
        }
    );

legend <- legend.grob(
    list(
        legend = list(
            title = expression(underline('P-value')),
            continuous = TRUE,
            colours = c('white', 'black'),
            total.colours = 100,
            labels = colourkey.labels,
            cex = 0.9,
            at = seq(0, 100, length.out = length(colourkey.labels)),
            height = 3
            )
        ),
    label.cex = 1,
    title.cex = 1,
    title.just = 'left',
    title.fontface = 'plain',
    between.row = 4
    );

# Spot size and color functions
spot.size.function <- function(x) {
    1 + (1.5 * abs(x));
    }
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[sign(x) == -1] <- default.colours(2, palette.type = 'dotmap')[1];
    colours[sign(x) == 1] <- default.colours(2, palette.type = 'dotmap')[2];
    return(colours);
    }

# Create dotmap for each subtype
dot.each <- create.dotmap(
    x = t(log2(metafor.stage.pseudo.odd.ci.p.data.each.subtype.odd.t1.no3)),
    main = expression('Association with the tumour stage'),
    main.cex = 1.4,
    yaxis.cex = 1.2,
    xaxis.rot = 90,
    xaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.lab = c('Basal', 'Her2', 'LumA', 'LumB', 'Normal'),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    text = list(
                        lab = expression('Odds Ratio'),
                        cex = 1
                        ),
                    padding.text = 4.5
                    )
                ),
            x = 1,
            y = 1
            ),
        inside = list(fun = legend, x = 1.07, y = 0)
        ),
    key = list(
        space = 'right',
        points = list(
            cex = spot.size.function(seq(-2, 2, 1)),
            col = spot.colour.function(seq(-2, 2, 1)),
            pch = 19
            ),
        text = list(
            lab = c('0.25', '0.5', '1', '2', '4'),
            # lab = c('0.5', '1', '2'),
            cex = 1,
            adj = 1,
            fontface = 'bold'
            ),
        padding.text = 8
        ),
    key.top = 1,
    right.padding = 2,
    pch = 21,
    pch.border.col = 'white',
    bg.data = t(-log10(metafor.stage.pseudo.odd.ci.p.data.each.subtype.p.t1.no3)),
    colourkey = FALSE,
    colourkey.cex = 0.95,
    bg.alpha = 1,
    colour.scheme = c('white', 'black'),
    at = seq(0, 2, 0.01),
    row.colour = 'white',
    col.colour = 'white',
    row.lwd = 1.3,
    col.lwd = 1.3
    );

# Create dotmap for all patients
dot.all <- create.dotmap(
    x = t(log2(metafor.stage.pseudo.odd.ci.p.data.t1$odd[1:2])),
    xaxis.lab = c('Stage II', 'Stage III'),
    yaxis.cex = 1.2,
    top.padding = 10,
    xaxis.rot = 90,
    xaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.lab = c('All patients'),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    key.top = 1,
    right.padding = 2,
    pch = 21,
    pch.border.col = 'white',
    # bg.data = t(-log10(p.adjust(metafor.stage.pseudo.odd.ci.p.data.t1[1:2, ]$p.value, method = 'BH'))),
    bg.data = t(-log10(metafor.stage.pseudo.odd.ci.p.data.t1[1:2, ]$p.value)),
    colourkey = FALSE,
    bg.alpha = 1,
    colour.scheme = c('white', 'black'),
    at = seq(0, 2, 0.01),
    row.colour = 'white',
    col.colour = 'white',
    row.lwd = 1.3,
    col.lwd = 1.3
    );

# Combine both dotmaps into a multipanel plot
dot.multi <- create.multipanelplot(
    list(dot.each, dot.all),
    layout.height = 2,
    layout.width = 1,
    layout.skip = c(FALSE, FALSE),
    plot.objects.heights = c(9.5, 5),
    x.spacing = -1,
    y.spacing = -10,
    bottom.padding = 0,
    top.padding = 2,
    right.padding = 0
    );


save.outlier.figure(
    dot.multi,
    c('Figure3f', 'tumour', 'stage', 'multipanel'),
    width = 4.8,
    height = 5.6
    );







# 1. TCGA-BRCA
# brca.clinic.order.origin <- brca.clinic.order;
brca.clinic.order <- brca.clinic.order[match(substr(colnames(outlier.patient.tag.01.brca), 1, 12), rownames(brca.clinic.order)),]
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
subtype.5.meta.status <- data.frame(table(meta.clinic.5.order.combine$pam50));
subtype.5.meta.status <- subtype.5.meta.status[match(c('Basal', 'Her2', 'LumA', 'LumB', 'Normal'), subtype.5.meta.status$Var1), ];

# METABRIC individual contingency table
### 1. subtype
outlier.num.subtype.5.patient.table.meta <- table(subtype.5.total.outlier.num.meta)[1:5,];
rownames(outlier.num.subtype.5.patient.table.meta) <- c('Basal', 'Her2', 'LuminalA', 'LuminalB', 'Normal');
outlier.num.subtype.5.patient.table.meta.5 <- cbind(outlier.num.subtype.5.patient.table.meta[,1:3],
                                             apply(outlier.num.subtype.5.patient.table.meta[,4:ncol(outlier.num.subtype.5.patient.table.meta)], 1, sum));
colnames(outlier.num.subtype.5.patient.table.meta.5) <- c(0:3);
outlier.num.subtype.5.patient.table.meta.5 <- as.table(outlier.num.subtype.5.patient.table.meta.5);
ylabel <- "Subtype";
xlabel <- "Number of outliers per patient";
outlier.num.subtype.5.patient.table.meta.5.only <- outlier.num.subtype.5.patient.table.meta.5[,1:4];
outlier.num.subtype.5.patient.table.meta.5.only <- outlier.num.subtype.5.patient.table.meta.5.only[c('Basal', 'Her2', 'LuminalB', 'Normal', 'LuminalA'),]

# Chi-square test
expected.p.meta <- prop.table(colSums(outlier.num.subtype.5.patient.table.meta.5.only));
expected.table.meta <- outer(as.vector(rowSums(outlier.num.subtype.5.patient.table.meta.5.only)), expected.p.meta);

# Chi-square test
expected.p.meta <- prop.table(colSums(outlier.num.subtype.5.patient.table.meta.5.only));
expected.table.meta <- outer(as.vector(rowSums(outlier.num.subtype.5.patient.table.meta.5.only)), expected.p.meta);

chisq.test(outlier.num.subtype.5.patient.table.meta.5.only)
p.value.outlier.subtype.5.chi.meta <- NULL;
for (i in 1:5) {
    p.chi <- chisq.test(outlier.num.subtype.5.patient.table.meta.5.only[i,], p = expected.p.meta)$p.value;
    p.value.outlier.subtype.5.chi.meta <- c(p.value.outlier.subtype.5.chi.meta, p.chi);
    }
p.value.outlier.subtype.5.chi.fdr.meta <- p.adjust(p.value.outlier.subtype.5.chi.meta, method = 'BH');


p.value.outlier.subtype.5.chi.sub.meta <- NULL;
for (i in 1:5) {
    expected.p <- prop.table(colSums(outlier.num.subtype.5.patient.table.meta.5.only[-i,]));
    p.chi <- chisq.test(outlier.num.subtype.5.patient.table.meta.5.only[i,], p = expected.p)$p.value;
    p.value.outlier.subtype.5.chi.sub.meta <- c(p.value.outlier.subtype.5.chi.sub.meta, p.chi);
    }


text.size <- 1;
main.heatmap2 <- create.heatmap(
	x = outlier.num.subtype.5.patient.table.meta.5.only / sum(outlier.num.subtype.5.patient.table.meta.5.only),
	clustering = 'none',
	colour.scheme = c('white', 'dodgerblue'),
	colour.alpha = 1,
	at = seq(0, 1, 0.1),
	cell.text = data.frame(outlier.num.subtype.5.patient.table.meta.5.only)$Freq,
	text.cex = 1.2,
	text.fontface = 1,
	xaxis.fontface = 1,
	yaxis.fontface = 1,
	grid.row = TRUE,
	grid.col = TRUE,
	ylab.label = expression("Subtype"),
	xlab.label = expression("Number of outliers per patient"),
	yaxis.lab = rownames(outlier.num.subtype.5.patient.table.meta.5.only),
	xaxis.tck = 0,
	yaxis.tck = 0,
	# xaxis.lab = colnames(plot.data.all),
	xaxis.cex = 1.1,
	yaxis.cex = 1.1,
    ylab.cex = 1.5,
	xlab.cex = 1.5,
	col.pos = rep(1:ncol(outlier.num.subtype.5.patient.table.meta.5.only), each = nrow(outlier.num.subtype.5.patient.table.meta.5.only)),
	row.pos = rep(nrow(outlier.num.subtype.5.patient.table.meta.5.only):1, times = ncol(outlier.num.subtype.5.patient.table.meta.5.only)),
	print.colour.key = FALSE,
	same.as.matrix = TRUE,
	yaxis.rot = 90,
	# xaxis.rot = 0,
	use.legacy.settings = FALSE,
	x.alternating = 2,
	xaxis.rot.top = 0,
	xaxis.lab.top = c(colnames(outlier.num.subtype.5.patient.table.meta.5.only)[1:3], "more than 3"),
	xlab.top.cex = 1.3,
	xlab.top.y = 2.4
	);
# row total heatmap
row.total <- create.heatmap(
	x = data.frame(rowSums(outlier.num.subtype.5.patient.table.meta.5.only), rowSums(outlier.num.subtype.5.patient.table.meta.5.only), rowSums(outlier.num.subtype.5.patient.table.meta.5.only)) / sum(outlier.num.subtype.5.patient.table.meta.5.only),
	clustering = 'none',
	colour.scheme = c('white', 'dodgerblue'),
	colour.alpha = 1,
	at = seq(0, 1, 0.1),
	cell.text = rev(rowSums(outlier.num.subtype.5.patient.table.meta.5.only)),
	text.cex = text.size,
	text.fontface = 1,
	xaxis.fontface = 1,
	xaxis.lab = c('    ', 'Total', rep('', 4)),
	xaxis.rot = 0, 
	yaxis.lab = rep('', 8),
	xaxis.cex = 1,
	grid.row = TRUE,
	grid.col = FALSE,
	col.pos = rep(2, nrow(outlier.num.subtype.5.patient.table.meta.5.only)),
	row.pos = 1:nrow(outlier.num.subtype.5.patient.table.meta.5.only),
	print.colour.key = FALSE,
	same.as.matrix = TRUE,
	xaxis.tck = c(0, 0),
	yaxis.tck = c(0, 0),
	use.legacy.settings = FALSE,
	);

# column total heatmap
col.total <- create.heatmap(
	x = data.frame(colSums(outlier.num.subtype.5.patient.table.meta.5.only)) / sum(outlier.num.subtype.5.patient.table.meta.5.only),
	clustering = 'none',
	colour.scheme = c('white', 'dodgerblue'),
	colour.alpha = 1,
	at = seq(0, 1, 0.1),
	cell.text = colSums(outlier.num.subtype.5.patient.table.meta.5.only),
	text.cex = text.size,
	text.fontface = 1,
	xaxis.fontface = 1,
	xaxis.lab = rep('', 8),
	yaxis.lab = expression('Total'),
	yaxis.cex = 1.1,
	yat = 1.5,
	grid.row = FALSE,
	grid.col = TRUE,
	col.pos = 1:ncol(outlier.num.subtype.5.patient.table.meta.5.only),
	row.pos = rep(1.5, ncol(outlier.num.subtype.5.patient.table.meta.5.only)),
	print.colour.key = FALSE,
	yaxis.tck = c(0, 0),
	xaxis.tck = c(0, 0),
	xaxis.top.tck = c(0, 0),
	same.as.matrix = TRUE,
	use.legacy.settings = FALSE,
	);


# Include fisher and odd ratio
p.value.outlier.subtype.5.fisher.meta <- NULL;
p.value.subtype.5.odd.sub.meta <- NULL;
for (i in 1:5) {
    total.subtype <- nrow(subtype.5.total.outlier.num.meta); # total number of patients
    target.subtype <- subtype.5.meta.status$Freq[i] # number of the subtypes of interest in the population
    total.outlier <- sum(subtype.5.total.outlier.num.1.meta$outlier > 0) # total sample size
    target.outlier <- outlier.subtype.5.meta.status$Freq[outlier.subtype.5.meta.status$outlier == 1][i] # number of patients with outliers 
    
    # Perform the hypergeometric test
    p_value <- fisher.test(matrix(c(target.outlier, total.outlier - target.outlier, target.subtype - target.outlier, total.subtype - total.outlier - target.subtype + target.outlier), nrow=2), alternative="two.sided")$p.value;
    p.value.outlier.subtype.5.fisher.meta <- c(p.value.outlier.subtype.5.fisher.meta, p_value);
    
    odd.ratio <- fisher.test(matrix(c(target.outlier, total.outlier - target.outlier, target.subtype - target.outlier, total.subtype - total.outlier - target.subtype + target.outlier), nrow=2), alternative="two.sided")$estimate
    p.value.subtype.5.odd.sub.meta <- c(p.value.subtype.5.odd.sub.meta, odd.ratio);
    }

names(p.value.subtype.5.odd.sub.meta) <- subtype.5.meta.status$Var1;
p.value.subtype.5.odd.sub.meta <- p.value.subtype.5.odd.sub.meta[c('Basal', 'Her2', 'LumB', 'Normal', 'LumA')]

# p.value.outlier.subtype.5.fisher.fdr.meta <- p.adjust(p.value.outlier.subtype.5.fisher.meta, method = 'BH');

# column for the odd ratio - color on FDR
row.total.4 <- create.heatmap(
	x = data.frame(cbind(log2(c(p.value.subtype.5.odd.sub.meta)),
	                     log2(c(p.value.subtype.5.odd.sub.meta)),
	                     log2(c(p.value.subtype.5.odd.sub.meta)))),
	clustering = 'none',
	colour.scheme = c('#107090','white', '#b2402b'),
	colour.alpha = 1,
	at = seq(-2, 2, 0.1),
	cell.text = c(rev(round(c(p.value.subtype.5.odd.sub.meta), digits = 2))),
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
	col.pos = rep(2, nrow(outlier.num.subtype.5.patient.table.meta.5.only)),
	row.pos = 1:nrow(outlier.num.subtype.5.patient.table.meta.5.only),
	print.colour.key = FALSE,
	same.as.matrix = TRUE,
	xaxis.tck = c(0, 0),
	yaxis.tck = c(0, 0),
	use.legacy.settings = FALSE,
	);


# Only with odd and fdr
text.pvalue.subtype <- display.statistical.result(
    x = chisq.test(outlier.num.subtype.5.patient.table.meta.5.only)$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

key.subtype <- list(
    text = list(
        lab = text.pvalue.subtype, 
        cex = 1
        ),
    x = 0.25,
    y = 0.95
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
    plot.objects.heights = c(10, 1.5),
	plot.objects.widths = c(10, 1.2, 1.2),
	x.spacing = -2.2, 
	y.spacing = -6,     # Legend
    legend = list(
        bottom = list(
            fun = draw.key,
            args = list(
                key = list(
                    text = list(
                        lab = text.pvalue.subtype
                        ),
                    cex = 1.1,
                just = 'right'
                    )
                )
            )
        ),
	bottom.padding = 0,
    bottom.legend.padding = -0.5,
	top.padding = 2,
	right.padding = 0
    );


save.outlier.figure(
    new.plot,
    c('Figure3a', 'meta_subtype', 'box'),
    width = 7,
    height = 6
    );





# 3. ICGC BRCA-EU
icgc.clinic.subtype.order <- icgc.clinic.order[match(colnames(outlier.patient.tag.01.icgc), icgc.clinic.order$sample), ];

icgc.clinic.subtype.order.data <- data.frame(as.character(icgc.clinic.subtype.order$subtype));
icgc.clinic.subtype.order.data[is.na(icgc.clinic.subtype.order.data$as.character.icgc.clinic.subtype.order.subtype.)| 
        icgc.clinic.subtype.order.data$as.character.icgc.clinic.subtype.order.subtype. == "", ] <- 6;
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

matador.clinic.order.data <- data.frame(matador.clinic.order$subtype);
matador.clinic.order.data[is.na(matador.clinic.order.data$matador.clinic.order.PAM50.Subtype), ] <- 6;
matador.clinic.order.data[matador.clinic.order.data$matador.clinic.order.PAM50.Subtype == 'Basal', ] <- 1;
matador.clinic.order.data[matador.clinic.order.data$matador.clinic.order.PAM50.Subtype == 'Her2', ] <- 2;
matador.clinic.order.data[matador.clinic.order.data$matador.clinic.order.PAM50.Subtype == 'LumA', ] <- 3;
matador.clinic.order.data[matador.clinic.order.data$matador.clinic.order.PAM50.Subtype == 'LumB', ] <- 4;
matador.clinic.order.data[matador.clinic.order.data$matador.clinic.order.PAM50.Subtype == 'Normal', ] <- 5;
matador.clinic.order.data.num <- data.frame(as.numeric(matador.clinic.order.data$matador.clinic.order.subtype));
rownames(matador.clinic.order.data.num) <- colnames(outlier.patient.tag.01.matador);

outlier.patient.tag.01.matador.sum <- apply(outlier.patient.tag.01.matador, 2, sum);
subtype.total.outlier.num.matador <- data.frame(cbind(
    subtype = matador.clinic.order.data.num,
    outlier = outlier.patient.tag.01.matador.sum
    ));
colnames(subtype.total.outlier.num.matador) <- c('subtype', 'outlier');

subtype.total.outlier.num.1.matador <- subtype.total.outlier.num.matador;
subtype.total.outlier.num.1.matador$outlier[subtype.total.outlier.num.1.matador$outlier > 0] <- 1;
outlier.subtype.matador.status <- data.frame(table(subtype.total.outlier.num.1.matador));
subtype.matador.status <- data.frame(table(subtype.total.outlier.num.matador$subtype));


# 6. Cheng
patient.cheng$subtype.genefu <- factor(
  patient.cheng$subtype.genefu,
  levels = c("Basal", "Her2", "LumA", "LumB", "Normal")
    )
patient.cheng <- data.frame(patient.cheng);
rownames(patient.cheng) <-  colnames(outlier.patient.tag.01.cheng);

patient.cheng.data <- data.frame(as.character(patient.cheng$subtype.genefu));
patient.cheng.data[is.na(patient.cheng.data$as.character.patient.cheng.subtype.genefu.),] <- 6;
patient.cheng.data[patient.cheng.data$as.character.patient.cheng.subtype.genefu. == 'Basal',] <- 1;
patient.cheng.data[patient.cheng.data$as.character.patient.cheng.subtype.genefu. == 'Her2',] <- 2;
patient.cheng.data[patient.cheng.data$as.character.patient.cheng.subtype.genefu. == 'LumA',] <- 3;
patient.cheng.data[patient.cheng.data$as.character.patient.cheng.subtype.genefu. == 'LumB',] <- 4;
patient.cheng.data[patient.cheng.data$as.character.patient.cheng.subtype.genefu. == 'Normal',] <- 5;
patient.cheng.data.num <- data.frame(as.numeric(patient.cheng.data$as.character.patient.cheng.subtype.genefu.));


outlier.patient.tag.01.cheng.sum <- apply(outlier.patient.tag.01.cheng, 2, sum);
subtype.total.outlier.num.cheng <- data.frame(cbind(subtype = patient.cheng.data.num,
                                   outlier = outlier.patient.tag.01.cheng.sum));
colnames(subtype.total.outlier.num.cheng) <- c("subtype", "outlier");
subtype.total.outlier.num.1.cheng <- subtype.total.outlier.num.cheng;
subtype.total.outlier.num.1.cheng$outlier[subtype.total.outlier.num.1.cheng$outlier > 0] <- 1;
outlier.subtype.cheng.status <- data.frame(table(subtype.total.outlier.num.1.cheng));
subtype.cheng.status <- data.frame(table(subtype.total.outlier.num.cheng$subtype));


# 7. Hatzis
patient.hatzis$subtype <- factor(
  patient.hatzis$subtype,
  levels = c("Basal", "Her2", "LumA", "LumB", "Normal")
    )
patient.hatzis <- data.frame(patient.hatzis);
rownames(patient.hatzis) <-  colnames(outlier.patient.tag.01.hatzis);

patient.hatzis.data <- data.frame(as.character(patient.hatzis$subtype));
patient.hatzis.data[is.na(patient.hatzis.data$as.character.patient.hatzis.subtype.),] <- 6;
patient.hatzis.data[patient.hatzis.data$as.character.patient.hatzis.subtype. == 'Basal',] <- 1;
patient.hatzis.data[patient.hatzis.data$as.character.patient.hatzis.subtype. == 'Her2',] <- 2;
patient.hatzis.data[patient.hatzis.data$as.character.patient.hatzis.subtype. == 'LumA',] <- 3;
patient.hatzis.data[patient.hatzis.data$as.character.patient.hatzis.subtype. == 'LumB',] <- 4;
patient.hatzis.data[patient.hatzis.data$as.character.patient.hatzis.subtype. == 'Normal',] <- 5;
patient.hatzis.data.num <- data.frame(as.numeric(patient.hatzis.data$as.character.patient.hatzis.subtype.));


outlier.patient.tag.01.hatzis.sum <- apply(outlier.patient.tag.01.hatzis, 2, sum);
subtype.total.outlier.num.hatzis <- data.frame(cbind(subtype = patient.hatzis.data.num,
                                   outlier = outlier.patient.tag.01.hatzis.sum));
colnames(subtype.total.outlier.num.hatzis) <- c("subtype", "outlier");
subtype.total.outlier.num.1.hatzis <- subtype.total.outlier.num.hatzis;
subtype.total.outlier.num.1.hatzis$outlier[subtype.total.outlier.num.1.hatzis$outlier > 0] <- 1;
outlier.subtype.hatzis.status <- data.frame(table(subtype.total.outlier.num.1.hatzis));
subtype.hatzis.status <- data.frame(table(subtype.total.outlier.num.hatzis$subtype));


# 8. Sjostrom
patient.sjostrom$subtype.genefu <- factor(
  patient.sjostrom$subtype.genefu,
  levels = c("Basal", "Her2", "LumA", "LumB", "Normal")
    )
patient.sjostrom <- data.frame(patient.sjostrom);
rownames(patient.sjostrom) <-  colnames(outlier.patient.tag.01.sjostrom);

patient.sjostrom.data <- data.frame(as.character(patient.sjostrom$subtype.genefu));
patient.sjostrom.data[is.na(patient.sjostrom.data$as.character.patient.sjostrom.subtype.genefu.),] <- 6;
patient.sjostrom.data[patient.sjostrom.data$as.character.patient.sjostrom.subtype.genefu. == 'Basal',] <- 1;
patient.sjostrom.data[patient.sjostrom.data$as.character.patient.sjostrom.subtype.genefu. == 'Her2',] <- 2;
patient.sjostrom.data[patient.sjostrom.data$as.character.patient.sjostrom.subtype.genefu. == 'LumA',] <- 3;
patient.sjostrom.data[patient.sjostrom.data$as.character.patient.sjostrom.subtype.genefu. == 'LumB',] <- 4;
patient.sjostrom.data[patient.sjostrom.data$as.character.patient.sjostrom.subtype.genefu. == 'Normal',] <- 5;
patient.sjostrom.data.num <- data.frame(as.numeric(patient.sjostrom.data$as.character.patient.sjostrom.subtype.genefu.));


outlier.patient.tag.01.sjostrom.sum <- apply(outlier.patient.tag.01.sjostrom, 2, sum);
subtype.total.outlier.num.sjostrom <- data.frame(cbind(subtype = patient.sjostrom.data.num,
                                   outlier = outlier.patient.tag.01.sjostrom.sum));
colnames(subtype.total.outlier.num.sjostrom) <- c("subtype", "outlier");
subtype.total.outlier.num.1.sjostrom <- subtype.total.outlier.num.sjostrom;
subtype.total.outlier.num.1.sjostrom$outlier[subtype.total.outlier.num.1.sjostrom$outlier > 0] <- 1;
outlier.subtype.sjostrom.status <- data.frame(table(subtype.total.outlier.num.1.sjostrom));
subtype.sjostrom.status <- data.frame(table(subtype.total.outlier.num.sjostrom$subtype));


# 9. Kao
patient.kao$subtype.genefu <- factor(
  patient.kao$subtype.genefu,
  levels = c("Basal", "Her2", "LumA", "LumB", "Normal")
    )
patient.kao <- data.frame(patient.kao);
rownames(patient.kao) <-  colnames(outlier.patient.tag.01.kao);

patient.kao.data <- data.frame(as.character(patient.kao$subtype.genefu));
patient.kao.data[is.na(patient.kao.data$as.character.patient.kao.subtype.genefu.),] <- 6;
patient.kao.data[patient.kao.data$as.character.patient.kao.subtype.genefu. == 'Basal',] <- 1;
patient.kao.data[patient.kao.data$as.character.patient.kao.subtype.genefu. == 'Her2',] <- 2;
patient.kao.data[patient.kao.data$as.character.patient.kao.subtype.genefu. == 'LumA',] <- 3;
patient.kao.data[patient.kao.data$as.character.patient.kao.subtype.genefu. == 'LumB',] <- 4;
patient.kao.data[patient.kao.data$as.character.patient.kao.subtype.genefu. == 'Normal',] <- 5;
patient.kao.data.num <- data.frame(as.numeric(patient.kao.data$as.character.patient.kao.subtype.genefu.));


outlier.patient.tag.01.kao.sum <- apply(outlier.patient.tag.01.kao, 2, sum);
subtype.total.outlier.num.kao <- data.frame(cbind(subtype = patient.kao.data.num,
                                   outlier = outlier.patient.tag.01.kao.sum));
colnames(subtype.total.outlier.num.kao) <- c("subtype", "outlier");
subtype.total.outlier.num.1.kao <- subtype.total.outlier.num.kao;
subtype.total.outlier.num.1.kao$outlier[subtype.total.outlier.num.1.kao$outlier > 0] <- 1;
outlier.subtype.kao.status <- data.frame(table(subtype.total.outlier.num.1.kao));
subtype.kao.status <- data.frame(table(subtype.total.outlier.num.kao$subtype));




# Function to perform Fisher's exact test and odds ratio calculation
perform.fisher.test <- function(total.patients, total.outliers, subtype.freq, outlier.freq) {
    fisher.test.result <- fisher.test(
        matrix(
            c(
                # outlier.freq + 1,
                # total.outliers - outlier.freq + 1,
                # subtype.freq - outlier.freq + 1,
                # total.patients - total.outliers - subtype.freq + outlier.freq + 1
                outlier.freq ,
                total.outliers - outlier.freq ,
                subtype.freq - outlier.freq ,
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

# Hatzis analysis
hatzis.results <- perform.subtype.analysis(subtype.total.outlier.num.hatzis, subtype.hatzis.status$Freq, outlier.subtype.hatzis.status[outlier.subtype.hatzis.status$outlier == 1, ]$Freq);

# cheng analysis
cheng.results <- perform.subtype.analysis(subtype.total.outlier.num.cheng, subtype.cheng.status$Freq, outlier.subtype.cheng.status[outlier.subtype.cheng.status$outlier == 1, ]$Freq);

# kao analysis
kao.results <- perform.subtype.analysis(subtype.total.outlier.num.kao, subtype.kao.status$Freq, outlier.subtype.kao.status[outlier.subtype.kao.status$outlier == 1, ]$Freq);

# sjostrom analysis
sjostrom.results <- perform.subtype.analysis(subtype.total.outlier.num.sjostrom, subtype.sjostrom.status$Freq, outlier.subtype.sjostrom.status[outlier.subtype.sjostrom.status$outlier == 1, ]$Freq);




all.odd.subtype <- cbind(
    meta.results$odd.ratio,
    brca.results$odd.ratio,
    ispy.results$odd.ratio,
    sjostrom.results$odd.ratio,
    cheng.results$odd.ratio,
    matador.results$odd.ratio,
    icgc.results$odd.ratio,
    kao.results$odd.ratio,
    hatzis.results$odd.ratio
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
    at = seq(log(0.2), log(5), 0.1),
    cell.text = round(data.frame(all.odd.subtype.table)$Freq, digits = 1),
    text.cex = 1.2,
    text.fontface = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    grid.row = TRUE,
    grid.col = TRUE,
    ylab.label = expression('Subtype'),
    yaxis.lab = rownames(all.odd.subtype.table),
    xaxis.lab.top = c(
        'METABRIC',
        'TCGA-BRCA',
        'I-SPY2',
        'Sjostrom', 
        'Cheng', 
        'matador',
        'ICGC BRCA-EU',
        'Kao', 
        'Hatzis'),
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
    log(meta.results$odd.ratio),
    log(brca.results$odd.ratio),
    log(ispy.results$odd.ratio),
    log(sjostrom.results$odd.ratio),
    log(cheng.results$odd.ratio),
    log(matador.results$odd.ratio),
    log(icgc.results$odd.ratio),
    log(kao.results$odd.ratio),
    log(hatzis.results$odd.ratio)
    );

se.odd.list <- list(
    (log(meta.results$ci[, 2]) - log(meta.results$ci[, 1])) / 3.92,
    (log(brca.results$ci[, 2]) - log(brca.results$ci[, 1])) / 3.92,
    (log(ispy.results$ci[, 2]) - log(ispy.results$ci[, 1])) / 3.92,
    (log(sjostrom.results$ci[, 2]) - log(sjostrom.results$ci[, 1])) / 3.92,
    (log(cheng.results$ci[, 2]) - log(cheng.results$ci[, 1])) / 3.92,
    (log(matador.results$ci[, 2]) - log(matador.results$ci[, 1])) / 3.92,
    (log(icgc.results$ci[, 2]) - log(icgc.results$ci[, 1])) / 3.92,
    (log(kao.results$ci[, 2]) - log(kao.results$ci[, 1])) / 3.92,
    (log(hatzis.results$ci[, 2]) - log(hatzis.results$ci[, 1])) / 3.92
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
    plot.objects.widths = c(3, 1.2),
    ylab.axis.padding = -10,
    x.spacing = 3,
    right.legend.padding = 0
    );

multi.gene;

save.outlier.figure(
    multi.gene,
    c('Figure3b', 'subtype', 'multi'),
    width = 7,
    height = 5
    );

### SAVE VARIABLES #############################################################
# Cache the important variables for later use
cache.multiple.computed.variables(c(
    'subtype.total.outlier.num.1.brca',
    'subtype.5.total.outlier.num.meta',
    'subtype.total.outlier.num.1.ispy',
    'subtype.total.outlier.num.1.cheng',
    'subtype.total.outlier.num.1.sjostrom',
    'subtype.total.outlier.num.1.matador',
    'subtype.total.outlier.num.1.icgc',
    'subtype.total.outlier.num.1.kao',
    'subtype.total.outlier.num.1.hatzis'
    ))

save.session.profile(file.path('output', 'Figure3abf.txt'));
