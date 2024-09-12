### HISTORY #####################################################################
# This script performs meta-analysis to compare the tumour stages across 
# multiple breast cancer subtypes. 
# Date: 2024-08-15


# Load required libraries
library(metafor);

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));

### 1. TCGA-BRCA
convert.stage <- function(stage) {
    result <- rep(NA, length(stage));
    result[stage %in% c("STAGE I", "STAGE IA", "STAGE IB")] <- 1;
    result[stage %in% c("STAGE II", "STAGE IIA", "STAGE IIB")] <- 2;
    result[stage %in% c("STAGE III", "STAGE IIIA", "STAGE IIIB", "STAGE IIIC")] <- 3;
    result[stage == "STAGE IV"] <- 4;
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
        
        fisher.result <- fisher.test(matrix(c(table.1+1, table.2+1, table.3+1, table.4+1), nrow=2), alternative="two.sided");
        
        p.values <- c(p.values, fisher.result$p.value);
        odd.ratios <- c(odd.ratios, fisher.result$estimate);
        ci.intervals <- rbind(ci.intervals, fisher.result$conf.int);
        }
    
    p.values.fdr <- p.adjust(p.values, method = 'BH');
    
    return(list(p.values = p.values, odd.ratios = odd.ratios, ci.intervals = ci.intervals, p.values.fdr = p.values.fdr));
    }



# brca.clinic <- read.csv(
#     file = '/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/data/brca_tcga_pan_can_atlas_2018_clinical_data.csv', 
#     header = TRUE, 
#     stringsAsFactors = F, 
#     sep = ',');


brca.sample.clinic <- gsub("-", ".", brca.clinic$Patient.ID, fixed = TRUE);
rownames(brca.clinic) <- brca.sample.clinic;
brca.clinic.sort <- brca.clinic[order(brca.clinic$Subtype),];
brca.clinic.order <- brca.clinic.sort[substr(colnames(brca.outlier.patient.tag.01.t.p.order), 1, 12),]

brca.outlier.patient.tag.01.t.p.order.sum <- apply(brca.outlier.patient.tag.01.t.p.order, 2, sum);

os.data.stage.pseudo.brca <- data.frame(
    status = substr(brca.clinic.order$Overall.Survival.Status, 1, 1),
    os = brca.clinic.order$Overall.Survival..Months.,
    patient = brca.outlier.patient.tag.01.t.p.order.sum,
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
p.value.outlier.stage.pseudo.fisher.t1.brca <- all.results$p.values;
p.value.stage.pseudo.odd.sub.t1.brca <- all.results$odd.ratios;
p.value.stage.pseudo.ci.sub.t1.brca <- all.results$ci.intervals;
p.value.outlier.stage.pseudo.fisher.fdr.t1.brca <- all.results$p.values.fdr;

# Basal
basal.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, "Basal");
p.value.outlier.stage.pseudo.basal.fisher.t1.brca <- basal.results$p.values;
p.value.stage.pseudo.basal.odd.sub.t1.brca <- basal.results$odd.ratios;
p.value.stage.pseudo.basal.ci.sub.t1.brca <- basal.results$ci.intervals;
p.value.outlier.stage.pseudo.basal.fisher.fdr.t1.brca <- basal.results$p.values.fdr;

# Her2
her2.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, "Her2");
p.value.outlier.stage.pseudo.her2.fisher.t1.brca <- her2.results$p.values;
p.value.stage.pseudo.her2.odd.sub.t1.brca <- her2.results$odd.ratios;
p.value.stage.pseudo.her2.ci.sub.t1.brca <- her2.results$ci.intervals;
p.value.outlier.stage.pseudo.her2.fisher.fdr.t1.brca <- her2.results$p.values.fdr;

# LumA
luma.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, "LumA");
p.value.outlier.stage.pseudo.luma.fisher.t1.brca <- luma.results$p.values;
p.value.stage.pseudo.luma.odd.sub.t1.brca <- luma.results$odd.ratios;
p.value.stage.pseudo.luma.ci.sub.t1.brca <- luma.results$ci.intervals;
p.value.outlier.stage.pseudo.luma.fisher.fdr.t1.brca <- luma.results$p.values.fdr;

# LumB
lumb.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, "LumB");
p.value.outlier.stage.pseudo.lumb.fisher.t1.brca <- lumb.results$p.values;
p.value.stage.pseudo.lumb.odd.sub.t1.brca <- lumb.results$odd.ratios;
p.value.stage.pseudo.lumb.ci.sub.t1.brca <- lumb.results$ci.intervals;
p.value.outlier.stage.pseudo.lumb.fisher.fdr.t1.brca <- lumb.results$p.values.fdr;

# Normal
normal.results <- perform.fisher.test.brca(os.data.stage.pseudo.brca, "Normal");
p.value.outlier.stage.pseudo.normal.fisher.t1.brca <- normal.results$p.values;
p.value.stage.pseudo.normal.odd.sub.t1.brca <- normal.results$odd.ratios;
p.value.stage.pseudo.normal.ci.sub.t1.brca <- normal.results$ci.intervals;
p.value.outlier.stage.pseudo.normal.fisher.fdr.t1.brca <- normal.results$p.values.fdr;




### 2. METABRIC
outlier.patient.tag.01.meta.sum <- apply(outlier.patient.tag.01.meta, 2, sum)

# meta.clinic.5 <- read.delim2(file = '/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/data/patient_combine.txt', row.names = 1, header = T);
# meta.clinic.5.order <- meta.clinic.5[names(outlier.patient.tag.01.meta.sum),];
# 
# # get subtype info from other clinical dataset
# meta.clinic <- read.delim2(file = '/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/data/data_clinical_patient.txt', row.names = 1, header = T);
# meta.clinic <- meta.clinic[5:nrow(meta.clinic),];
# rownames(meta.clinic) <- gsub("-", ".", rownames(meta.clinic));
# meta.clinic.order <- meta.clinic[names(outlier.patient.tag.01.meta.sum),];

# meta.clinic.5.order.combine <- meta.clinic.order;
meta.clinic.5.order.combine$pam50 <- meta.clinic.5.order$Pam50Subtype;


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
        
        fisher.result <- fisher.test(matrix(c(table.1+1, table.2+1, table.3+1, table.4+1), nrow=2), alternative="two.sided")
        
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
    stage = meta.clinic.5.order$stage
    );

os.data.stage.pseudo.meta[,1] <- as.numeric(os.data.stage.pseudo.meta[,1]);
os.data.stage.pseudo.meta[,2] <- as.numeric(os.data.stage.pseudo.meta[,2]);
os.data.stage.pseudo.meta[,3] <- as.numeric(os.data.stage.pseudo.meta[,3]);
os.data.stage.pseudo.meta[,5] <- as.numeric(os.data.stage.pseudo.meta[,5]);

os.data.stage.pseudo.meta.stage <- data.frame(table(os.data.stage.pseudo.meta$stage));
os.data.stage.pseudo.meta.stage <- os.data.stage.pseudo.meta.stage[2:5,];


# Metabric - All subtypes
meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta);
p.value.outlier.stage.pseudo.fisher.t1.meta <- meta.results$p.values;
p.value.stage.pseudo.odd.sub.t1.meta <- meta.results$odd.ratios;
p.value.stage.pseudo.ci.sub.t1.meta <- meta.results$ci.intervals;
p.value.outlier.stage.pseudo.fisher.fdr.t1.meta <- meta.results$p.values.fdr;

# Basal subtype
basal.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, "Basal");
p.value.outlier.stage.pseudo.basal.fisher.t1.meta <- basal.meta.results$p.values;
p.value.stage.pseudo.basal.odd.sub.t1.meta <- basal.meta.results$odd.ratios;
p.value.stage.pseudo.basal.ci.sub.t1.meta <- basal.meta.results$ci.intervals;
p.value.outlier.stage.pseudo.basal.fisher.fdr.t1.meta <- basal.meta.results$p.values.fdr;

# Her2 subtype
her2.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, "Her2");
p.value.outlier.stage.pseudo.her2.fisher.t1.meta <- her2.meta.results$p.values;
p.value.stage.pseudo.her2.odd.sub.t1.meta <- her2.meta.results$odd.ratios;
p.value.stage.pseudo.her2.ci.sub.t1.meta <- her2.meta.results$ci.intervals;
p.value.outlier.stage.pseudo.her2.fisher.fdr.t1.meta <- her2.meta.results$p.values.fdr;   

# Luminal A subtype
luma.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, "LumA");
p.value.outlier.stage.pseudo.luma.fisher.t1.meta <- luma.meta.results$p.values;
p.value.stage.pseudo.luma.odd.sub.t1.meta <- luma.meta.results$odd.ratios;
p.value.stage.pseudo.luma.ci.sub.t1.meta <- luma.meta.results$ci.intervals;
p.value.outlier.stage.pseudo.luma.fisher.fdr.t1.meta <- luma.meta.results$p.values.fdr;   

# Luminal B subtype
lumb.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, "LumB");
p.value.outlier.stage.pseudo.lumb.fisher.t1.meta <- lumb.meta.results$p.values;
p.value.stage.pseudo.lumb.odd.sub.t1.meta <- lumb.meta.results$odd.ratios;
p.value.stage.pseudo.lumb.ci.sub.t1.meta <- lumb.meta.results$ci.intervals;
p.value.outlier.stage.pseudo.lumb.fisher.fdr.t1.meta <- lumb.meta.results$p.values.fdr;   

# Normal subtype
normal.meta.results <- perform.fisher.test.meta(os.data.stage.pseudo.meta, "Normal");
p.value.outlier.stage.pseudo.normal.fisher.t1.meta <- normal.meta.results$p.values;
p.value.stage.pseudo.normal.odd.sub.t1.meta <- normal.meta.results$odd.ratios;
p.value.stage.pseudo.normal.ci.sub.t1.meta <- normal.meta.results$ci.intervals;
p.value.outlier.stage.pseudo.normal.fisher.fdr.t1.meta <- normal.meta.results$p.values.fdr;





# ### 3. ICGC BRCA-EU
# icgc.clinic <- read.delim2(file = '/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/data/Supplementary.Table.1 CLINICAL.PATHOLOGY.DATA.FREEZE.ANALYSIS.v4.032015.csv', header = T, sep =',');
# icgc.clinic.order <- icgc.clinic[match(colnames(outlier.patient.tag.01.icgc), icgc.clinic$sample),];
# colnames(icgc.clinic) <- icgc.clinic[1,];
# icgc.clinic <- icgc.clinic[-1,];
icgc.sample.num <- substr(icgc.clinic$sample_name, 3, nchar(icgc.clinic$sample_name));

# Prepare ICGC sample numbers
icgc.sample.num.nra <- gsub("PR(\\d+)a(\\.RNA|\\.2)?", "\\1", colnames(outlier.patient.tag.01.icgc));
icgc.sample.num.nra <- gsub("PR(\\d+)b(\\.RNA|\\.2)?", "\\1", icgc.sample.num.nra);
icgc.sample.num.nra <- gsub("PR(\\d+)c(\\.RNA|\\.2)?", "\\1", icgc.sample.num.nra);

icgc.clinic.all.order <- icgc.clinic[match(icgc.sample.num.nra, icgc.sample.num), ];
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
        } 
    else if ((T == 'T1' && N == 'N0' && M == 'M0') || (T %in% c('T0', 'T1') && N == 'N1mi' && M == 'M0')) {
        stage <- 1;
        } 
    else if ((T == 'T0' && N == 'N1' && M == 'M0') || (T == 'T1' && N == 'N1' && M == 'M0') || (T == 'T2' && N == 'N0' && M == 'M0')) {
        stage <- 2;
        } 
    else if ((T == 'T2' && N == 'N1' && M == 'M0') || (T == 'T3' && N == 'N0' && M == 'M0')) {
        stage <- 2;
        } 
    else if ((T %in% c('T0', 'T1', 'T2', 'T3') && N == 'N2' && M == 'M0') || (T == 'T3' && N == 'N1' && M == 'M0')) {
        stage <- 3;
        } 
    else if (T == 'T4' && (N %in% c('N0', 'N1', 'N2')) && M == 'M0') {
        stage <- 3
        } 
    else if ((T %in% c('T0', 'T1', 'T2', 'T3', 'T4') && N == 'N3' && M == 'M0')) {
        stage <- 3
        } 
    else if (M %in% c('M1', 'M2')) {
        stage <- 4
        } 
    else {
        stage <- NA  # Set stage as NA if none of the conditions match
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
        
        fisher.result <- fisher.test(matrix(c(table.1 + 1, table.2 + 1, table.3 + 1, table.4 + 1), nrow = 2), alternative = "two.sided");
        
        p.values <- c(p.values, fisher.result$p.value);
        odd.ratios <- c(odd.ratios, fisher.result$estimate);
        ci.intervals <- rbind(ci.intervals, fisher.result$conf.int);
        }
    
    p.values.fdr <- p.adjust(p.values, method = 'BH');
    
    return(list(p.values = p.values, odd.ratios = odd.ratios, ci.intervals = ci.intervals, p.values.fdr = p.values.fdr));
    }


# Data preparation
os.data.stage.pseudo.icgc <- data.frame(cbind(subtype.total.outlier.num.icgc,
                                        stage = icgc.stage.convert.new.x));
os.data.stage.pseudo.icgc <- na.omit(os.data.stage.pseudo.icgc);

# All ICGC
all.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc);
p.value.outlier.stage.pseudo.fisher.t1.icgc <- all.results.icgc$p.values;
p.value.stage.pseudo.odd.sub.t1.icgc <- all.results.icgc$odd.ratios;
p.value.stage.pseudo.ci.sub.t1.icgc <- all.results.icgc$ci.intervals;
p.value.outlier.stage.pseudo.fisher.fdr.t1.icgc <- all.results.icgc$p.values.fdr;

# Basal
basal.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 1);
p.value.outlier.stage.pseudo.basal.fisher.t1.icgc <- basal.results.icgc$p.values;
p.value.stage.pseudo.basal.odd.sub.t1.icgc <- basal.results.icgc$odd.ratios;
p.value.stage.pseudo.basal.ci.sub.t1.icgc <- basal.results.icgc$ci.intervals;
p.value.outlier.stage.pseudo.basal.fisher.fdr.t1.icgc <- basal.results.icgc$p.values.fdr;

# Her2
her2.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 2);
p.value.outlier.stage.pseudo.her2.fisher.t1.icgc <- her2.results.icgc$p.values;
p.value.stage.pseudo.her2.odd.sub.t1.icgc <- her2.results.icgc$odd.ratios;
p.value.stage.pseudo.her2.ci.sub.t1.icgc <- her2.results.icgc$ci.intervals;
p.value.outlier.stage.pseudo.her2.fisher.fdr.t1.icgc <- her2.results.icgc$p.values.fdr;

# Lum A
luma.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 3);
p.value.outlier.stage.pseudo.luma.fisher.t1.icgc <- luma.results.icgc$p.values;
p.value.stage.pseudo.luma.odd.sub.t1.icgc <- luma.results.icgc$odd.ratios;
p.value.stage.pseudo.luma.ci.sub.t1.icgc <- luma.results.icgc$ci.intervals;
p.value.outlier.stage.pseudo.luma.fisher.fdr.t1.icgc <- luma.results.icgc$p.values.fdr;

# Lum B
lumb.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 4);
p.value.outlier.stage.pseudo.lumb.fisher.t1.icgc <- lumb.results.icgc$p.values;
p.value.stage.pseudo.lumb.odd.sub.t1.icgc <- lumb.results.icgc$odd.ratios;
p.value.stage.pseudo.lumb.ci.sub.t1.icgc <- lumb.results.icgc$ci.intervals;
p.value.outlier.stage.pseudo.lumb.fisher.fdr.t1.icgc <- lumb.results.icgc$p.values.fdr;

# Normal
normal.results.icgc <- perform.fisher.test.icgc(os.data.stage.pseudo.icgc, subtype = 5);
p.value.outlier.stage.pseudo.normal.fisher.t1.icgc <- normal.results.icgc$p.values;
p.value.stage.pseudo.normal.odd.sub.t1.icgc <- normal.results.icgc$odd.ratios;
p.value.stage.pseudo.normal.ci.sub.t1.icgc <- normal.results.icgc$ci.intervals;
p.value.outlier.stage.pseudo.normal.fisher.fdr.t1.icgc <- normal.results.icgc$p.values.fdr;





### Meta-analysis
# Calculate log odds and standard error
calculate.ln.odd.se <- function(odd, ci) {
    ln.odd <- log(odd);
    se.odd <- (log(ci[,2]) - log(ci[,1])) / 3.92;
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
            } 
        else {
            metafor.all <- c(NA, NA, NA, NA);
            }
        
        metafor.stage.pseudo.odd.ci.p.t1 <- rbind(metafor.stage.pseudo.odd.ci.p.t1, metafor.all);
        }
    
    metafor.stage.pseudo.odd.ci.p.data.t1 <- data.frame(
        p.value = metafor.stage.pseudo.odd.ci.p.t1[,4],
        odd = metafor.stage.pseudo.odd.ci.p.t1[,1],
        ci.min = metafor.stage.pseudo.odd.ci.p.t1[,2],
        ci.max = metafor.stage.pseudo.odd.ci.p.t1[,3]
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

basal.results <- perform.meta.analysis(
    list(ln.odd.brca.basal$ln.odd, ln.odd.meta.basal$ln.odd, ln.odd.icgc.basal$ln.odd),
    list(ln.odd.brca.basal$se.odd, ln.odd.meta.basal$se.odd, ln.odd.icgc.basal$se.odd)
    );

# 2. Her2 (excluding ICGC)
ln.odd.brca.her2 <- calculate.ln.odd.se(p.value.stage.pseudo.her2.odd.sub.t1.brca, p.value.stage.pseudo.her2.ci.sub.t1.brca);
ln.odd.meta.her2 <- calculate.ln.odd.se(p.value.stage.pseudo.her2.odd.sub.t1.meta, p.value.stage.pseudo.her2.ci.sub.t1.meta);

her2.results <- perform.meta.analysis(
    list(ln.odd.brca.her2$ln.odd, ln.odd.meta.her2$ln.odd),
    list(ln.odd.brca.her2$se.odd, ln.odd.meta.her2$se.odd)
    );

# 3. LumA
ln.odd.brca.luma <- calculate.ln.odd.se(p.value.stage.pseudo.luma.odd.sub.t1.brca, p.value.stage.pseudo.luma.ci.sub.t1.brca);
ln.odd.meta.luma <- calculate.ln.odd.se(p.value.stage.pseudo.luma.odd.sub.t1.meta, p.value.stage.pseudo.luma.ci.sub.t1.meta);
ln.odd.icgc.luma <- calculate.ln.odd.se(p.value.stage.pseudo.luma.odd.sub.t1.icgc, p.value.stage.pseudo.luma.ci.sub.t1.icgc);

luma.results <- perform.meta.analysis(
    list(ln.odd.brca.luma$ln.odd, ln.odd.meta.luma$ln.odd, ln.odd.icgc.luma$ln.odd),
    list(ln.odd.brca.luma$se.odd, ln.odd.meta.luma$se.odd, ln.odd.icgc.luma$se.odd)
    );

# 4. LumB
ln.odd.brca.lumb <- calculate.ln.odd.se(p.value.stage.pseudo.lumb.odd.sub.t1.brca, p.value.stage.pseudo.lumb.ci.sub.t1.brca);
ln.odd.meta.lumb <- calculate.ln.odd.se(p.value.stage.pseudo.lumb.odd.sub.t1.meta, p.value.stage.pseudo.lumb.ci.sub.t1.meta);
ln.odd.icgc.lumb <- calculate.ln.odd.se(p.value.stage.pseudo.lumb.odd.sub.t1.icgc, p.value.stage.pseudo.lumb.ci.sub.t1.icgc);

lumb.results <- perform.meta.analysis(
    list(ln.odd.brca.lumb$ln.odd, ln.odd.meta.lumb$ln.odd, ln.odd.icgc.lumb$ln.odd),
    list(ln.odd.brca.lumb$se.odd, ln.odd.meta.lumb$se.odd, ln.odd.icgc.lumb$se.odd)
    );

# 5. Normal (excluding ICGC)
ln.odd.brca.normal <- calculate.ln.odd.se(p.value.stage.pseudo.normal.odd.sub.t1.brca, p.value.stage.pseudo.normal.ci.sub.t1.brca);
ln.odd.meta.normal <- calculate.ln.odd.se(p.value.stage.pseudo.normal.odd.sub.t1.meta, p.value.stage.pseudo.normal.ci.sub.t1.meta);

normal.results <- perform.meta.analysis(
    list(ln.odd.brca.normal$ln.odd, ln.odd.meta.normal$ln.odd),
    list(ln.odd.brca.normal$se.odd, ln.odd.meta.normal$se.odd)
    );


# 6. All datasets
ln.odd.brca <- calculate.ln.odd.se(p.value.stage.pseudo.odd.sub.t1.brca, p.value.stage.pseudo.ci.sub.t1.brca);
ln.odd.meta <- calculate.ln.odd.se(p.value.stage.pseudo.odd.sub.t1.meta, p.value.stage.pseudo.ci.sub.t1.meta);
ln.odd.icgc <- calculate.ln.odd.se(p.value.stage.pseudo.odd.sub.t1.icgc, p.value.stage.pseudo.ci.sub.t1.icgc);

all.results <- perform.meta.analysis(
    list(ln.odd.brca$ln.odd, ln.odd.meta$ln.odd, ln.odd.icgc$ln.odd),
    list(ln.odd.brca$se.odd, ln.odd.meta$se.odd, ln.odd.icgc$se.odd)
    );
# Combine the results for heatmap generation
metafor.stage.pseudo.odd.ci.p.data.each.subtype.odd.t1.no3 <- data.frame(
    basal = basal.results$data$odd[1:2],
    her2 = her2.results$data$odd[1:2],
    luma = luma.results$data$odd[1:2],
    lumb = lumb.results$data$odd[1:2],
    normal = normal.results$data$odd[1:2]
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
            } 
        else if (x != background.cutoff) {
            return(as.expression(bquote('10'^-.(as.character(x)))));
            } 
        else {
            return(as.expression(bquote('<10'^-.(as.character(x)))));
            }
        }
    );

legend <- legend.grob(
    list(
        legend = list(
            title = expression(underline('FDR')),
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
spot.size.function <- function(x) { 0.1 + (1.5 * abs(x)); }
spot.colour.function <- function(x) {
    colours <- rep("white", length(x));
    colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1]; 
    colours[sign(x) == 1] <- default.colours(2, palette.type = "dotmap")[2]; 
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
        space = "right",
        points = list(
            cex = spot.size.function(seq(-2, 2, 1)),
            col = spot.colour.function(seq(-2, 2, 1)),
            pch = 19
            ),
        text = list(
            lab = c("0.25", "0.5", "1", "2", "4"),
            cex = 1,
            adj = 1,
            fontface = "bold"
            ),
        padding.text = 8
        ),
    key.top = 1,
    right.padding = 2,
    pch = 21,
    pch.border.col = "white",
    bg.data = t(-log10(metafor.stage.pseudo.odd.ci.p.data.each.subtype.fdr.t1.no3)),
    colourkey = FALSE,
    colourkey.cex = 0.95,
    bg.alpha = 1,
    colour.scheme = c("white", "black"),
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
    pch.border.col = "white",
    bg.data = t(-log10(p.adjust(metafor.stage.pseudo.odd.ci.p.data.t1[1:2,]$p.value, method = 'BH'))),
    colourkey = FALSE,
    bg.alpha = 1,
    colour.scheme = c("white", "black"),
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
    c('tumour', 'stage', 'multipanel'),
    width = 4.5,
    height = 5.6
    );
