### TCGA-BRCA Heatmapping #########################################################################

### Preamble ######################################################################################
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(BoutrosLab.statistics.survival);

#loading TCGA_BRCA
load(
    file = '2023-07-07_TCGA_BRCA_Outlier.rda',
    );
#loading Metabric
load(
    file = '2023-07-07_Metabric_Outlier.rda',
    );

### fpkm.normilization.zscore #####################################################################
#Description
    #normalizing the fpkm values for use in heatmaps
#Input Variables
    # (1) data
    # the name of the data frame that contains non-normalized fpkm values
#Output Variables
    # (1) norm.fpkm
    #data frame with normalized fpkm values
fpkm.normalization.zscore <- function(data) {
    number.rows <- nrow(data);
    number.cols <- ncol(data);
    fpkm.stats <- data.frame(
        Genes = row.names(data),
        Means.by.row = rep(NA,number.rows),
        SD.by.row = rep(NA,number.rows)
        );
    for (i in 1:number.rows) {
        fpkm.stats[i,2]<- mean(as.numeric(c(t(data[i,]))),na.rm = TRUE);
        fpkm.stats[i,3]<- sd(as.numeric(c(t(data[i,]))),na.rm = TRUE);
        }
    norm.fpkm <- data.frame(matrix(
        ncol = number.cols,nrow = number.rows,
        dimnames = list(row.names(data),colnames(data))));
    for (r in 1:number.rows) {
        for (c in 1:number.cols) {
            norm.fpkm[r,c] <- (as.numeric(data[r,c]) - fpkm.stats[r,2]) / fpkm.stats[r,3];
            }
        }
    result <- list(fpkm.stats,norm.fpkm)
    return(norm.fpkm)
    }

### fpkm.robust.zscore ############################################################################
#Description
    #normalizing the fpkm values with median and MAD
#Input Variables
    # (1) data:
    #the name of the data frame that contains non-normalized fpkm values
#Output Variables
    # (2) median.norm.fpkm
    #data frame with normalized fpkm values utilizing the median and IQR

fpkm.robust.zscore <- function(data) {
    number.rows <- nrow(data);
    number.cols <- ncol(data);
    fpkm.stats <- data.frame(
        Genes = row.names(data),
        Median.by.row = rep(NA,number.rows),
        IQR.by.row = rep(NA,number.rows)
        );
    for (i in 1:number.rows) {
        fpkm.stats[i,2]<- median(as.numeric(c(t(data[i,]))),na.rm = TRUE);
        fpkm.stats[i,3]<- IQR(as.numeric(c(t(data[i,]))),na.rm = TRUE) * 0.7413;
        }
    #empty df for results
    median.norm.fpkm <- data.frame(
        matrix(ncol = number.cols,
               nrow = number.rows,
               dimnames = list(row.names(data),colnames(data))));
    for (r in 1:number.rows) {
        for (c in 1:number.cols) {
            median.norm.fpkm[r,c] <- (as.numeric(data[r,c]) - fpkm.stats[r,2]) / fpkm.stats[r,3]
            }
        }
    #result <- list(fpkm.stats,median.norm.fpkm)
    return(median.norm.fpkm)
    }

### Data Analysis #################################################################################
sub.fpkm.TCGA <- fpkm.tumor.symbol.filter[rownames(outlier.patient.tag.01.t.p.order),];
sub.fpkm.TCGA <-sub.fpkm.TCGA[,colnames(outlier.patient.tag.01.t.p.order)];

### Fixing Column names in fpkm data to match clinical id
colnames(sub.fpkm.TCGA) <- substr(
    x = colnames(sub.fpkm.TCGA),
    start = 1,
    stop = nchar(colnames(sub.fpkm.TCGA))-4
    );
#replacing '.' with '-'
colnames(sub.fpkm.TCGA) <- gsub(
    pattern = '\\.',
    replacement = '-',
    x = colnames(sub.fpkm.TCGA)
    );

sub.fpkm.META <- fpkm.tumor.symbol.filter.symbol[rownames(outlier.patient.tag.01),];
sub.fpkm.META <- sub.fpkm.META[,colnames(outlier.patient.tag.01)];

#mean TCGA
mean.sub.TCGA.zscore.fpkm <- fpkm.normalization.zscore(data = sub.fpkm.TCGA);
saveRDS(mean.sub.TCGA.fpkm,file = 'TCGAZscore.RDS')
#mean.sub.TCGA.zscore.fpkm <- readRDS(
    #"C:/Users/jre83/OneDrive/UCLA B.I.G. Summer/Paul_Lab/Survival_Analysis_TCGA-BRCA/TCGAZscore.RDS")

#median TCGA
median.sub.TCGA.robust.zscore.fpkm <- fpkm.robust.zscore(data = sub.fpkm.TCGA);
saveRDS(median.META.robust,file = 'TCGArobustzscore.RDS')
#median.sub.TCGA.robust.zscore.fpkm <- readRDS(
    #"C:/Users/jre83/OneDrive/UCLA B.I.G. Summer/Paul_Lab/Survival_Analysis_TCGA-BRCA/TCGArobustzscore.RDS")

#mean META
mean.sub.META.zscore.fpkm <- fpkm.normalization.zscore(data = sub.fpkm.META)
saveRDS(median.META.robust,file = 'METAzscore.RDS')
#mean.sub.META.zscore.fpkm <- readRDS(
    #"C:/Users/jre83/OneDrive/UCLA B.I.G. Summer/Paul_Lab/Survival_Analysis_TCGA-BRCA/METAzscore.RDS")

#median META
median.sub.META.robust.zscore.fpkm <- fpkm.robust.zscore(data = sub.fpkm.META)
saveRDS(median.META.robust,file = 'METArobustzscore.RDS')
#median.sub.META.robust.zscore.fpkm <- readRDS(
    #"C:/Users/jre83/OneDrive/UCLA B.I.G. Summer/Paul_Lab/Survival_Analysis_TCGA-BRCA/METArobustzscore.RDS")




save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Normalization-Function',
    extension = 'txt')
    );