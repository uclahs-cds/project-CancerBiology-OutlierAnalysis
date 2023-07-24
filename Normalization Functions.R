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
    #data
    # the name of the data frame that contains non-normalized fpkm values
#Output Variables
    #norm.fpkm
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
    #data:
    #the name of the data frame that contains non-normalized fpkm values
#Output Variables
    #median.norm.fpkm
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
#replicate data set for analysis
fpkm.nonnorm <- fpkm.tumor.symbol.filter;

### Fixing Column names in fpkm data to match clinical id
colnames(fpkm.nonnorm) <- substr(
    x = colnames(fpkm.nonnorm),
    start = 1,
    stop = nchar(colnames(fpkm.nonnorm))-4
    );
#replacing '.' with '-'
colnames(fpkm.nonnorm) <- gsub(
    pattern = '\\.',
    replacement = '-',
    x = colnames(fpkm.nonnorm)
    );

#mean
mean.zscore <- fpkm.normalization.zscore(fpkm.nonnorm);
#saveRDS(fpkm.norm,file = 'Z.score.RDS')

#median
median.robust.zscore <- fpkm.robust.zscore(fpkm.nonnorm);
#saveRDS(median.normalized.fpkm,file = "Robust-Median-Zscore.RDS")

save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Normalization-Function',
    extension = 'txt')
    );