#!/usr/bin/env Rscript

### 7.Significant_Outlier_Detection.R ####################################################
# Compute p-values 


# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of dataset
dataset.name <- 'CCLE';

# Load the R environment
#   - 1. File from script 1: short version
load(file = '2023-09-18_CCLE_final_outlier_rank_bic.short.rda');
#   - 2. File from script 6
load(file = paste('6.Simulated_Data_5method_combine.', dataset.name, '.rda', sep = ''));

# Required R package
install.packages('parallel', repo = 'http://cran.us.r-project.org');
install.packages('foreach', repo = 'http://cran.us.r-project.org');
install.packages('doParallel', repo = 'http://cran.us.r-project.org');
library(parallel);
library(foreach);
library(doParallel);


gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- data.frame(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M);


# Run 1000 genes at once
#   - array number should be ceiling(nrow(fpkm.tumor.symbol.filter))
args <- commandArgs(trailingOnly = TRUE);
row.num.args <- as.numeric(args);

# This will be used identify the number of outlier patients per gene
#   - if '0', use whole patients (first step)
#   - if '1', use n-1 patients (exclude the patient having the largest value)
#   - repeat this '2', '3', '4'... until there is no outlier genes
data.args <- 0;





### Rank each methods #####
# Function
outlier.rank <- function(x) {
    methods <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean', 'cosine');
    rank.matrix <- NULL;
    # Give rank for each methods based on z-score range/fraction of kmean
    for (i in 1:3) {
        methods.column <- methods[i];
        rank.methods <- rank(-x[,methods.column], ties.method = 'max', na.last = 'keep');
        rank.matrix <- cbind(rank.matrix, rank.methods);
        }
    for (i in 4:5) {
        methods.column <- methods[i];
        rank.methods <- rank(x[,methods.column], ties.method = 'max', na.last = 'keep');
        rank.matrix <- cbind(rank.matrix, rank.methods);
        }
    rownames(rank.matrix) <-rownames(x);
    colnames(rank.matrix) <- methods;
    rank.matrix <- data.frame(rank.matrix);
    }


### Rank product to determine Top ranked genes #####
# Function
# x: ranked matrix
# NA.number = Number of methods with non-NA should be more than assigned number
outlier.rank.product <- function(x, NA.number = 0) {
    rank <- as.numeric(x[1:5]);
    num <- length(which(!is.na(rank)));
    if (NA.number >= num) {
        NA;
        }
    else {
        prod(rank, na.rm = TRUE)^(1/num);
        }
    }



### Combine matrix
# - relabel the null data
gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relable <- cbind(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M,
                                                                            gene = rownames(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M));
rownames(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relable) <- paste(rep("ND", 1000000), c(1:1000000), sep = '');


# Assign the row number from start to end
gene.number.start.end.matrix <- NULL;
for (i in 1:ceiling(nrow(fpkm.tumor.symbol.filter)/1000)) {
    if (i == ceiling(nrow(fpkm.tumor.symbol.filter)/1000)) {
        range.number <- c((i-1)*1000 + 1, nrow(fpkm.tumor.symbol.filter));
        gene.number.start.end.matrix <- data.frame(rbind(
            gene.number.start.end.matrix,
            i = range.number
            ));
        } 
    else {
        range.number<- c((i-1)*1000 + 1, i*1000)
        gene.number.start.end.matrix <- data.frame(rbind(
            gene.number.start.end.matrix,
            i = range.number
            ));
        }
    }



cl <- makeCluster(20);
# register the cluster with the parallel package
registerDoParallel(cl);
clusterExport(cl, c("outlier.rank", "outlier.rank.product"));

gene.zrange.fraction.fpkm.bic.5method.1M.data <- get(paste('gene.zrange.fraction.cosine.last.point.bic', sep = ''));

gene.rank.p.value.one.gene <- NULL;
gene.rank.p.value.one.gene <- foreach(i = as.numeric(gene.number.start.end.matrix[row.num.args,1]):as.numeric(gene.number.start.end.matrix[row.num.args,2]), .combine=rbind) %dopar% {
  observed.gene <-  gene.zrange.fraction.fpkm.bic.5method.1M.data[i,1:5];
  combine.matrix <- rbind(observed.gene,
                          gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relable[,1:5]);
  # get ranks
  data.rank.bic <- outlier.rank(combine.matrix);
  rank.product.bic <- apply(data.rank.bic, 1, outlier.rank.product, NA.number = 3);
  gene.rank.poduct.bic <- cbind(data.rank.bic,
                                rank.product.bic);
  obs <- rank.product.bic[1]
  null <- rank.product.bic[2:(nrow(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relable) + 1)]
  length.null <- nrow(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relable);
  obs.p.value <- (sum(obs >= null) + 1) / (length.null + 1)
  
  obs.p.value.rank <- cbind(gene.rank.poduct.bic[1,], obs.p.value);
  p.value.one.gene <- data.frame(x = obs.p.value.rank, i = i);
  p.value.one.gene;
}


p.value.one <- paste('gene.rank.p.value.one.gene.', data.args, sep = '');
assign(p.value.one, gene.rank.p.value.one.gene);


stopImplicitCluster();




save(
    list = paste0('gene.rank.p.value.one.gene.', data.args, sep = ''),
    file = paste('7.Significant_Outlier_Detection.', dataset.name, '.', row.num.args, '.',data.args, '.rda', sep = '')
    );


