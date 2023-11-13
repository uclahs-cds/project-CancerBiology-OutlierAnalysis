#!/usr/bin/env Rscript

### 4.Simulated_Data_generation_2.R ####################################################
# Generate the simulated data: part2

# Run parallel: 10 chucnks
args <- commandArgs(trailingOnly = TRUE);

# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of dataset
dataset.name <- 'CCLE';

# load the R environment file saved from 2.Distribution_Identification.R and 3.Simulated_Data_generation_1.R
#load(file = '2023-09-06_CCLE_final_outlier_rank_bic_distribution.short.rda');
load(file = paste('2.Distribution_Identification.', dataset.name, '.short.rda', sep = ''));

# load(file = paste('2023_09-06_parallel_null_distribution_100k_45sample_bic_only1_', args, '.rda', sep = ''));
load(file = paste('3.Simulated_Data_generation_1.', dataset.name, '.', args, '.rda', sep = ''));

# Required R package
install.packages('extraDistr', repo = 'http://cran.us.r-project.org');
install.packages('truncnorm', repo = 'http://cran.us.r-project.org');
library(extraDistr);
library(truncnorm);
install.packages('parallel', repo = 'http://cran.us.r-project.org');
install.packages('foreach', repo = 'http://cran.us.r-project.org');
install.packages('doParallel', repo = 'http://cran.us.r-project.org');
library(parallel);
library(foreach);
library(doParallel);






# sample size
patient.part <- 1:ncol(fpkm.tumor.symbol.filter);
sample.number <- 1:ncol(fpkm.tumor.symbol.filter);
molecular.data.filter <- fpkm.tumor.symbol.filter[,patient.part];


fpkm.tumor.symbol.filter.bic <- fpkm.tumor.symbol.filter[match(names(bic.trim.distribution.fit), rownames(fpkm.tumor.symbol.filter)),];


random.col <- sample(patient.part, 1)
decimal.number.max <- lapply(na.omit(fpkm.tumor.symbol.filter[,random.col]), function(x) {
    decimal.numbers <- sapply(x, function(y) {
        nchar(as.character(y)) - nchar(as.integer(y)) - 1
        })
    return(decimal.numbers)
    })    
add.minimum.value <- 1 / 10^as.numeric(max(unlist(decimal.number.max)));


# 2. residual
residue.negative.random.number.bic <- obs.residue.quantile.trim[match(substr(rownames(negative.random.number.bic), 1, 15), substr(rownames(obs.residue.quantile.trim), 1, 15)),];

noise.min.off.bic.distribution.residue <- noise.min.off.bic.distribution.fit[match(substr(rownames(negative.random.number.bic), 1, 15), substr(names(noise.min.off.bic.distribution.fit), 1, 15))];


# create a cluster of CPU cores to use in parallel processing
cl <- makeCluster(10)
# register the cluster with the parallel package
registerDoParallel(cl);
clusterExport(cl, "add.minimum.value");
clusterEvalQ(cl, library(extraDistr));

random.col <- sample(patient.part, 1)
decimal.number.max <- lapply(na.omit(fpkm.tumor.symbol.filter[,random.col]), function(x) {
    decimal.numbers <- sapply(x, function(y) {
        nchar(as.character(y)) - nchar(as.integer(y)) - 1
        })
    return(decimal.numbers)
    })    
add.minimum.value <- 1 / 10^as.numeric(max(unlist(decimal.number.max)));

negative.random.number.noise.bic <- foreach(i=1:nrow(residue.negative.random.number.bic), .combine='rbind') %dopar% {

    sample.fpkm.qq <- as.numeric(residue.negative.random.number.bic[i,]);
    sample.fpkm.qq.sort <- sort(sample.fpkm.qq);
    if (min(sample.fpkm.qq.sort) < 0) {
        sample.fpkm.qq.nozero <- sample.fpkm.qq.sort - min(sample.fpkm.qq.sort) + add.minimum.value;
        }
    else {
        sample.fpkm.qq.nozero <- sample.fpkm.qq.sort + add.minimum.value;
        }    

        if (1 == noise.min.off.bic.distribution.residue[i]){
            ### 1) Normal distribution
            norm.mean <- mean(sample.fpkm.qq.nozero);
            norm.sd <- sd(sample.fpkm.qq.nozero);
            simulated.sample <- rtnorm(length(sample.number), mean = norm.mean, sd = norm.sd, a = 0);
            }
        else if (2 == noise.min.off.bic.distribution.residue[i]){
            mean.log <- mean(sample.fpkm.qq.nozero);
            sd.log <- sd(sample.fpkm.qq.nozero);
            m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
            sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
            simulated.sample <- rlnorm(n=length(sample.number), m2, sd2);
            }
        else if (3 == noise.min.off.bic.distribution.residue[i]){
            ### 4) exponential distribution
            exp.rate <- 1 / mean(sample.fpkm.qq.nozero);
            simulated.sample <- rexp(n=length(sample.number), rate = exp.rate);
            }
        else if (4 == noise.min.off.bic.distribution.residue[i]){
            ### 5) gamma distribution
           mean.gamma <- mean(sample.fpkm.qq.nozero);
           sd.gamma <- sd(sample.fpkm.qq.nozero);
           gamma.shape <- (mean.gamma/sd.gamma)^2;
           gamma.rate <- mean.gamma/(sd.gamma^2);
           simulated.sample <- rgamma(n=length(sample.number), gamma.shape, gamma.rate);
            }
        
        if (min(sample.fpkm.qq.sort) < 0) {
            simulated.sample.min <- simulated.sample + min(sample.fpkm.qq.sort) - add.minimum.value;
            }
        else {
            simulated.sample.min <- simulated.sample - add.minimum.value; 
            }
        simulated.sample.min;
        
        }

rownames(negative.random.number.noise.bic) <- rownames(residue.negative.random.number.bic);    


stopImplicitCluster()

# 3. sum
negative.simulated.sum <- abs(negative.random.number.bic + negative.random.number.noise.bic);


# save the R environment
save(
    fpkm.tumor.symbol.filter,
    patient.part,
    bic.trim.distribution.fit,
    negative.random.number.bic,
    obs.residue.quantile.trim,
    residue.negative.random.number.bic,
    noise.min.off.bic.distribution.fit,
    noise.min.off.bic.distribution.residue,
    negative.random.number.noise.bic,
    negative.simulated.sum,
    file = paste('4.Simulated_Data_generation_2.', dataset.name, '.', args, '.rda', sep = '')
    )






