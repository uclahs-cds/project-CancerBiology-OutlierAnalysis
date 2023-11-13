#!/usr/bin/env Rscript

### 3.Simulated_Data_generation_1.R ####################################################
# Generate the simulated data: part1


# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of dataset
dataset.name <- 'CCLE';

# load the R environment file saved from 2.Distribution_Identification.R
load(file = paste('2.Distribution_Identification.', dataset.name, '.long.rda', sep = ''));

# Required R packages
install.packages('extraDistr', repo = 'http://cran.us.r-project.org');
install.packages('truncnorm', repo = 'http://cran.us.r-project.org');
install.packages('SnowballC', repo = 'http://cran.us.r-project.org');
install.packages('lsa', repo = 'http://cran.us.r-project.org');
install.packages('parallel', repo = 'http://cran.us.r-project.org');
install.packages('foreach', repo = 'http://cran.us.r-project.org');
install.packages('doParallel', repo = 'http://cran.us.r-project.org');
library(extraDistr);
library(truncnorm);
library(SnowballC);
library(lsa);
library(parallel);
library(foreach);
library(doParallel);

# Run parallel: 10 chucnks
args <- commandArgs(trailingOnly = TRUE)

# sample size
patient.part <- 1:ncol(fpkm.tumor.symbol.filter);
sample.number <- 1:ncol(fpkm.tumor.symbol.filter);
molecular.data.filter <- fpkm.tumor.symbol.filter[,patient.part];
### Make null distribution (BIC)
# 1. random numbers

##### 1. Make random numbers
# - Using BIC
fpkm.tumor.symbol.filter.bic <- fpkm.tumor.symbol.filter[match(names(bic.trim.distribution.fit), rownames(fpkm.tumor.symbol.filter)),];



# create a cluster of CPU cores to use in parallel processing
cl <- makeCluster(detectCores()-2)
# register the cluster with the parallel package
registerDoParallel(cl)

clusterEvalQ(cl, library(extraDistr));
    
simulated.generation.negative <- function (x = fpkm.tumor.symbol.filter.distribution.fit, 
                                           distribution = bic.trim.distribution.fit, 
                                           num.negative = 10000,
                                           sample.size = sample.number) {
    
    # Define a minimum value
    random.col <- sample(sample.size, 1)
    decimal.number.max <- lapply(na.omit(x[,random.col]), function(x) {
        decimal.numbers <- sapply(x, function(y) {
            nchar(as.character(y)) - nchar(as.integer(y)) - 1
            })
        return(decimal.numbers)
        })    
    add.minimum.value <- 1 / 10^as.numeric(max(unlist(decimal.number.max)));
    
    # function: Trim 5% of samples from each side
    trim.sample <- function(x, trim.portion = 5) {
        trim.sample.number <- length(x) * (trim.portion/100);
        trim.sample.number.integer <- round(trim.sample.number, digits = 0);
        patient.trimr.value <- (trim.sample.number.integer + 1):(length(x)-trim.sample.number.integer);
        patient.trimr.value;
        }

    
        # Validate input parameters
    if(!is.data.frame(x)) stop("x should be a data frame.")
    if(!is.numeric(distribution)) stop("distribution should be numeric.")
    if(!is.numeric(num.negative)) stop("num.negative should be numeric.")
    if(!is.numeric(sample.size)) stop("sample.size should be numeric.")
    
    random.number.negative <- sample(length(distribution), num.negative, replace = TRUE);
    names(random.number.negative) <- names(distribution)[random.number.negative]
    
    # use the foreach function to parallelize the sapply loop
    simulated.negative <- foreach(i = random.number.negative, .combine = rbind) %dopar% {
        
        sample.fpkm <- x[names(distribution)[i], sample.size];
        sample.fpkm.qq <- na.omit(as.numeric(sample.fpkm));
        sample.fpkm.qq.nozero <- sample.fpkm.qq + add.minimum.value;
        
        # Trimmed samples -Trim 5% of each side
        sample.trim.number <- trim.sample(sample.fpkm.qq.nozero, 5);
        sample.fpkm.qq.trim <- sort(sample.fpkm.qq)[sample.trim.number];
        sample.fpkm.qq.nozero.trim <- sample.fpkm.qq.trim + add.minimum.value;
                                   
        if (1 == distribution[i]){
            ### 1) Normal distribution
            norm.mean <- mean(sample.fpkm.qq.nozero.trim);
            norm.sd <- sd(sample.fpkm.qq.nozero.trim);
            rtnorm(length(sample.size), mean = norm.mean, sd = norm.sd, a = 0);
            }
        else if (2 == distribution[i]){
            mean.log <- mean(sample.fpkm.qq.nozero.trim);
            sd.log <- sd(sample.fpkm.qq.nozero.trim);
            m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
            sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
            rlnorm(n=length(sample.size), m2, sd2);
        }
        else if (3 == distribution[i]){
            ### 4) exponential distribution
            exp.rate <- 1 / mean(sample.fpkm.qq.nozero.trim);
            rexp(n=length(sample.size), rate = exp.rate);
            }
        else if (4 == distribution[i]){
            ### 5) gamma distribution
           mean.gamma <- mean(sample.fpkm.qq.nozero.trim);
           sd.gamma <- sd(sample.fpkm.qq.nozero.trim);
           gamma.shape <- (mean.gamma/sd.gamma)^2;
           gamma.rate <- mean.gamma/(sd.gamma^2);
           rgamma(n=length(sample.size), gamma.shape, gamma.rate);
            }
        }

    rownames(simulated.negative) <- names(distribution)[random.number.negative];
    simulated.negative;
}


negative.random.number.bic <- simulated.generation.negative(x = fpkm.tumor.symbol.filter.bic, 
                                           distribution = bic.trim.distribution.fit, 
                                           num.negative = 100000,
                                           sample.size = sample.number);

stopImplicitCluster()


# save the R environment
save(
    fpkm.tumor.symbol.filter,
    patient.part,
    sample.number,
    bic.trim.distribution.fit,
    negative.random.number.bic,
    file = paste('3.Simulated_Data_generation_1.', dataset.name, '.', args, '.rda', sep = '')
    );




