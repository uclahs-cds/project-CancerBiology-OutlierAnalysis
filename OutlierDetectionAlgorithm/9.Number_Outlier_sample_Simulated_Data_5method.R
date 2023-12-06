#!/usr/bin/env Rscript

### 9.Number_Outlier_sample_Simulated_Data_5method.R ####################################################
# Compute the 5 statistics of simulated data with patient removal

# Run parallel: 10 chucnks
args <- commandArgs(trailingOnly = TRUE)

# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of dataset
dataset.name <- 'CCLE';

# load the R environment file saved from 4.Simulated_Data_generation_2.R and 2.Distribution_Identfication.R
load(file = paste('4.Simulated_Data_generation_2.', dataset.name, '.', args, '.rda', sep = ''));
load(file = paste('5.Simulated_Data_5method.', dataset.name, '.', args, '.short.rda', sep = ''));


# Required R package
install.packages('extraDistr', repo = 'http://cran.us.r-project.org');
install.packages('truncnorm', repo = 'http://cran.us.r-project.org');
install.packages('SnowballC', repo = 'http://cran.us.r-project.org');
install.packages('lsa', repo = 'http://cran.us.r-project.org');
library(extraDistr);
library(truncnorm);
library(SnowballC);
library(lsa);
install.packages('parallel', repo = 'http://cran.us.r-project.org');
install.packages('foreach', repo = 'http://cran.us.r-project.org');
install.packages('doParallel', repo = 'http://cran.us.r-project.org');
library(parallel);
library(foreach);
library(doParallel);








# Manually set the number of patients to be excluded
#    - First round should be '1'
patient.arg <- 1;
patient.part.arg <- patient.part[1:(length(patient.part)-patient.arg )];
sample.number <- patient.part.arg;

# Remove the patients from the simulated data
negative.simulated.sum.arg <- negative.simulated.sum[,patient.part.arg];
rownames(negative.simulated.sum.arg) <- rownames(negative.simulated.sum);

rm(negative.simulated.sum);








# Define a minimum value
random.col <- sample(patient.part, 1)
decimal.number.max <- lapply(na.omit(fpkm.tumor.symbol.filter[,random.col]), function(x) {
    decimal.numbers <- sapply(x, function(y) {
        nchar(as.character(y)) - nchar(as.integer(y)) - 1
        })
    return(decimal.numbers)
    })    
add.minimum.value <- 1 / 10^as.numeric(max(unlist(decimal.number.max)));



# function: Compute the cosine similarity of the largest data point
cosine.similarity.large.value.percent <- function(x, y, large.value.percent) {

    # rounding function
    roundToInteger <- function(z) round(z, digits = 0)

    # check if large value percent is zero
    if (0 == large.value.percent) {
        large.value.number.integer <- 1;
        }
    else {
        large.value.number <- length(x) * (large.value.percent/100);
        large.value.number.integer <- roundToInteger(large.value.number);
        }
    
    # subset the largest values
    patient.larger.value <- (length(x)-large.value.number.integer + 1):length(x);
    observed.value <- sort(y);
    theoretical.value <- sort(x);
    mid.value <- c(1, 1);
    value.x.y <- data.frame(theoretical.value, observed.value);

    # calculate cosine similarity
    cosine.large.value <- NULL;
    cosine.large.value <- sapply(patient.larger.value, function(i) {
        cosine(as.numeric(value.x.y[i,]), c(1, 1))
        })
    cosine.large.value;
    }



# function: Trim 5% of samples from each side
trim.sample <- function(x, trim.portion = 5) {
    if (length(x) <= 10) {
        patient.trim.value <- 2:(length(x)-1);
    } else {
        trim.sample.number <- length(x) * (trim.portion/100);
        trim.sample.number.integer <- round(trim.sample.number, digits = 0);
        patient.trim.value <- (trim.sample.number.integer + 1):(length(x)-trim.sample.number.integer);
        }
    patient.trim.value;
    }

outlier.detection.cosine <- function (x, value.portion = 1) {

        # Define a minimum value
    decimal.number.max <- lapply(na.omit(x), function(x) {
        decimal.numbers <- sapply(x, function(y) {
            nchar(as.character(y)) - nchar(as.integer(y)) - 1
            })
        return(decimal.numbers)
        })    
    add.minimum.value <- 1 / 10^as.numeric(max(unlist(decimal.number.max)));


    trim.sample <- function(x, trim.portion = 5) {
        if (length(x) <= 10) {
            patient.trim.value <- 2:(length(x)-1);
        } else {
            trim.sample.number <- length(x) * (trim.portion/100);
            trim.sample.number.integer <- round(trim.sample.number, digits = 0);
            patient.trim.value <- (trim.sample.number.integer + 1):(length(x)-trim.sample.number.integer);
            }
        patient.trim.value;
        }

    # function: Compute the cosine similarity of the largest data point
    cosine.similarity.large.value.percent <- function(x, y, large.value.percent) {

        # rounding function
        roundToInteger <- function(z) round(z, digits = 0)

        # check if large value percent is zero
        if (0 == large.value.percent) {
            large.value.number.integer <- 1;
            }
        else {
            large.value.number <- length(x) * (large.value.percent/100);
            large.value.number.integer <- roundToInteger(large.value.number);
            }
        
        # subset the largest values
        patient.larger.value <- (length(x)-large.value.number.integer + 1):length(x);
        observed.value <- sort(y);
        theoretical.value <- sort(x);
        mid.value <- c(1, 1);
        value.x.y <- data.frame(theoretical.value, observed.value);

        # calculate cosine similarity
        cosine.large.value <- NULL;
        cosine.large.value <- sapply(patient.larger.value, function(i) {
            cosine(as.numeric(value.x.y[i,]), c(1, 1))
            })
        cosine.large.value;
        }

    sample.fpkm.qq <- na.omit(as.numeric(x[sample.number]))
    sample.fpkm.qq.nozero <- sample.fpkm.qq + add.minimum.value;
    

    # Trimmed samples -Trim 5% of each side
    sample.trim.number <- trim.sample(seq(length(sample.fpkm.qq.nozero)), 5);
    sample.fpkm.qq.trim <- sort(sample.fpkm.qq)[sample.trim.number];
    sample.fpkm.qq.nozero.trim <- sample.fpkm.qq.trim + add.minimum.value;

    
    # Quantile
    p <- ppoints(sample.fpkm.qq.nozero);
    
    # Distribution
    distribution.fit <- as.numeric(x[length(x)]);
    
    if (1 == distribution.fit){
        # 1. Normal distribution
        norm.mean <- mean(sample.fpkm.qq.nozero.trim);
        norm.sd <- sd(sample.fpkm.qq.nozero.trim);
        # Use truncated norm
        norm.quantiles <- qtruncnorm(p, a=0, b=Inf, mean = norm.mean, sd = norm.sd);
        obs.quantile.norm <- quantile(sample.fpkm.qq.nozero, prob = p);
        last.cos <- cosine.similarity.large.value.percent(norm.quantiles, obs.quantile.norm, large.value.percent = value.portion);
        }
    else if (2 == distribution.fit) {
        # 2. Log-normal distribution
        mean.log <- mean(sample.fpkm.qq.nozero.trim);
        sd.log <- sd(sample.fpkm.qq.nozero.trim);
        m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
        sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
        lnorm.quantile <- qlnorm(p, meanlog = m2, sdlog = sd2);
        obs.quantile.lnorm <- quantile(sample.fpkm.qq.nozero, prob = p);
        last.cos <- cosine.similarity.large.value.percent(lnorm.quantile, obs.quantile.lnorm, large.value.percent = value.portion);
        }
    else if (3 == distribution.fit) {
        # 3. Exponential distribution
        exp.rate <- 1 / mean(sample.fpkm.qq.nozero.trim);
        exp.quantile <- qexp(p, rate = exp.rate);
        obs.quantile.exp <- quantile(sample.fpkm.qq.nozero, prob = p);
        last.cos <- cosine.similarity.large.value.percent(exp.quantile, obs.quantile.exp, large.value.percent = value.portion);
        }
    else if (4 == distribution.fit) {
        ### 4 gamma distribution
        mean.gamma <- mean(sample.fpkm.qq.nozero.trim);
        sd.gamma <- sd(sample.fpkm.qq.nozero.trim);
        gamma.shape <- (mean.gamma/sd.gamma)^2;
        gamma.rate <- mean.gamma/(sd.gamma^2);
        gamma.quantile <- qgamma(p, shape = gamma.shape, rate = gamma.rate);
        obs.quantile.gamma <- quantile(sample.fpkm.qq.nozero, prob = p);
        last.cos <- cosine.similarity.large.value.percent(gamma.quantile, obs.quantile.gamma, large.value.percent = value.portion);
        }

    cosine.sum.distribution.fit <- c(last.cos, distribution.fit);
    cosine.sum.distribution.fit;
    }




# Check the cosine similarity
negative.simulated.sum.fit <- cbind(negative.simulated.sum.arg, distribution = data.cosine.negative.t$distribution);
# run it parallel
cl <- makeCluster(20);
# register the cluster with the parallel package
registerDoParallel(cl);
clusterExport(cl, "outlier.detection.cosine");
clusterEvalQ(cl, c(library(lsa), library(SnowballC)));

data.cosine.negative <- apply(negative.simulated.sum.fit, 
                         1, 
                         outlier.detection.cosine, 
                         value.portion = 0);

stopImplicitCluster();


data.cosine.negative.t <- t(data.cosine.negative);
data.cosine.negative.t.arg <- data.frame(data.cosine.negative.t);
colnames(data.cosine.negative.t.arg) <- c('cosine', 'distribution');




# 1,2,3,4 
quantify.outliers <- function(x, methods = 'mean', trim = 0, exclude.zero = FALSE) {
    x.na <- na.omit(as.numeric(x));
    if (methods == 'median') {
        if (exclude.zero) { 
            x.nonzero <- x.na[0 != x.na]; 
            data.median <- median(x.nonzero);
            data.mad <- mad(x.nonzero);
            } 
        else {
            data.median <- median(x.na);
            data.mad <- mad(x.na);
            }
        result.na <- (x.na - data.median) / data.mad;
        x[which(!is.na(x))] <- result.na;
        x;
        }
    else if (methods == 'kmean') {
        if (exclude.zero) {
            if (length(unique(as.numeric(x.na))) == 1) {
                kmean.matrix <- rep(NA, length(x.na));
                names(kmean.matrix) <- names(x.na);
                } 
            else {
                data.order <- sort(x.na, decreasing = TRUE);
                non.zero <- data.order[data.order > 0];
                if (length(unique(as.numeric(non.zero))) <= 2) {
                    na.matrix <- rep(NA, length(non.zero));
                    cluster.zero <- c(na.matrix, rep(0, length(x.na[x.na == 0])));
                    kmean.matrix <- cluster.zero[match(x.na, data.order)];
                    names(kmean.matrix) <- names(x.na);  
                    } 
                else {
                    kmean <- kmeans(non.zero, 2, nstart = 1000);
                    cluster <- kmean$cluster;
                    cluster.zero <- c(cluster, rep(0, length(x[x == 0])));
                    kmean.matrix <- cluster.zero[match(x.na, data.order)];
                    names(kmean.matrix) <- names(x.na);   
                    }
                }
            } 
    
        else {
            if (length(unique(as.numeric(x.na))) == 1) {
                kmean.matrix <- rep(NA, length(x.na));
                names(kmean.matrix) <- names(x.na);  
                } 
            else {
                kmean <- kmeans(x.na, 2, nstart = 1000);
                cluster <- kmean$cluster;
                kmean.matrix <- cluster;
                names(kmean.matrix) <- names(x.na);  
                }
            }
        result.na <- kmean.matrix;
        x[which(!is.na(x))] <- result.na;
        x;
        }
    else {
        gene.order <- x.na[order(x.na, decreasing = TRUE)];
        if (exclude.zero) { 
            gene.order.nonzero <- gene.order[0 != gene.order]; 
            top.patient <- round(length(gene.order.nonzero) * (trim / 100), digit = 0);
            low.patient <- round(length(gene.order.nonzero) * (1 - (trim / 100)), digit = 0);
            data.mean <- mean(gene.order.nonzero, trim = (trim / 100));
            data.sd <- sd(gene.order.nonzero[(top.patient+1):(low.patient)]);
            } 
        else {
            top.patient <- round(length(x.na) * (trim / 100), digit = 0);
            low.patient <- round(length(x.na) * (1 - (trim / 100)), digit = 0);
            data.mean <- mean(gene.order, trim = (trim / 100));
            data.sd <- sd(gene.order[(top.patient+1):(low.patient)]);
            }
        result.na <- (x.na - data.mean) / data.sd;
        x[which(!is.na(x))] <- result.na;
        x;
        }
    }






# Parallel running
cl <- makeCluster(20);
# register the cluster with the parallel package
registerDoParallel(cl);

# 1. MEAN and SD : method = 'mean', trim = 0
data.mean <- foreach(i=1:nrow(negative.simulated.sum.arg), .combine = rbind) %dopar% quantify.outliers(negative.simulated.sum.arg[i,]);
data.mean <- data.frame(data.mean);

# 2. TRIMMED MEAN and TRIMMED SD : method = 'mean', trim = 5
data.trimmean <- foreach(i=1:nrow(negative.simulated.sum.arg), .combine = rbind) %dopar% quantify.outliers(negative.simulated.sum.arg[i,], trim = 5);
data.trimmean <- data.frame(data.trimmean);

# 3. MEDIAN and MAD : method = 'median'
data.median <- foreach(i=1:nrow(negative.simulated.sum.arg), .combine = rbind) %dopar% quantify.outliers(negative.simulated.sum.arg[i,], methods = 'median');
data.median <- data.frame(data.median);

# 4. KMEAN : method = 'kmean'
data.kmean <- foreach(i=1:nrow(negative.simulated.sum.arg), .combine = rbind) %dopar% quantify.outliers(negative.simulated.sum.arg[i,], methods = 'kmean')
data.kmean <- data.frame(data.kmean);

stopCluster(cl)



outlier.detection.zrange <- function(x) {
  x.na <- na.omit(x)
  zrange <- max(x.na) - min(x.na);
  zrange.matrix <- c(x, zrange);
  names(zrange.matrix) <- c(names(x), 'zrange');
  zrange.matrix;
    }




# 1. MEAN and SD
data.zrange.mean <- apply(data.mean, 1, outlier.detection.zrange);
mean.simulated.negative.1M.arg <- data.frame(t(data.zrange.mean));

# 2. TRIMMED MEAN and TRIMMED SD
data.zrange.trimmean <- apply(data.trimmean, 1, outlier.detection.zrange);
trimmean.simulated.negative.1M.arg <- data.frame(t(data.zrange.trimmean));

# 4. MEDIAN and MAD
data.zrange.median <- apply(data.median, 1, outlier.detection.zrange);
median.simulated.negative.1M.arg<- data.frame(t(data.zrange.median));




### Calculate the kmean fraction #####
# Function
outlier.detection.kmean <- function(x) {
    if (1== length(unique(as.numeric(x)))) {
        fraction <- NA;
        }
    else {
        cluster.one <- length(x[x == 1]);
        cluster.two <- length(x[x == 2]);
        cluster.sum <- cluster.one + cluster.two;
        smaller.value <- min(cluster.one, cluster.two);
        fraction <- round(smaller.value/cluster.sum, digit = 4);
        }
    fraction.matrix <- c(x, fraction);
    names(fraction.matrix) <- c(names(x), 'fraction');
    fraction.matrix;
    }

# 4. KMEAN fraction
data.fraction.kmean <- apply(data.kmean, 1, outlier.detection.kmean);
kmean.simulated.negative.1M.arg <- data.frame(t(data.fraction.kmean));




### Final gene-wise matrix #####
gene.zrange.fraction.negative.simulated.sum.1M.arg <- cbind(mean.simulated.negative.1M.arg$zrange,
                     median.simulated.negative.1M.arg$zrange,
                     trimmean.simulated.negative.1M.arg$zrange,
                     kmean.simulated.negative.1M.arg$fraction);
rownames(gene.zrange.fraction.negative.simulated.sum.1M.arg) <- rownames(negative.simulated.sum.arg);
colnames(gene.zrange.fraction.negative.simulated.sum.1M.arg) <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean');


# Final statistic matrix
gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.arg <- cbind(gene.zrange.fraction.negative.simulated.sum.1M.arg[,c(1,2,3,4)],
                                                    data.cosine.negative.t.arg$cosine,
                                                    data.cosine.negative.t.arg$distribution);
colnames(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.arg) <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean', 'cosine', 'distribution');





zrange.mean <- paste('mean.simulated.negative.1M.', args, sep = '');
assign(zrange.mean, mean.simulated.negative.1M.arg);

zrange.trimmean <- paste('trimmean.simulated.negative.1M.', args, sep = '');
assign(zrange.trimmean, trimmean.simulated.negative.1M.arg);

zrange.median <- paste('median.simulated.negative.1M.', args, sep = '');
assign(zrange.median, median.simulated.negative.1M.arg);

fraction.kmean <- paste('kmean.simulated.negative.1M.', args, sep = '');
assign(fraction.kmean, kmean.simulated.negative.1M.arg);

cosine.bic <- paste('data.cosine.negative.t.', args, sep = '');
assign(cosine.bic, data.cosine.negative.t.arg);

gene.zrange.fraction <- paste('gene.zrange.fraction.negative.simulated.sum.1M.', args, sep = '');
assign(gene.zrange.fraction, gene.zrange.fraction.negative.simulated.sum.1M.arg);

gene.zrange.fraction.cosine <- paste('gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.', args, sep = '');
assign(gene.zrange.fraction.cosine, gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.arg);





save(
    fpkm.tumor.symbol.filter,
    sample.number,
    bic.trim.distribution.fit.obs,
    bic.trim.distribution.fit,
    list = c(
        paste0('gene.zrange.fraction.negative.simulated.sum.1M.', args, sep = ''),
        paste0('data.cosine.negative.t.', args, sep = ''),
        paste0('gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.', args, sep = '')),
    file = paste('9.Number_Outlier_sample_Simulated_Data_5method.', dataset.name, '.', args, '.', patient.arg, '.short.rda', sep = '')
    );


save(
    fpkm.tumor.symbol.filter,
    patient.part,
    bic.trim.distribution.fit.obs,
    bic.trim.distribution.fit,
    list = c(
        paste0('mean.simulated.negative.1M.', args, sep = ''),
        paste0('median.simulated.negative.1M.', args, sep = ''),
        paste0('trimmean.simulated.negative.1M.', args, sep = ''),
        paste0('kmean.simulated.negative.1M.', args, sep = ''),
        paste0('gene.zrange.fraction.negative.simulated.sum.1M.', args, sep = ''),
        paste0('data.cosine.negative.t.', args, sep = ''),
        paste0('gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.', args, sep = '')),
    file = paste('9.Number_Outlier_sample_Simulated_Data_5method.', dataset.name, '.', args, '.', patient.arg, '.long.rda', sep = '')
    );


