#!/usr/bin/env Rscript


### 10.Number_Outlier_sample_5method.R ####################################################
# Basically, this is the same script as 1.Outlier_detection.R
#   - exclude the largest value (exclude one patient) and run the outlier detection algorithm


# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of the dataset
dataset.name <- 'CCLE';

# Required R packages
# Install and load the 'gamlss' package
install.packages('gamlss', repo = 'http://cran.us.r-project.org');
library(gamlss);
# Install and load the 'doParallel' package
install.packages('doParallel', repo = 'http://cran.us.r-project.org');
library(doParallel);
# Install and load the 'foreach' package
install.packages('foreach', repo = 'http://cran.us.r-project.org');
library(foreach);
# Install and load the 'parallel' package
install.packages('parallel', repo = 'http://cran.us.r-project.org');
library(parallel);
# Install and load the 'extraDistr' package
install.packages('extraDistr', repo = 'http://cran.us.r-project.org');
library(extraDistr);
# Install and load the 'truncnorm' package
install.packages('truncnorm', repo = 'http://cran.us.r-project.org');
library(truncnorm);
# Install and load the 'lsa' package
install.packages('lsa', repo = 'http://cran.us.r-project.org');
library(lsa);
# Install and load the 'SnowballC' package
install.packages('SnowballC', repo = 'http://cran.us.r-project.org');
library(SnowballC);


# Load the R environment
#   - 1. File from script 1: short version
load(file = '2023-09-18_CCLE_final_outlier_rank_bic.short.rda');


# Set the number of patients to be excluded
#    - First round should be '1'
args <- commandArgs(trailingOnly = TRUE);
args.num <- as.numeric(args);

# Sample number
patient.part.arg <- patient.part[1:(length(patient.part)-args.num)];
sample.number <- patient.part.arg;

# Remove the largest value of each gene
fpkm.tumor.symbol.filter.arg <- NULL;
for(i in 1:nrow(fpkm.tumor.symbol.filter)) {
    fpkm.sort.1 <- sort(as.numeric(fpkm.tumor.symbol.filter[i,patient.part]))[seq(length(patient.part.arg))];
    fpkm.tumor.symbol.filter.arg <- rbind(fpkm.tumor.symbol.filter.arg, fpkm.sort.1);
    }
rownames(fpkm.tumor.symbol.filter.arg) <- rownames(fpkm.tumor.symbol.filter);



# Same script from #1
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





outlier.detection.cosine <- function (x, value.portion = 1) {

        # Define a minimum value
    decimal.number.max <- lapply(na.omit(x), function(x) {
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

    sample.fpkm.qq <- as.numeric(x[sample.number])
    sample.fpkm.qq.nozero <- sample.fpkm.qq + add.minimum.value;
    

    # Trimmed samples -Trim 5% of each side
    sample.trim.number <- trim.sample(sample.number, 5);
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


# Determine the distribution
# Find the best fitted distribution

trim.sample <- function(x, trim.portion = 5) {
    trim.sample.number <- length(x) * (trim.portion/100);
    trim.sample.number.integer <- round(trim.sample.number, digits = 0);
    patient.trimr.value <- (trim.sample.number.integer + 1):(length(x)-trim.sample.number.integer);
    patient.trimr.value;
    }





cl <- makeCluster(25);
# register the cluster with the parallel package
registerDoParallel(cl);
clusterExport(cl, "outlier.detection.cosine");
clusterEvalQ(cl, c(library(lsa), library(SnowballC)));

fpkm.tumor.symbol.filter.bic.fit <- cbind(fpkm.tumor.symbol.filter.arg, distribution = bic.trim.distribution.fit);
data.cosine.bic <- apply(fpkm.tumor.symbol.filter.bic.fit, 
                         1, 
                         outlier.detection.cosine, 
                         value.portion = 0);

stopImplicitCluster();


data.cosine.bic.t.arg <- data.frame(t(data.cosine.bic));
colnames(data.cosine.bic.t.arg) <- c('cosine', 'distribution');




# 1,2,3,4 
quantify.outliers <- function(x, methods = 'mean', trim = 0, exclude.zero = FALSE) {
    x <- as.numeric(x);
    if (methods == 'median') {
        if (exclude.zero) { 
            x.nonzero <- x[0 != x]; 
            data.median <- median(x.nonzero);
            data.mad <- mad(x.nonzero);
            } 
        else {
            data.median <- median(x);
            data.mad <- mad(x);
            }
        (x - data.median) / data.mad;
        }
    else if (methods == 'kmean') {
        if (exclude.zero) {
            if (length(unique(as.numeric(x))) == 1) {
                kmean.matrix <- rep(NA, length(x));
                names(kmean.matrix) <- names(x);
                } 
            else {
                data.order <- sort(x, decreasing = TRUE);
                non.zero <- data.order[data.order > 0];
                if (length(unique(as.numeric(non.zero))) <= 2) {
                    na.matrix <- rep(NA, length(non.zero));
                    cluster.zero <- c(na.matrix, rep(0, length(x[x == 0])));
                    kmean.matrix <- cluster.zero[match(x, data.order)];
                    names(kmean.matrix) <- names(x);  
                    } 
                else {
                    kmean <- kmeans(non.zero, 2, nstart = 1000);
                    cluster <- kmean$cluster;
                    cluster.zero <- c(cluster, rep(0, length(x[x == 0])));
                    kmean.matrix <- cluster.zero[match(x, data.order)];
                    names(kmean.matrix) <- names(x);   
                    }
                }
            } 
    
        else {
            if (length(unique(as.numeric(x))) == 1) {
                kmean.matrix <- rep(NA, length(x));
                names(kmean.matrix) <- names(x);  
                } 
            else {
                kmean <- kmeans(x, 2, nstart = 1000);
                cluster <- kmean$cluster;
                kmean.matrix <- cluster;
                names(kmean.matrix) <- names(x);  
                }
            }
        kmean.matrix;
        }
    else {
        gene.order <- x[order(x, decreasing = TRUE)];
        if (exclude.zero) { 
            gene.order.nonzero <- gene.order[0 != gene.order]; 
            top.patient <- round(length(gene.order.nonzero) * (trim / 100), digit = 0);
            low.patient <- round(length(gene.order.nonzero) * (1 - (trim / 100)), digit = 0);
            data.mean <- mean(gene.order.nonzero, trim = (trim / 100));
            data.sd <- sd(gene.order.nonzero[(top.patient+1):(low.patient)]);
            } 
        else {
            top.patient <- round(length(x) * (trim / 100), digit = 0);
            low.patient <- round(length(x) * (1 - (trim / 100)), digit = 0);
            data.mean <- mean(gene.order, trim = (trim / 100));
            data.sd <- sd(gene.order[(top.patient+1):(low.patient)]);
            }
        (x - data.mean) / data.sd;
        }
    }


# Parallel running
cl <- makeCluster(detectCores()-1);
# register the cluster with the parallel package
registerDoParallel(cl);

# 1. MEAN and SD : method = 'mean', trim = 0
data.mean <- foreach(i=1:nrow(fpkm.tumor.symbol.filter.arg[,patient.part.arg]), .combine = rbind) %dopar% quantify.outliers(fpkm.tumor.symbol.filter.arg[,patient.part.arg][i,]);
data.mean <- data.frame(data.mean);

# 2. TRIMMED MEAN and TRIMMED SD : method = 'mean', trim = 5
data.trimmean <- foreach(i=1:nrow(fpkm.tumor.symbol.filter.arg[,patient.part.arg]), .combine = rbind) %dopar% quantify.outliers(fpkm.tumor.symbol.filter.arg[,patient.part.arg][i,], trim = 5);
data.trimmean <- data.frame(data.trimmean);

# 3. MEDIAN and MAD : method = 'median'
data.median <- foreach(i=1:nrow(fpkm.tumor.symbol.filter.arg[,patient.part.arg]), .combine = rbind) %dopar% quantify.outliers(fpkm.tumor.symbol.filter.arg[,patient.part.arg][i,], methods = 'median');
data.median <- data.frame(data.median);

# 4. KMEAN : method = 'kmean'
data.kmean <- foreach(i=1:nrow(fpkm.tumor.symbol.filter.arg[,patient.part.arg]), .combine = rbind) %dopar% quantify.outliers(fpkm.tumor.symbol.filter.arg[,patient.part.arg][i,], methods = 'kmean')
data.kmean <- data.frame(data.kmean);

stopCluster(cl)





outlier.detection.zrange <- function(x) {
    zrange <- max(x) - min(x);
    zrange.matrix <- c(x, zrange);
    names(zrange.matrix) <- c(names(x), 'zrange');
    zrange.matrix;
    }




# 1. MEAN and SD
data.zrange.mean <- apply(data.mean, 1, outlier.detection.zrange);
data.zrange.mean.t.arg <- data.frame(t(data.zrange.mean));

# 2. TRIMMED MEAN and TRIMMED SD
data.zrange.trimmean <- apply(data.trimmean, 1, outlier.detection.zrange);
data.zrange.trimmean.t.arg <- data.frame(t(data.zrange.trimmean));

# 3. MEDIAN and MAD
data.zrange.median <- apply(data.median, 1, outlier.detection.zrange);
data.zrange.median.t.arg <- data.frame(t(data.zrange.median));




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
data.fraction.kmean.t.arg <- data.frame(t(data.fraction.kmean));




### Final gene-wise matrix #####
### Final gene-wise matrix #####
gene.zrange.fraction.cosine.last.point.bic.arg <- data.frame(cbind(data.zrange.mean.t.arg$zrange,
                     data.zrange.median.t.arg$zrange,
                     data.zrange.trimmean.t.arg$zrange,
                     data.fraction.kmean.t.arg$fraction,
                     data.cosine.bic.t.arg$cosine,
                     data.cosine.bic.t.arg$distribution));
rownames(gene.zrange.fraction.cosine.last.point.bic.arg) <- rownames(fpkm.tumor.symbol.filter);
colnames(gene.zrange.fraction.cosine.last.point.bic.arg) <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean', 'cosine', 'distribution');



zrange.mean <- paste('data.zrange.mean.t.', args.num, sep = '');
assign(zrange.mean, data.zrange.mean.t.arg);

zrange.trimmean <- paste('data.zrange.trimmean.t.', args.num, sep = '');
assign(zrange.trimmean, data.zrange.trimmean.t.arg);

zrange.median <- paste('data.zrange.median.t.', args.num, sep = '');
assign(zrange.median, data.zrange.median.t.arg);

fraction.kmean <- paste('data.fraction.kmean.t.', args.num, sep = '');
assign(fraction.kmean, data.fraction.kmean.t.arg);

cosine.bic <- paste('data.cosine.bic.t.', args.num, sep = '');
assign(cosine.bic, data.cosine.bic.t.arg);

gene.zrange.fraction <- paste('gene.zrange.fraction.cosine.last.point.bic.', args.num, sep = '');
assign(gene.zrange.fraction, gene.zrange.fraction.cosine.last.point.bic.arg);




save(
    fpkm.tumor.symbol.filter,
    patient.part,
    sample.number,
    bic.trim.distribution.fit,
    list = paste0('gene.zrange.fraction.cosine.last.point.bic.', args.num, sep = ''),
    file = paste0('10.Number_Outlier_sample_5method.', args.num, '.short.rda', sep = '')
    );


save(
    fpkm.tumor.symbol.filter,
    patient.part,
    sample.number,
    bic.trim.distribution.fit,
    list = c(
        paste0('data.zrange.mean.t.', args.num, sep = ''),
        paste0('data.zrange.trimmean.t.', args.num, sep = ''),
        paste0('data.zrange.median.t.', args.num, sep = ''),
        paste0('data.fraction.kmean.t.', args.num, sep = ''),
        paste0('data.cosine.bic.t.', args.num, sep = ''),
        paste0('gene.zrange.fraction.cosine.last.point.bic.', args.num, sep = '')),
    file = paste0('10.Number_Outlier_sample_5method.', args.num, '.long.rda', sep = '')
    );




