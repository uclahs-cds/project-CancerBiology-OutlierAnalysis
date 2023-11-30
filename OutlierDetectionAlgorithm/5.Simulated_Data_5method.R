#!/usr/bin/env Rscript

# Rscript 5.Simulated_Data_5method.R --dataset.name BRCA_EU --working.directory /hot/users/jlivingstone/outlier/run_method \
# --distribution.identification.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_Distribution_Identification_short_BRCA_EU.rda \
# --simulated.data.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_Simulated_data_generation_2_BRCA_EU.1.rda

### 5.Simulated_Data_5method.R ####################################################
# Required R package
library(BoutrosLab.utilities)
library(doParallel);
library(extraDistr);
library(foreach);
library(gamlss);
library(getopt)
library(lsa);
library(truncnorm);
library(parallel);

params <- matrix(
        data = c(
                'dataset.name', 'd', '0', 'character',
                'working.directory', 'w', '0', 'character',
                'distribution.identification.file', 'o', '0', 'character',
                'simulated.data.file', 's', '0', 'character'
                ),
        ncol = 4,
        byrow = TRUE
        );

opt <- getopt(params);
dataset.name <- opt$dataset.name
working.directory <- opt$working.directory
distribution.identification.file <- opt$distribution.identification.file
simulated.data.file <- opt$simulated.data.file

#working.directory <- '/hot/users/jlivingstone/outlier/run_method'
#dataset.name <- 'BRCA_EU'
#distribution.identification.file <- '/hot/users/jlivingstone/outlier/run_method/2023-11-20_Distribution_Identification_short_BRCA_EU.rda'
#simulated.data.file <- '/hot/users/jlivingstone/outlier/run_method/2023-11-21_Simulated_data_generation_2_BRCA_EU.1.rda'

# Set the working directory
setwd(working.directory);

pattern <- "\\d+"
parsed.file <- substr(simulated.data.file, nchar(simulated.data.file) - 5, nchar(simulated.data.file))
index <- gregexpr(pattern = pattern, text = parsed.file)
replicate <- regmatches(parsed.file, index)[[1]]

# load the R environment file saved from 4.Simulated_Data_generation_2.R and 2.Distribution_Identfication.R
load(
	file = simulated.data.file
	)

load(
	file = distribution.identification.file
	)

# sample size
patient.part <- 1:ncol(fpkm.tumor.symbol.filter);
sample.number <- 1:ncol(fpkm.tumor.symbol.filter);
molecular.data.filter <- fpkm.tumor.symbol.filter[, patient.part];
bic.trim.distribution.fit.obs <- bic.trim.distribution.fit;

# Define a minimum value
random.col <- sample(patient.part, 1)
decimal.number.max <- lapply(na.omit(fpkm.tumor.symbol.filter[,random.col]), function(x) {
    decimal.numbers <- sapply(x, function(y) {
        nchar(as.character(y)) - nchar(as.integer(y)) - 1
        })
    return(decimal.numbers)
    })    
add.minimum.value <- 1 / 10 ^ as.numeric(max(unlist(decimal.number.max)));


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
    patient.larger.value <- (length(x) - large.value.number.integer + 1):length(x);
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
        patient.trim.value <- 2:(length(x) - 1);
    } else {
        trim.sample.number <- length(x) * (trim.portion / 100);
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
    add.minimum.value <- 1 / 10 ^ as.numeric(max(unlist(decimal.number.max)));

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
        norm.quantiles <- qtruncnorm(p, a = 0, b = Inf, mean = norm.mean, sd = norm.sd);
        obs.quantile.norm <- quantile(x = sample.fpkm.qq.nozero, prob = p);
        last.cos <- cosine.similarity.large.value.percent(norm.quantiles, obs.quantile.norm, large.value.percent = value.portion);
        }
    else if (2 == distribution.fit) {
        # 2. Log-normal distribution
        mean.log <- mean(sample.fpkm.qq.nozero.trim);
        sd.log <- sd(sample.fpkm.qq.nozero.trim);
        m2 <-  log(mean.log ^ 2 / sqrt(sd.log ^ 2 + mean.log ^ 2));
        sd2 <- sqrt(log(1 + (sd.log ^ 2 / mean.log ^ 2)));
        lnorm.quantile <- qlnorm(p, meanlog = m2, sdlog = sd2);
        obs.quantile.lnorm <- quantile(x = sample.fpkm.qq.nozero, prob = p);
        last.cos <- cosine.similarity.large.value.percent(lnorm.quantile, obs.quantile.lnorm, large.value.percent = value.portion);
        }
    else if (3 == distribution.fit) {
        # 3. Exponential distribution
        exp.rate <- 1 / mean(sample.fpkm.qq.nozero.trim);
        exp.quantile <- qexp(p, rate = exp.rate);
        obs.quantile.exp <- quantile(x = sample.fpkm.qq.nozero, prob = p);
        last.cos <- cosine.similarity.large.value.percent(exp.quantile, obs.quantile.exp, large.value.percent = value.portion);
        }
    else if (4 == distribution.fit) {
        ### 4 gamma distribution
        mean.gamma <- mean(sample.fpkm.qq.nozero.trim);
        sd.gamma <- sd(sample.fpkm.qq.nozero.trim);
        gamma.shape <- (mean.gamma / sd.gamma) ^ 2;
        gamma.rate <- mean.gamma / (sd.gamma ^ 2);
        gamma.quantile <- qgamma(p, shape = gamma.shape, rate = gamma.rate);
        obs.quantile.gamma <- quantile(sample.fpkm.qq.nozero, prob = p);
        last.cos <- cosine.similarity.large.value.percent(gamma.quantile, obs.quantile.gamma, large.value.percent = value.portion);
        }

    cosine.sum.distribution.fit <- c(last.cos, distribution.fit);
    cosine.sum.distribution.fit;
    }

# Determine the distribution
# Find the best fitted distribution

cl <- makeCluster(spec = detectCores() - 2);

# register the cluster with the parallel package
registerDoParallel(cl = cl);
clusterExport(
	cl = cl,
	varlist = 'trim.sample'
	)
clusterEvalQ(
	cl = cl,
	expr = library(gamlss)
	)

# Define a minimum value
random.col <- sample(patient.part, 1)
decimal.number.max <- lapply(na.omit(fpkm.tumor.symbol.filter[,random.col]), function(x) {
    decimal.numbers <- sapply(x, function(y) {
        nchar(as.character(y)) - nchar(as.integer(y)) - 1
        })
    return(decimal.numbers)
    })    

bic.trim.distribution <- NULL;

# Use foreach to iterate over the rows of fpkm.tumor.symbol.filter in parallel
bic.trim.distribution <- foreach (j = 1:nrow(negative.simulated.sum), .combine = rbind) %dopar% {
    sample.fpkm.qq <- round(as.numeric(negative.simulated.sum[j,sample.number]), digits = 6);
    sample.trim.number <- trim.sample(sample.number, 5);
    sample.fpkm.qq.sort <- sort(sample.fpkm.qq)[sample.trim.number];
    sample.fpkm.qq.nozero <- sample.fpkm.qq.sort + add.minimum.value;
    
    
    glm.norm <- gamlss(sample.fpkm.qq.nozero ~ 1, family = NO);
    glm.lnorm <- gamlss(sample.fpkm.qq.nozero ~ 1, family = LNO);
    glm.gamma <- gamlss(sample.fpkm.qq.nozero ~ 1, family = GA);
    glm.exp <- gamlss(sample.fpkm.qq.nozero ~ 1, family = EXP);

    glm.bic <- c(glm.norm$sbc,
                 glm.lnorm$sbc,
                 glm.exp$sbc,
                 glm.gamma$sbc);
    glm.bic;
    }

stopCluster(cl = cl);


# check the distribution again
rownames(bic.trim.distribution) <- rownames(negative.simulated.sum);
bic.trim.distribution.fit <- apply(bic.trim.distribution, 1, which.min);

# Check the cosine similarity
negative.simulated.sum.fit <- cbind(negative.simulated.sum, distribution = bic.trim.distribution.fit);
# run it parallel
cl <- makeCluster(spec = detectCores() - 2);
# register the cluster with the parallel package
registerDoParallel(cl = cl);
clusterExport(
	cl = cl,
	varlist = 'outlier.detection.cosine'
	);
clusterEvalQ(
	cl = cl,
	expr = c(library(lsa), library(SnowballC))
	)

data.cosine.negative <- apply(
	X = negative.simulated.sum.fit,
        MARGIN = 1,
        FUN = outlier.detection.cosine,
        value.portion = 0
	)

stopCluster(cl = cl);


data.cosine.negative.t <- t(data.cosine.negative);
data.cosine.negative.t <- data.frame(data.cosine.negative.t);
colnames(data.cosine.negative.t) <- c('cosine', 'distribution');


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
            data.sd <- sd(gene.order[(top.patient + 1):(low.patient)]);
            }
        result.na <- (x.na - data.mean) / data.sd;
        x[which(!is.na(x))] <- result.na;
        x;
        }
    }

# Parallel running
cl <- makeCluster(detectCores() - 2);
# register the cluster with the parallel package
registerDoParallel(cl = cl);

# 1. MEAN and SD : method = 'mean', trim = 0
data.mean <- foreach (i = 1:nrow(negative.simulated.sum), .combine = rbind) %dopar% quantify.outliers(negative.simulated.sum[i,]);
data.mean <- data.frame(data.mean);

# 2. TRIMMED MEAN and TRIMMED SD : method = 'mean', trim = 5
data.trimmean <- foreach (i = 1:nrow(negative.simulated.sum), .combine = rbind) %dopar% quantify.outliers(negative.simulated.sum[i,], trim = 5);
data.trimmean <- data.frame(data.trimmean);

# 3. MEDIAN and MAD : method = 'median'
data.median <- foreach (i = 1:nrow(negative.simulated.sum), .combine = rbind) %dopar% quantify.outliers(negative.simulated.sum[i,], methods = 'median');
data.median <- data.frame(data.median);

# 4. KMEAN : method = 'kmean'
data.kmean <- foreach (i = 1:nrow(negative.simulated.sum), .combine = rbind) %dopar% quantify.outliers(negative.simulated.sum[i,], methods = 'kmean')
data.kmean <- data.frame(data.kmean);

stopCluster(cl = cl)

outlier.detection.zrange <- function(x) {
    x.na <- na.omit(x)
    zrange <- max(x.na) - min(x.na);
    zrange.matrix <- c(x, zrange);
    names(zrange.matrix) <- c(names(x), 'zrange');
    zrange.matrix;
    }

# 1: MEAN and SD
data.zrange.mean <- apply(data.mean, 1, outlier.detection.zrange);
mean.simulated.negative.1M <- data.frame(t(data.zrange.mean));

# 2. TRIMMED MEAN and TRIMMED SD
data.zrange.trimmean <- apply(data.trimmean, 1, outlier.detection.zrange);
trimmean.simulated.negative.1M <- data.frame(t(data.zrange.trimmean));

# 3. MEDIAN and MAD
data.zrange.median <- apply(data.median, 1, outlier.detection.zrange);
median.simulated.negative.1M <- data.frame(t(data.zrange.median));

### Calculate the kmean fraction #####
# Function
outlier.detection.kmean <- function(x) {
    if (1 == length(unique(as.numeric(x)))) {
        fraction <- NA;
        }
    else {
        cluster.one <- length(x[x == 1]);
        cluster.two <- length(x[x == 2]);
        cluster.sum <- cluster.one + cluster.two;
        smaller.value <- min(cluster.one, cluster.two);
        fraction <- round(smaller.value / cluster.sum, digit = 4);
        }
    fraction.matrix <- c(x, fraction);
    names(fraction.matrix) <- c(names(x), 'fraction');
    fraction.matrix;
    }

# 4. KMEAN fraction
data.fraction.kmean <- apply(data.kmean, 1, outlier.detection.kmean);
kmean.simulated.negative.1M <- data.frame(t(data.fraction.kmean));

### Final gene-wise matrix #####
gene.zrange.fraction.negative.simulated.sum.1M <- cbind(
	mean.simulated.negative.1M$zrange,
        median.simulated.negative.1M$zrange,
        trimmean.simulated.negative.1M$zrange,
        kmean.simulated.negative.1M$fraction
	)
rownames(gene.zrange.fraction.negative.simulated.sum.1M) <- rownames(negative.simulated.sum);
colnames(gene.zrange.fraction.negative.simulated.sum.1M) <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean');


# Final statistic matrix
gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- cbind(
	gene.zrange.fraction.negative.simulated.sum.1M[,c(1, 2, 3, 4)],
	data.cosine.negative.t$cosine,
	data.cosine.negative.t$distribution
	)
colnames(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M) <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean', 'cosine', 'distribution');

save(
    fpkm.tumor.symbol.filter,
    sample.number,
    bic.trim.distribution.fit.obs,
    bic.trim.distribution.fit,
    noise.min.off.bic.distribution.fit,
    gene.zrange.fraction.negative.simulated.sum.1M,
    data.cosine.negative.t,
    gene.zrange.fraction.negative.simulated.sum.bic.5method.1M,
    file = generate.filename('Simulated_Data_5method', paste(dataset.name, replicate, 'short', sep = '.'), 'rda')
    );

save(
    fpkm.tumor.symbol.filter,
    patient.part,
    bic.trim.distribution.fit.obs,
    bic.trim.distribution.fit,
    negative.random.number.bic,
    obs.residue.quantile.trim,
    residue.negative.random.number.bic,
    noise.min.off.bic.distribution.fit,
    noise.min.off.bic.distribution.residue,
    negative.simulated.sum,
    data.mean,
    data.trimmean,
    data.median,
    data.kmean,
    gene.zrange.fraction.negative.simulated.sum.1M,
    data.cosine.negative.t,
    gene.zrange.fraction.negative.simulated.sum.bic.5method.1M,
    file = generate.filename('Simulated_Data_5method', paste(dataset.name, replicate, 'long', sep = '.'), 'rda')
    )
