### 1.Outlier_Detection.R ####################################################

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




# Load the RNA abundance file
#   - any input data is available
# fpkm.tumor.symbol: gene x sample matrix
# example:
#   fpkm.tumor.symbol <- read.csv(file = 'METADOR/GSE167977_20210223_matadorCounts.csv', check.names = FALSE, stringsAsFactors = F, sep = ',', row.names = 1);


# The original value is log2(FPKM+1)
#   - make it non log format
fpkm.tumor.symbol <- 2^fpkm.tumor.symbol-1;

# Number of samples
#    - if the last column has symbol, it should be 1:(ncol(fpkm.tumor.symbol))-1)
patient.part <- 1:ncol(fpkm.tumor.symbol);
sample.number <- 1:ncol(fpkm.tumor.symbol);


# Get the genes with less than 1% of zero values
#   - excludes genes which have zero values more than 99%
zero.portion <- apply(fpkm.tumor.symbol[,patient.part], 1, function(x) {length(x[0 == x]) / length(patient.part)});
fpkm.tumor.symbol.filter <- fpkm.tumor.symbol[rownames(fpkm.tumor.symbol) %in% names(zero.portion[0.01 > zero.portion]),];
molecular.data.filter <- fpkm.tumor.symbol.filter[, patient.part];




### Outlier detection function
# Default : methods = 'mean', trim = 0
# 1. MEAN and SD : methods = 'mean', trim = 0
# 2. TRIMMED MEAN and TRIMMED SD : methods = 'mean', trim = 5
# 3. MEDIAN and MAD : methods = 'median'
# 4. KMEAN : methods = 'kmean'
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
# cl <- makeCluster(detectCores()-1);
cl <- 2;
# register the cluster with the parallel package
registerDoParallel(cl);

# 1. MEAN and SD : method = 'mean', trim = 0
data.mean <- foreach(i=1:nrow(molecular.data.filter), .combine = rbind) %dopar% quantify.outliers(molecular.data.filter[i,]);
data.mean <- data.frame(data.mean);

# 2. TRIMMED MEAN and TRIMMED SD : method = 'mean', trim = 5
data.trimmean <- foreach(i=1:nrow(molecular.data.filter), .combine = rbind) %dopar% quantify.outliers(molecular.data.filter[i,], trim = 5);
data.trimmean <- data.frame(data.trimmean);

# 3. MEDIAN and MAD : method = 'median'
data.median <- foreach(i=1:nrow(molecular.data.filter), .combine = rbind) %dopar% quantify.outliers(molecular.data.filter[i,], methods = 'median');
data.median <- data.frame(data.median);

# 4. KMEAN : method = 'kmean'
data.kmean <- foreach(i=1:nrow(molecular.data.filter), .combine = rbind) %dopar% quantify.outliers(molecular.data.filter[i,], methods = 'kmean')
data.kmean <- data.frame(data.kmean);

stopCluster(cl)




### Calculate the range of z-score #####
# Function
outlier.detection.zrange <- function(x) {
  x.na <- na.omit(x)
  zrange <- max(x.na) - min(x.na);
  zrange.matrix <- c(x, zrange);
  names(zrange.matrix) <- c(names(x), 'zrange');
  zrange.matrix;
    }

# 1. MEAN and SD
data.zrange.mean <- apply(data.mean, 1, outlier.detection.zrange);
data.zrange.mean.t <- data.frame(t(data.zrange.mean));

# 2. TRIMMED MEAN and TRIMMED SD
data.zrange.trimmean <- apply(data.trimmean, 1, outlier.detection.zrange);
data.zrange.trimmean.t <- data.frame(t(data.zrange.trimmean));

# 3. MEDIAN and MAD
data.zrange.median <- apply(data.median, 1, outlier.detection.zrange);
data.zrange.median.t <- data.frame(t(data.zrange.median));




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
data.fraction.kmean.t <- data.frame(t(data.fraction.kmean));





# 5. Cosine similarity
# function: Compute the cosine similarity of the largest data point
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


# Trimming function
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

# Determine the distribution
# Find the best fitted distribution
cl <- makeCluster(2);
# register the cluster with the parallel package
registerDoParallel(cl);
clusterExport(cl, "trim.sample");
clusterEvalQ(cl, library(gamlss));

# Define a minimum value
random.col <- sample(patient.part, 1)
decimal.number.max <- lapply(na.omit(fpkm.tumor.symbol.filter[,random.col]), function(x) {
    decimal.numbers <- sapply(x, function(y) {
        nchar(as.character(y)) - nchar(as.integer(y)) - 1
        })
    return(decimal.numbers)
    })    
add.minimum.value <- 1 / 10^as.numeric(max(unlist(decimal.number.max)));


bic.trim.distribution <- NULL;

# Use foreach to iterate over the rows of fpkm.tumor.symbol.filter in parallel
bic.trim.distribution <- foreach(j = 1:nrow(fpkm.tumor.symbol.filter), .combine = rbind) %dopar% {
    sample.fpkm.qq <- round(as.numeric(fpkm.tumor.symbol.filter[j,patient.part]), digits = 6);
    sample.trim.number <- trim.sample(sample.number, 5);
    sample.fpkm.qq.sort <- sort(sample.fpkm.qq)[sample.trim.number];
    sample.fpkm.qq.nozero <- sample.fpkm.qq.sort + add.minimum.value;
    
    
    glm.norm <- gamlss(sample.fpkm.qq.nozero ~ 1, family=NO);
    glm.lnorm <- gamlss(sample.fpkm.qq.nozero ~ 1, family=LNO);
    glm.gamma <- gamlss(sample.fpkm.qq.nozero ~ 1, family=GA);
    glm.exp <- gamlss(sample.fpkm.qq.nozero ~ 1, family=EXP);

    glm.bic <- c(glm.norm$sbc,
                 glm.lnorm$sbc,
                 glm.exp$sbc,
                 glm.gamma$sbc);
    glm.bic;
    }

stopImplicitCluster();


# Find the best fitted distribution
# - BIC
rownames(bic.trim.distribution) <- rownames(fpkm.tumor.symbol.filter);
bic.trim.distribution.fit <- apply(bic.trim.distribution, 1, which.min);

# Check the cosine similarity
fpkm.tumor.symbol.filter.bic.fit <- cbind(fpkm.tumor.symbol.filter, distribution = bic.trim.distribution.fit);

# run it parallel
cl <- makeCluster(2);
# register the cluster with the parallel package
registerDoParallel(cl);
clusterExport(cl, "outlier.detection.cosine");
clusterEvalQ(cl, c(library(lsa), library(SnowballC)));

data.cosine.bic <- apply(fpkm.tumor.symbol.filter.bic.fit, 
                         1, 
                         outlier.detection.cosine, 
                         value.portion = 0);

stopImplicitCluster();

data.cosine.bic.t <- data.frame(t(data.cosine.bic));
colnames(data.cosine.bic.t) <- c('cosine', 'distribution');







### Final gene-wise matrix #####
gene.zrange.fraction.cosine.last.point.bic <- data.frame(cbind(data.zrange.mean.t$zrange,
                     data.zrange.median.t$zrange,
                     data.zrange.trimmean.t$zrange,
                     data.fraction.kmean.t$fraction,
                     data.cosine.bic.t$cosine,
                     data.cosine.bic.t$distribution));
rownames(gene.zrange.fraction.cosine.last.point.bic) <- rownames(fpkm.tumor.symbol.filter);
colnames(gene.zrange.fraction.cosine.last.point.bic) <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean', 'cosine', 'distribution');






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



# Compute the rank product
data.rank.bic <- outlier.rank(gene.zrange.fraction.cosine.last.point.bic);

rank.product.bic <- apply(data.rank.bic, 1, outlier.rank.product, NA.number = 3);
gene.rank.poduct.bic <- cbind(data.rank.bic, 
                              rank.product.bic,
                              distribution = gene.zrange.fraction.cosine.last.point.bic$distribution);
gene.rank.order.5method.cosine.last.point.bic <- gene.rank.poduct.bic[order(gene.rank.poduct.bic$rank.product, decreasing = FALSE),];


# save the R environment
#   - short version
save(
    fpkm.tumor.symbol.filter,
    patient.part,
    sample.number,
    bic.trim.distribution.fit,
    gene.zrange.fraction.cosine.last.point.bic,
    gene.rank.order.5method.cosine.last.point.bic,
    file = generate.filename(dataset.name, 'final_outlier_rank_bic.short', 'rda')
    );

#   - long version
save(
    fpkm.tumor.symbol.filter,
    patient.part,
    sample.number,
    bic.trim.distribution.fit,
    data.zrange.mean.t,
    data.zrange.median.t,
    data.zrange.trimmean.t,
    data.fraction.kmean.t,
    data.cosine.bic.t,
    gene.zrange.fraction.cosine.last.point.bic,
    gene.rank.order.5method.cosine.last.point.bic,
    file = generate.filename(dataset.name, 'final_outlier_rank_bic.long', 'rda')
    );


