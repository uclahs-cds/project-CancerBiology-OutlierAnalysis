### 2.Distribution_Identification.R ####################################################
# Identify the distribution for generating simulated data


# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of dataset
dataset.name <- 'CCLE';

# load the R environment file saved from 1.Outlier_detection_5method.R
# load(file = '2023-09-05_CCLE_final_outlier_rank_bic.long.rda')
load(file = '2023-09-18_CCLE_four_zero_final_outlier_rank_bic.long.rda')

# Required R packages
install.packages('gamlss', repo = 'http://cran.us.r-project.org');
library(gamlss);
install.packages('doParallel', repo = 'http://cran.us.r-project.org');
library(doParallel);
install.packages('foreach', repo = 'http://cran.us.r-project.org');
library(foreach);
install.packages('parallel', repo = 'http://cran.us.r-project.org');
library(parallel);


### Function ###
# Define a minimum value
least.significant.digit <- function(x) {
    decimal.number.max <- sapply(
        X = na.omit(x),
        FUN = function(y) {
            decimal.numbers <- sapply(
                X = y,
                FUN = function(z) {
                    nchar(as.character(z)) - nchar(as.integer(z)) - 1
                    }
                );
            decimal.numbers;
            }
        );
    1 / 10^as.numeric(max(decimal.number.max));
    }

calculate.residuals.observed.data <- function(x, distr) {
    # Define a minimum value
    add.minimum.value <- least.significant.digit(x);
    # Add the minimum value so that no values in the data are zero.
    # Do the same for a 5% trimmed sample of the data.  The trimmed
    # data will be used to calculate the parameters of the
    # distributions, while the untrimmed data will be used to generate
    # residuals.
    x.nozero <- x + add.minimum.value;
    x.trim <- trim.sample(x, 0.05);
    x.nozero.trim <- x.trim + add.minimum.value;

    # Quantile
    p <- ppoints(x.nozero);
    # Calculate residuals between the quantiles of the observed data
    # and the quantiles of the optimal distribution for the observed
    # data.
    if (1 == distr) {
        # 1. Normal distribution
        norm.mean <- mean(x.nozero.trim);
        norm.sd <- sd(x.nozero.trim);
        norm.quantiles <- qnorm(p, mean = norm.mean, sd = norm.sd);
        obs.quantile.norm <- quantile(x.nozero, prob = p);
        obs.residual.non.trim <- obs.quantile.norm - norm.quantiles;
        } else if (2 == distr) {
        # 2. Log-normal distribution
        mean.log <- mean(x.nozero.trim);
        sd.log <- sd(x.nozero.trim);
        m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
        sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
        lnorm.quantile <- qlnorm(p, meanlog = m2, sdlog = sd2);
        obs.quantile.lnorm <- quantile(x.nozero, prob = p);
        obs.residual.non.trim <- obs.quantile.lnorm - lnorm.quantile;
        } else if (3 == distr) {
        # 3. Exponential distribution
        exp.rate <- 1 / mean(x.nozero.trim);
        exp.quantile <- qexp(p, rate = exp.rate);
        obs.quantile.exp <- quantile(x.nozero, prob = p);
        obs.residual.non.trim <- obs.quantile.exp - exp.quantile;
        } else if (4 == distr) {
        ### 4. Gamma distribution
        mean.gamma <- mean(x.nozero.trim);
        sd.gamma <- sd(x.nozero.trim);
        gamma.shape <- (mean.gamma/sd.gamma)^2;
        gamma.rate <- mean.gamma/(sd.gamma^2);
        gamma.quantile <- qgamma(p, shape = gamma.shape, rate = gamma.rate);
        obs.quantile.gamma <- quantile(x.nozero, prob = p);
        obs.residual.non.trim <- obs.quantile.gamma - gamma.quantile;
        }

    obs.residual.non.trim;
    }

# Identifying the distribution of each residual
# -1. Get residuals
obs.residual.quantile <- lapply(
    X = seq_len(nrow(fpkm.tumor.symbol.filter)),
    FUN = function(i) {
        calculate.residuals.observed.data(
            x = fpkm.tumor.symbol.filter[i, ],
            distr = bic.trim.distribution.fit[i]
            );
        }
    );
obs.residual.quantile <- do.call(
    what = rbind,
    args = obs.residual.quantile
    );
obs.residual.quantile <- data.frame(obs.residual.quantile);
rownames(obs.residual.quantile) <- rownames(fpkm.tumor.symbol.filter);

# Trim each 5%
obs.residue.quantile.trim <- apply(
    X = obs.residue.quantile,
    MARGIN = 1,
    FUN = sort
    );
obs.residue.quantile.trim <- data.frame(t(obs.residue.quantile.trim));
obs.residue.quantile.trim <- obs.residue.quantile.trim[
    ,
    trim.sample(seq_len(ncol(obs.residue.quantile.trim)), trim = 0.05)
    ];
obs.residue.quantile.trim.max <- apply(obs.residue.quantile.trim, 1, max);

# -2. Get distribution of residues
# Using "min off" the values
cl <- makeCluster(detectCores()-1) # create a cluster with all available cores except one
registerDoParallel(cl) # register the cluster with the parallel package

clusterExport(cl, "trim.sample");
clusterEvalQ(cl, library(gamlss));

noise.min.off.bic.distribution <- NULL

noise.min.off.bic.distribution <- foreach(j = 1:nrow(obs.residue.quantile.trim), .combine = rbind) %dopar% {
    sample.fpkm.qq <- round(as.numeric(obs.residue.quantile.trim[j,]), digits = 6);
    sample.fpkm.qq.sort <- sort(sample.fpkm.qq);
    if (min(sample.fpkm.qq.sort) < 0) {
        sample.fpkm.qq.nozero <- sample.fpkm.qq.sort - min(sample.fpkm.qq.sort) + add.minimum.value;
    } else {
        sample.fpkm.qq.nozero <- sample.fpkm.qq.sort + add.minimum.value ;
    }

    glm.norm <- gamlss(sample.fpkm.qq.nozero ~ 1, family=NO);
    glm.lnorm <- gamlss(sample.fpkm.qq.nozero ~ 1, family=LNO);
    glm.gamma <- gamlss(sample.fpkm.qq.nozero ~ 1, family=GA);
    glm.exp <- gamlss(sample.fpkm.qq.nozero ~ 1, family=EXP);

    glm.bic <- c(glm.norm$sbc,
                 glm.lnorm$sbc,
                 glm.exp$sbc,
                 glm.gamma$sbc)
    glm.bic;

    }

stopCluster(cl) # stop the cluster after the computation is done

# Find the best fitted distribution
# - BIC
rownames(noise.min.off.bic.distribution) <- rownames(obs.residue.quantile.trim);
noise.min.off.bic.distribution.fit <- apply(noise.min.off.bic.distribution, 1, which.min);


# save the R environment
#   - short version
save(
    fpkm.tumor.symbol.filter,
    patient.part,
    sample.number,
    bic.trim.distribution.fit,
    gene.zrange.fraction.cosine.last.point.bic,
    gene.rank.order.5method.cosine.last.point.bic,
    obs.residue.quantile.trim,
    noise.min.off.bic.distribution.fit,
    # file = '2023-09-06_CCLE_final_outlier_rank_bic_distribution.short.rda'
    file = paste('2.Distribution_Identification.', dataset.name, '.short.rda', sep = '')
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
    obs.residue.quantile.trim,
    noise.min.off.bic.distribution.fit,
    #file = '2023-09-06_CCLE_final_outlier_rank_bic_distribution.long.rda'
    file = paste('2.Distribution_Identification.', dataset.name, '.long.rda', sep = '')
    );
