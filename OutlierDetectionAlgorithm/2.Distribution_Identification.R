### 2.Distribution_Identification.R ####################################################
# Identify the distribution for generating simulated data

# Rscript 2.Distribution_Identification.R --dataset.name BRCA_EU --working.directory /hot/users/jlivingstone/outlier/run_method --outlier.rank.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_BRCA-EU_final_outlier_rank_bic.long.rda

# Required R packages
library(BoutrosLab.utilities)
library(doParallel);
library(foreach);
library(gamlss);
library(getopt);
library(parallel);

params <- matrix(
        data = c(
                'dataset.name', 'd', '0', 'character',
                'working.directory', 'w', '0', 'character',
                'outlier.rank.file', 'f', '0', 'character'
                ),
        ncol = 4,
        byrow = TRUE
        );

opt <- getopt(params);
dataset.name <- opt$dataset.name
working.directory <- opt$working.directory
outlier.rank.file <- opt$outlier.rank.file

working.directory <- '/hot/users/jlivingstone/outlier/run_method'
dataset.name <- 'BRCA_EU'
outlier.rank.file <- '/hot/users/jlivingstone/outlier/run_method/2023-11-20_BRCA-EU_final_outlier_rank_bic.long.rda'

# Set the working directory
setwd(working.directory);

# load the R environment file saved from 1.Outlier_detection_5method.R
load(
    file = outlier.rank.file
    )

### Function ###
# Define a minimum value
print('define a minimum value')
random.col <- sample(patient.part, 1)
decimal.number.max <- lapply(na.omit(fpkm.tumor.symbol.filter[,random.col]), function(x) {
    decimal.numbers <- sapply(x, function(y) {
        nchar(as.character(y)) - nchar(as.integer(y)) - 1
        })
    return(decimal.numbers)
    })
add.minimum.value <- 1 / 10 ^ as.numeric(max(unlist(decimal.number.max)));
    
# function: Trim 5% of samples from each side
trim.sample <- function(x, trim.portion = 5) {
    trim.sample.number <- length(x) * (trim.portion / 100);
    trim.sample.number.integer <- round(trim.sample.number, digits = 0);
    patient.trimr.value <- (trim.sample.number.integer + 1):(length(x) - trim.sample.number.integer);
    patient.trimr.value;
    }


# Identifying the distribution of each residue
# -1. Get residues
print('Get residues')

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

obs.residue.quantile <- foreach(i = 1:nrow(fpkm.tumor.symbol.filter), .combine = rbind) %dopar% {

    sample.fpkm.qq <- round(as.numeric(fpkm.tumor.symbol.filter[i,patient.part]), digits = 6);
    sample.fpkm.qq.sort <- sort(sample.fpkm.qq);
    sample.fpkm.qq.nozero <- sample.fpkm.qq.sort + add.minimum.value;

    sample.trim.number <- trim.sample(seq(sample.fpkm.qq.sort), 5);
    sample.fpkm.qq.trim <- sort(sample.fpkm.qq)[sample.trim.number];
    sample.fpkm.qq.nozero.trim <- sample.fpkm.qq.trim + add.minimum.value;

    # Quantile
    # p <- seq(0.001, 0.patient.number, 0.001);
    p <- ppoints(length(patient.part));

    if (1 == bic.trim.distribution.fit[i]) {
        # 1. Normal distribution
        norm.mean <- mean(sample.fpkm.qq.nozero.trim);
        norm.sd <- sd(sample.fpkm.qq.nozero.trim);
        norm.quantiles <- qnorm(p, mean = norm.mean, sd = norm.sd);
        obs.quantile.norm <- quantile(sample.fpkm.qq.nozero, prob = p);
        obs.residue.non.trim <- obs.quantile.norm - norm.quantiles;
        }

    else if (2 == bic.trim.distribution.fit[i]) {
        # 2. Log-normal distribution
        mean.log <- mean(sample.fpkm.qq.nozero.trim);
        sd.log <- sd(sample.fpkm.qq.nozero.trim);
        m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
        sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
        lnorm.quantile <- qlnorm(p, meanlog = m2, sdlog = sd2);
        obs.quantile.lnorm <- quantile(sample.fpkm.qq.nozero, prob = p);
        obs.residue.non.trim <- obs.quantile.lnorm - lnorm.quantile;
        }

    else if (3 == bic.trim.distribution.fit[i]) {
        # 3. Exponential distribution
        exp.rate <- 1 / mean(sample.fpkm.qq.nozero.trim);
        exp.quantile <- qexp(p, rate = exp.rate);
        obs.quantile.exp <- quantile(sample.fpkm.qq.nozero, prob = p);
        obs.residue.non.trim <- obs.quantile.exp - exp.quantile;
        }

    else if (4 == bic.trim.distribution.fit[i]) {
        ### 4 gamma distribution
        mean.gamma <- mean(sample.fpkm.qq.nozero.trim);
        sd.gamma <- sd(sample.fpkm.qq.nozero.trim);
        gamma.shape <- (mean.gamma / sd.gamma) ^ 2;
        gamma.rate <- mean.gamma / (sd.gamma ^ 2);
        gamma.quantile <- qgamma(p, shape = gamma.shape, rate = gamma.rate);
        obs.quantile.gamma <- quantile(sample.fpkm.qq.nozero, prob = p);
        obs.residue.non.trim <- obs.quantile.gamma - gamma.quantile;
        }

    obs.residue.non.trim
    }

stopCluster(cl = cl)

obs.residue.quantile <- data.frame(obs.residue.quantile);
rownames(obs.residue.quantile) <- rownames(fpkm.tumor.symbol.filter);
# Trim each 5%
obs.residue.quantile.trim <- apply(
    X = obs.residue.quantile,
    MARGIN = 1,
    FUN - function(x) {
        sort(as.numeric(x))
        }
    );
obs.residue.quantile.trim <- data.frame(t(obs.residue.quantile.trim));

sample.trim.number <- trim.sample(patient.part, 5);
obs.residue.quantile.trim <- obs.residue.quantile.trim[,sample.trim.number];
obs.residue.quantile.trim.max <- apply(obs.residue.quantile.trim, 1, max);


# -2. Get distribution of residues
print('get distribution of resides')
# Using "min off" the values
cl <- makeCluster(spec = detectCores() - 2)
registerDoParallel(cl = cl)

clusterExport(
    cl = cl,
    varlist = 'trim.sample'
    )
clusterEvalQ(
    cl = cl,
    expr = library(gamlss)
    )

noise.min.off.bic.distribution <- NULL

noise.min.off.bic.distribution <- foreach(j = 1:nrow(obs.residue.quantile.trim), .combine = rbind) %dopar% {
    sample.fpkm.qq <- round(as.numeric(obs.residue.quantile.trim[j,]), digits = 6);
    sample.fpkm.qq.sort <- sort(sample.fpkm.qq);
    if (min(sample.fpkm.qq.sort) < 0) {
        sample.fpkm.qq.nozero <- sample.fpkm.qq.sort - min(sample.fpkm.qq.sort) + add.minimum.value;
    } else {
        sample.fpkm.qq.nozero <- sample.fpkm.qq.sort + add.minimum.value ;
    }

    glm.norm <- gamlss(sample.fpkm.qq.nozero ~ 1, family = NO);
    glm.lnorm <- gamlss(sample.fpkm.qq.nozero ~ 1, family = LNO);
    glm.gamma <- gamlss(sample.fpkm.qq.nozero ~ 1, family = GA);
    glm.exp <- gamlss(sample.fpkm.qq.nozero ~ 1, family = EXP);

    glm.bic <- c(glm.norm$sbc,
                 glm.lnorm$sbc,
                 glm.exp$sbc,
                 glm.gamma$sbc)
    glm.bic;
    }

stopCluster(cl = cl)

# Find the best fitted distribution
# - BIC
rownames(noise.min.off.bic.distribution) <- rownames(obs.residue.quantile.trim);
noise.min.off.bic.distribution.fit <- apply(noise.min.off.bic.distribution, 1, which.min);

print('save objects')
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
    file = generate.filename('Distribution_Identification_short', dataset.name, 'rda')
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
    file = generate.filename('Distribution_Identification_long', dataset.name, 'rda')
    );
