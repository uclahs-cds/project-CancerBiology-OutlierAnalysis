### 4.Simulated_Data_generation_2.R ####################################################
# Generate the simulated data: part2

# Rscript 4.Simulated_Data_generation_2.R --dataset.name BRCA_EU --working.directory /hot/users/jlivingstone/outlier/run_method \
# --distribution.identification.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_Distribution_Identification_short_BRCA_EU.rda \
# --simulated.data.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_Simulated_data_generation_1_BRCA_EU.1.rda

# Required R package
library(BoutrosLab.utilities)
library(doParallel);
library(extraDistr);
library(foreach);
library(getopt)
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
#simulated.data.file <- '/hot/users/jlivingstone/outlier/run_method/2023-11-21_Simulated_data_generation_1_BRCA_EU.1.rda'

# replicate number is parsed from the input file
pattern <- '\\d+'
parsed.file <- substr(simulated.data.file, nchar(simulated.data.file) - 5, nchar(simulated.data.file))
index <- gregexpr(pattern = pattern, text = parsed.file)
replicate <- regmatches(parsed.file, index)[[1]]

# Set the working directory
setwd(working.directory);

# load the R environment file saved from 2.Distribution_Identification.R and 3.Simulated_Data_generation_1.R
load(
    file = distribution.identification.file
    )

load(
    file = simulated.data.file
    )

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
add.minimum.value <- 1 / 10 ^ as.numeric(max(unlist(decimal.number.max)));

# 2. residual ** using magic numbers
residue.negative.random.number.bic <- obs.residue.quantile.trim[match(substr(rownames(negative.random.number.bic), 1, 15), substr(rownames(obs.residue.quantile.trim), 1, 15)),];

noise.min.off.bic.distribution.residue <- noise.min.off.bic.distribution.fit[match(substr(rownames(negative.random.number.bic), 1, 15), substr(names(noise.min.off.bic.distribution.fit), 1, 15))];

# create a cluster of CPU cores to use in parallel processing
cl <- makeCluster(spec = detectCores() - 2)
# register the cluster with the parallel package
registerDoParallel(cl = cl);
clusterExport(
    cl = cl,
    varlist = 'add.minimum.value'
    )
clusterEvalQ(
    cl = cl,
    expr = library(extraDistr)
    )

random.col <- sample(patient.part, 1)
decimal.number.max <- lapply(na.omit(fpkm.tumor.symbol.filter[,random.col]), function(x) {
    decimal.numbers <- sapply(x, function(y) {
        nchar(as.character(y)) - nchar(as.integer(y)) - 1
        })
    return(decimal.numbers)
    })
add.minimum.value <- 1 / 10 ^ as.numeric(max(unlist(decimal.number.max)));

negative.random.number.noise.bic <- foreach(i = 1:nrow(residue.negative.random.number.bic), .combine = 'rbind') %dopar% {

    sample.fpkm.qq <- as.numeric(residue.negative.random.number.bic[i,]);
    sample.fpkm.qq.sort <- sort(sample.fpkm.qq);
    if (min(sample.fpkm.qq.sort) < 0) {
        sample.fpkm.qq.nozero <- sample.fpkm.qq.sort - min(sample.fpkm.qq.sort) + add.minimum.value;
        }
    else {
        sample.fpkm.qq.nozero <- sample.fpkm.qq.sort + add.minimum.value;
        }

    if (1 == noise.min.off.bic.distribution.residue[i]) {
            # 1. Normal distribution
            norm.mean <- mean(sample.fpkm.qq.nozero);
            norm.sd <- sd(sample.fpkm.qq.nozero);
            simulated.sample <- rtnorm(length(sample.number), mean = norm.mean, sd = norm.sd, a = 0);
            }
    else if (2 == noise.min.off.bic.distribution.residue[i]) {
            # 2. Log-normal distribution
            mean.log <- mean(sample.fpkm.qq.nozero);
            sd.log <- sd(sample.fpkm.qq.nozero);
            m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
            sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
            simulated.sample <- rlnorm(n = length(sample.number), m2, sd2);
            }
    else if (3 == noise.min.off.bic.distribution.residue[i]) {
            # 3. Exponential distribution
            exp.rate <- 1 / mean(sample.fpkm.qq.nozero);
            simulated.sample <- rexp(n = length(sample.number), rate = exp.rate);
            }
    else if (4 == noise.min.off.bic.distribution.residue[i]) {
           # 4. Gamma distribution
           mean.gamma <- mean(sample.fpkm.qq.nozero);
           sd.gamma <- sd(sample.fpkm.qq.nozero);
           gamma.shape <- (mean.gamma / sd.gamma) ^ 2;
           gamma.rate <- mean.gamma / (sd.gamma ^ 2)
           simulated.sample <- rgamma(n = length(sample.number), gamma.shape, gamma.rate);
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

stopCluster(cl = cl)

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
    file = generate.filename('Simulated_Data_generation_2', paste(dataset.name, replicate, sep = '.'), 'rda')
    )
