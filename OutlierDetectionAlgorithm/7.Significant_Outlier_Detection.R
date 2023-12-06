#!/usr/bin/env Rscript

# Rscript 7.Significant_Outlier_Detection.R --dataset.name BRCA_EU --working.directory /hot/user/jlivingstone/outlier/run_method \
# --rank.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_BRCA-EU_final_outlier_rank_bic.short.rda \
# --combined.file /hot/users/jlivingstone/outlier/run_method/2023-11-30_Simulated_Data_5method_combine_BRCA_EU.rda \
# --row.chunks 1000 --matrix.chunk 1 --method.iteration 0

### 7.Significant_Outlier_Detection.R ####################################################
# Compute p-values 

# Required R package
library(BoutrosLab.utilities)
library(doParallel);
library(foreach);
library(getopt)
library(parallel);

params <- matrix(
    data = c(
        'dataset.name', 'd', '0', 'character',
        'working.directory', 'w', '0', 'character',
        'rank.file.', 'f', '0', 'character',
        'combined.file', 'c', '0', 'character',
        'row.chunks', 'r', '0', 'character',
	'matrix.chunk', 'm', '0', 'character',
	'method.iteration', 'i', '0', 'character'
        ),
    ncol = 4,
    byrow = TRUE
    );

opt <- getopt(params);
dataset.name <- opt$dataset.name
working.directory <- opt$working.directory
rank.file <- opt$rank.file
combined.file <- opt$combined.file
row.chunks <- as.numeric(opt$row.chunks)
matrix.chunk <- as.numeric(opt$matrix.chunk)

# This will be used identify the number of outlier patients per gene
#   - if '0', use whole patients (first step)
#   - if '1', use n-1 patients (exclude the patient having the largest value)
#   - repeat this '2', '3', '4'... until there is no outlier genes
method.iteration <- opt$method.iteration

# Set the working directory
setwd(working.directory);

# Load the R environment
#   - 1. File from script 1: short version
load(
    file = rank.file
    )

#   - 2. File from script 6
load(
    file = combined.file
    )

### Rank each methods #####
# Function
outlier.rank <- function(outlier.matrix) {
    methods <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean', 'cosine')
    direction <- c(-1, -1, -1, 1, 1)
    rank.matrix <- NULL
    # Give rank for each methods based on z-score range/fraction of kmean
    for (i in 1:length(methods)) {
        rank.methods <- rank(
            x = outlier.matrix[,i] * direction[i],
            ties.method = 'max',
            na.last = 'keep'
            )
        rank.matrix <- cbind(rank.matrix, rank.methods);
        }
    rownames(rank.matrix) <- rownames(outlier.matrix);
    colnames(rank.matrix) <- methods;
    rank.matrix <- data.frame(rank.matrix, stringsAsFactors = FALSE);
    }

### Rank product to determine Top ranked genes #####
# Function
# data.rank: ranked matrix
# NA.number = Number of methods with non-NA should be more than assigned number
outlier.rank.product <- function(data.rank, NA.number = 0) {
    num <- length(which(!is.na(data.rank)));
    if (NA.number >= num) {
        NA;
        }
    else {
        prod(data.rank, na.rm = TRUE) ^ (1 / num);
        }
    }

### Combine matrix
# - relabel the null data
gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relabel <- data.frame(
    gene.zrange.fraction.negative.simulated.sum.bic.5method.1M,
    gene = rownames(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M)
    )
rownames(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relabel) <- paste0('ND', 1:nrow(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M))

# gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- data.frame(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M);

# Assign the row number from start to end
gene.number.start.end.matrix <- data.frame(
    start = numeric(),
    end = numeric(),
    stringsAsFactors = FALSE
    )

number.of.rows <- ceiling(nrow(fpkm.tumor.symbol.filter) / row.chunks)
for (i in 1:number.of.rows) {
    gene.number.start.end.matrix[i,'start'] <- (i - 1) * 1000 + 1
    if (i == number.of.rows) {
        gene.number.start.end.matrix[i,'end'] <- nrow(fpkm.tumor.symbol.filter)
        }
    else {
        gene.number.start.end.matrix[i,'end'] <- i * 1000
        }
    }

cl <- makeCluster(spec = detectCores() - 2);
# register the cluster with the parallel package
registerDoParallel(cl = cl);
clusterExport(
    cl = cl,
    varlist = c('outlier.rank', 'outlier.rank.product')
    )

gene.rank.p.value.one.gene <- NULL;
gene.rank.p.value.one.gene <- foreach (i = gene.number.start.end.matrix[matrix.chunk,'start']:gene.number.start.end.matrix[matrix.chunk,'end'], .combine = rbind) %dopar% {
    methods <- c('zrange.mean', 'zrange.median', 'zrange.trimmean', 'fraction.kmean', 'cosine')
    observed.gene <- gene.zrange.fraction.cosine.last.point.bic[i, methods];
    combine.matrix <- rbind(
        observed.gene,
        gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relabel[, methods]
        )
    # get ranks
    data.rank.bic <- outlier.rank(outlier.matrix = combine.matrix);
    rank.product.bic <- apply(data.rank.bic, 1, outlier.rank.product, NA.number = 3);
    gene.rank.poduct.bic <- data.frame(
        data.rank.bic,
        rank.product.bic
        )
    # gene we are testing is always first
    obs <- rank.product.bic[1]
    # all the rest is simulated data
    null <- rank.product.bic[2:nrow(gene.rank.poduct.bic)]
    length.null <- nrow(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.relabel);
    # permutation test
    obs.p.value <- (sum(obs >= null) + 1) / (length.null + 1)

    p.value.one.gene <- data.frame(gene.rank.poduct.bic[1,], obs.p.value, i = i, gene = rownames(observed.gene));
    p.value.one.gene;
    }

stopCluster(cl = cl);

# to update the object name with the interation number
p.value.one <- paste0('gene.rank.p.value.one.gene.', method.iteration);
assign(p.value.one, gene.rank.p.value.one.gene);

save(
    list = paste0('gene.rank.p.value.one.gene.', method.iteration),
    file = generate.filename('Significant_Outlier_Detection', paste(dataset.name, matrix.chunk, method.iteration, sep = '.'), 'rda')
    );
