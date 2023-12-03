#!/usr/bin/env Rscript

### 11.Simulated_Data_5method_combine.R ####################################################
# Combine the 10 chunks of statistics results


# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of dataset
dataset.name <- 'CCLE';

# Combine all 10 chunks
args <- commandArgs(trailingOnly = TRUE)


for (i in 1:10) {
    load(paste('9.Number_Outlier_sample_Simulated_Data_5method.', dataset.name, '.', i, '.', args, '.short.rda', sep = ''));
    gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- get(paste('gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.', i, sep = ''));
    all.statistics <- paste('gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.', args, i, sep = '');
    assign(all.statistics, gene.zrange.fraction.negative.simulated.sum.bic.5method.1M);    
    }

#1. residue.negative.random.number.bic
gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- NULL;
for (i in 1:10) {
    p.value <- get(paste('gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.', args, i, sep = ''));
    gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- rbind(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M, p.value);
    }


gene.zrange.fraction <- paste('gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.', args, sep = '');
assign(gene.zrange.fraction, gene.zrange.fraction.negative.simulated.sum.bic.5method.1M);

save(
    fpkm.tumor.symbol.filter,
    bic.trim.distribution.fit,
    list = paste0('gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.', args, sep = ''),
    file = paste('11.Number_Outlier_sample_Simulated_Data_5method_combine.', dataset.name, '.', args, '.rda', sep = '')
    );
