#!/usr/bin/env Rscript

### 8.Significant_Outlier_Pvalue_Calculation.R ####################################################
# Compute p-values 


# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of dataset
dataset.name <- 'CCLE';

# Manually enter 'ceiling(nrow(fpkm.tumor.symbol.filter))'
#   - should be changed depending on the dataset
row.chunk.num <- 14;

# args <- commandArgs(trailingOnly = TRUE)

for (i in 1:row.chunk.num) {
    load(file = paste('7.Significant_Outlier_Detection.', dataset.name, '.', i, '.', '0', '.rda', sep = ''));
    p.value.set <- paste('gene.rank.p.value.', i, sep = '');
    assign(p.value.set, get(paste('gene.rank.p.value.one.gene.', '0', sep = '')));    
    }

#1. residue.negative.random.number.bic
gene.p.value.each.null <- NULL;
for (i in 1:row.chunk.num) {
    p.value <- get(paste('gene.rank.p.value.', i, sep = ''));
    gene.p.value.each.null <- rbind(gene.p.value.each.null, p.value);
    }


p.value.all <- paste('gene.rank.p.value.one.gene.p', '0', sep = '');
assign(p.value.all, gene.p.value.each.null);


save(
  list = paste0('gene.rank.p.value.one.gene.p', '0', sep = ''),
  file = paste('8.Significant_Outlier_Pvalue_Calculation.', dataset.name, '.', '0', '.rda', sep = '')
);


