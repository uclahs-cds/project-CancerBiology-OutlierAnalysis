#!/usr/bin/env Rscript

### 13.Number_Outlier_sample_Significant_Outlier_Pvalue_Calculation.R ####################################################
# Compute p-values 


# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

# Set the name of dataset
dataset.name <- 'CCLE';

# Number of excluded patients
args <- commandArgs(trailingOnly = TRUE)


# Manually enter 'ceiling(nrow(fpkm.tumor.symbol.filter))'
#   - should be changed depending on the dataset
row.chunk.num <- 14;


for (i in 1:row.chunk.num) {
    load(file = paste('12.Number_Outlier_sample_Significant_Outlier_Detection.', dataset.name, '.', i, '.', args, '.rda', sep = ''));
    p.value.set <- paste('gene.rank.p.value.', i, sep = '');
    assign(p.value.set, get(paste('gene.rank.p.value.one.gene.', args, sep = '')));    
    }

#1. residue.negative.random.number.bic
gene.p.value.each.null <- NULL;
for (i in 1:row.chunk.num) {
    p.value <- get(paste('gene.rank.p.value.', i, sep = ''));
    gene.p.value.each.null <- rbind(gene.p.value.each.null, p.value);
    }


p.value.all <- paste('gene.rank.p.value.one.gene.p', args, sep = '');
assign(p.value.all, gene.p.value.each.null);


save(
  list = paste0('gene.rank.p.value.one.gene.p', args, sep = ''),
  file = paste('13.Number_Outlier_sample_Significant_Outlier_Pvalue_Calculation.', dataset.name, '.', args, '.rda', sep = '')
);


