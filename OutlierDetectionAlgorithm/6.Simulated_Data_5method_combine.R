#!/usr/bin/env Rscript

### 6.Simulated_Data_5method_combine.R ####################################################
# Combine the 10 chunks of statistics results


# Set the working directory
setwd('RNA-seq/CCLE/four_zero/');

dataset.name <- 'CCLE';


for (i in 1:10) {
    load(file = paste('5.Simulated_Data_5method.', dataset.name, '.', i, '.short.rda', sep = ''));
    gene.zrange <- paste('gene.zrange.fraction.', i, sep = '');
    assign(gene.zrange, gene.zrange.fraction.negative.simulated.sum.1M);    
    gene.zrange.5method <- paste('gene.zrange.fraction.5method.', i, sep = '');
    assign(gene.zrange.5method, gene.zrange.fraction.negative.simulated.sum.bic.5method.1M);  
    }


gene.zrange.fraction.negative.simulated.sum.1M <- NULL;
gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- NULL;
for (i in 1:10) {
    gene.zrange.fraction <- get(paste('gene.zrange.fraction.', i, sep = ''));
    gene.zrange.fraction.negative.simulated.sum.1M <- rbind(gene.zrange.fraction.negative.simulated.sum.1M, gene.zrange.fraction);
    gene.zrange.fraction.5method <- get(paste('gene.zrange.fraction.5method.', i, sep = ''));
   gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- rbind(gene.zrange.fraction.negative.simulated.sum.bic.5method.1M, gene.zrange.fraction.5method);
    }





    save(
    fpkm.tumor.symbol.filter,
    bic.trim.distribution.fit.obs,
    noise.min.off.bic.distribution.fit,
    gene.zrange.fraction.negative.simulated.sum.1M,
    gene.zrange.fraction.negative.simulated.sum.bic.5method.1M,
    file = paste('6.Simulated_Data_5method_combine.', dataset.name, '.rda', sep = '')
    );


