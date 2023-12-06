#!/usr/bin/env Rscript

# Rscript Simulated_Data_5method_combine.R --dataset.name BRCA_EU --working.directory /hot/user/jlivingstone/outlier/run_method --file.date 2023-11-28

### 6.Simulated_Data_5method_combine.R ####################################################
# Combine the 10 chunks of statistics results
library(BoutrosLab.utilities)
library(getopt)

params <- matrix(
    data = c(
        'dataset.name', 'd', '0', 'character',
        'working.directory', 'w', '0', 'character',
        'file.date', 'd', '0', 'character'
        ),
    ncol = 4,
    byrow = TRUE
    );

opt <- getopt(params);
dataset.name <- opt$dataset.name
working.directory <- opt$working.directory
file.date <- opt$file.date

# Set the working directory
setwd(working.directory);

# will load and append data to previous
gene.zrange.fraction.negative.simulated.sum.1M.combined <- NULL
gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.combined <- NULL
for (i in 1:10) {
    load(
        file = paste(file.date, '_Simulated_Data_5method_', dataset.name, '.', i, '.short.rda', sep = '')
        )
    gene.zrange.fraction.negative.simulated.sum.1M.combined <- rbind(
        gene.zrange.fraction.negative.simulated.sum.1M.combined,
        gene.zrange.fraction.negative.simulated.sum.1M
        )

    gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.combined <- rbind(
        gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.combined,
        gene.zrange.fraction.negative.simulated.sum.bic.5method.1M
        )
    }
gene.zrange.fraction.negative.simulated.sum.1M <- gene.zrange.fraction.negative.simulated.sum.1M.combined
gene.zrange.fraction.negative.simulated.sum.bic.5method.1M <- gene.zrange.fraction.negative.simulated.sum.bic.5method.1M.combined

save(
    fpkm.tumor.symbol.filter,
    bic.trim.distribution.fit.obs,
    noise.min.off.bic.distribution.fit,
    gene.zrange.fraction.negative.simulated.sum.bic.5method.1M,
    gene.zrange.fraction.negative.simulated.sum.1M,
    file = generate.filename('Simulated_Data_5method_combine', dataset.name, 'rda')
    )
