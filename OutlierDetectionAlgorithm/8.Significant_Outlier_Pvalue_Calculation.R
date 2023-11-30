#!/usr/bin/env Rscript

# Rscript 8.Significant_Outlier_Pvalue_Calculation.R --dataset.name BRCA_EU --working.directory /hot/user/jlivingstone/outlier/run_method --row.chunk 18
### 8.Significant_Outlier_Pvalue_Calculation.R ####################################################
# Compute p-values 
library(BoutrosLab.utilities)
library(getopt)

params <- matrix(
    data = c(
	'dataset.name', 'd', '0', 'character',
        'working.directory', 'w', '0', 'character',
	'row.chunk', 'r', '0', 'character'
        ),
    ncol = 4,
    byrow = TRUE
    );

opt <- getopt(params);
dataset.name <- opt$dataset.name
working.directory <- opt$working.directory
row.chunk.num <- opt$row.chunk

# Set the working directory
setwd(working.directory)


for (i in 1:row.chunk.num) {
    load(
	file = paste('Significant_Outlier_Detection', dataset.name, i, '0', 'rda', sep = '.')
	)
    p.value.set <- paste('gene.rank.p.value', i, sep = '.');
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
  file = generate.filename('Significant_Outlier_Pvalue_Calculation', paste(dataset.name, '0', sep = '.'), 'rda')
  )


