#!/usr/bin/env Rscript

# Rscript 8.Significant_Outlier_Pvalue_Calculation.R --dataset.name BRCA_EU \
# --working.directory /hot/user/jlivingstone/outlier/run_method \
# --method.iteration 0

### 8.Significant_Outlier_Pvalue_Calculation.R ####################################################
# Compute p-values 
library(BoutrosLab.utilities)
library(getopt)

params <- matrix(
    data = c(
	'dataset.name', 'd', '0', 'character',
        'working.directory', 'w', '0', 'character',
	'method.iteration', 'i', '0', 'character'
        ),
    ncol = 4,
    byrow = TRUE
    );

opt <- getopt(params);
dataset.name <- opt$dataset.name
working.directory <- opt$working.directory
method.iteration <- opt$method.iteration

# Set the working directory
setwd(working.directory)

files <- list.files(
	pattern = 'Significant_Outlier_Detection'
	)

p.value.all <- NULL
for (i in 1:length(files)) {
    load(
	file = files[i]
	)
    assign(
	x = 'variable.name',
	value = paste('gene.rank.p.value.one.gene', method.iteration, sep = '.')
	)
    p.value.all <- rbind(
	p.value.all,
	get(x = variable.name)
	)
    }
p.value.all <- p.value.all[order(p.value.all$i),]
p.value.all$q.value <- p.adjust(
	p = p.value.all$obs.p.value,
	method = 'fdr'
	)

# assign back to original variable name
assign(x = variable.name, value = p.value.all)

save(
    list = paste('gene.rank.p.value.one.gene', method.iteration, sep = '.'),
    file = generate.filename('Significant_Outlier_Pvalue_Calculation', paste(dataset.name, method.iteration, sep = '.'), 'rda')
    )
