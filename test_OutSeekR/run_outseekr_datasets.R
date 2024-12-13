## run_outseekr_datasets.R #############################################################
# Description
# script to run the OutSeekR package on the datasets in dataset-breast-cancer

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-12-03	jlivingstone	copied from run_outseekr

### PREAMBLE ###########################################################################
library(BoutrosLab.datasets.breast.cancer)
library(BoutrosLab.utilities)
library(getopt)
library(OutSeekR)

params <- matrix(
	data = c(
		'dataset', 'd', '0', 'character',
		'workers', 'w', '0', 'numeric'
		),
	ncol = 4,
	byrow = TRUE
	)

opt <- getopt(params)
dataset.name <- opt$dataset
workers <- opt$workers

set.seed(12345)
options(future.globals.maxSize = 125000 * 1024^2)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr/datasets_breast_cancer')

future::plan(future::multisession, workers = workers)

# Sjostrom, Hatzis, Cheng, Kao
# values need to be un-logged
print('load dataset')
data <- load.breast.cancer.datasets(
	datasets.to.load = dataset.name,
	with.survival.only = FALSE
	)

logged <- data$all.data[[1]]
unlogged <- 2 ^ logged

# need to adjust extremely small numbers; set to zero
# issues with kmeans clustering; failed Error: empty cluster: try a better set of initial centers

print('detect outliers')
outliers <- detect.outliers(
	data = unlogged,
	num.null = 1000000,
	initial.screen.method = 'fdr',
	p.value.threshold = 0.001,
	fdr.threshold = 0.01,
	kmeans.nstart = 1
	)

save(
	outliers,
	file = generate.filename(dataset.name, 'results_one_million', 'rda')
	)

future::plan(future::sequential)
