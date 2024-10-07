## run_outseekr.R ##################################################################
# Description
# script to run the OutSeekR package

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-07-30	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(OutSeekR)
#library(devtools)

set.seed(12345)
options(future.globals.maxSize = 16 * 1024^3)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr')

#load_all('/hot/code/jlivingstone/GitHub/uclahs-cds/package-OutSeekR/')

future::plan(future::multisession, workers = 8)

# TCGA Uterine (17247 x 57) - filtered for genes detected in 50% of samples
exprs <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr',
		'2024-08-22_outlier_GDC_TCGA-UCS_processed_TPM.tsv'
		),
	as.is = TRUE
	)

ucs.outliers <- detect.outliers(
	data = exprs,
	num.null = 1000,
	initial.screen.method = 'fdr',
#	p.value.threshold = 0.05,
	fdr.threshold = 0.01,
	kmeans.nstart = 1
	)

save(
	ucs.outliers,
	file = genereate.filename('outliers', 'GDC_TCGA-UCS_results_using_package_function', 'rda')
	)

# values are non-logged
brca <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/NikZainal_2016/processed',
		'2024-08-14_BRCA_EU_processed_unlogged.tsv'
		),
	as.is = TRUE
	)
#dim(brca)
#[1] 17696   342

# need to adjust extremely small numbers; set to zero
# issues with kmeans clustering; failed Error: empty cluster: try a better set of initial centers
brca[brca > 0 & brca < 1 * 10^(-50)] <- 0

# remove one gene that has the same value in 225 samples (weird)
eu <- brca[-which(rownames(brca) == 'ENSG00000244428'),]

outliers <- detect.outliers(
	data = eu,
	num.null = 1000000,
#	initial.screen.method = 'p.value',
	initial.screen.method = 'fdr',
	p.value.threshold = 0.001,
	fdr.threshold = 0.01,
	kmeans.nstart = 1
	)

save(
	outliers,
	file = genereate.filename('outliers', 'BRCA_EU_results_one_million', 'rda')
	)

future::plan(future::sequential)
