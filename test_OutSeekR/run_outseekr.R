## run_outseekr.R ##################################################################
# Description
# script to run the OutSeekR package

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-07-30	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(OutSeekR)
library(devtools)

set.seed(12345)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr')

load_all('/hot/code/jlivingstone/GitHub/uclahs-cds/package-OutSeekR/')

future::plan(future::multisession)

# TCGA Uterine (5245 x 58)
exprs <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr',
		'2024-07-30_outlier_GDC_TCGA-UCS_processed_TPM.tsv'
		),
	as.is = TRUE,
	row.names = 1
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

outliers <- detect.outliers(
	data = brca,
	num.null = 1000,
	p.value.threshold = 0.05,
	fdr.threshold = 0.1
	)

future::plan(future::sequential)
