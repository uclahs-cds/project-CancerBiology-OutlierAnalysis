## run_outseekr.R ##################################################################
# Description
# script to run the OutSeekR package

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-07-30	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(OutSeekR)

set.seed(12345)

future::plan(future::multisession)

# TCGA Uterine (5245 x 58)
#exprs <- read.delim(
#	file = file.path(
#		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr',
#		'2024-07-30_outlier_GDC_TCGA-UCS_processed_TPM.tsv'
#		),
#	as.is = TRUE,
#	row.names = 1
#	)

brca <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/NikZainal_2016/processed',
		'2023-11-14_BRCA_EU_transcriptomics_342.tsv'
		),
	as.is = TRUE,
	row.names = 1
	)

gene.var <- apply(brca, 1, function(x) { var(as.numeric(x), na.rm = TRUE) })
exprs <- brca[which(gene.var > 3),]
exprs[is.na(exprs)] <- 0

outliers <- detect.outliers(
	data = exprs,
	num.null = 1000,
	p.value.threshold = 0.05,
	fdr.threshold = 0.1
	)

future::plan(future::sequential)
