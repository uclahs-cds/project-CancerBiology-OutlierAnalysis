## parse_gdc_exprs.R ##################################################################
# Description
# parse UCS (Uterine) GDC expression data 

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-07-29	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(BoutrosLab.utilities)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr')

data <- read.delim(
	file = file.path(
		'/hot/ref/database/GDC-34.0/TCGA-UCS/harmonized/Transcriptome_Profiling',
		'Gene_Expression_Quantification.tsv'
		),
	as.is = TRUE
	)

tpm <- data[,grep('tpm', colnames(data))]
fpkm <- data[,grep('uq_unstranded', colnames(data))]

# in the same order
tpm.samples <- sub('tpm_unstranded_TCGA.', '', colnames(tpm))
fpkm.sample <- sub('fpkm_uq_unstranded_TCGA.', '', colnames(fpkm))
all(tpm.samples == fpkm.sample)

normalized.cor <- rep(NA, ncol(tpm))
for (i in 1:ncol(tpm)) {
	normalized.cor[i] <- cor(tpm[,i], fpkm[,i], use = 'complete.obs', method = 'pearson')
	}

# use TPM
colnames(tpm) <- tpm.samples
exprs <- data.frame(
	ensembl = data$gene_id,
	log2(tpm + 1)
	)

write.table(
	x = exprs,
	file = generate.filename('outlier', 'GDC_TCGA-UCS_processed_TPM', 'tsv'),
	quote = FALSE,
	sep = '\t',
	row.names = FALSE
	)
