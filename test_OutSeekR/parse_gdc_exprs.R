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

tpm <- data[grep('ENSG', data$gene_id), grep('tpm', colnames(data))]
fpkm <- data[grep('ENSG', data$gene_id), grep('uq_unstranded', colnames(data))]

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

# since this is a test and we want to reduce the number of genes
# lets apply a variation cutoff > 1
gene.var <- apply(log2(tpm + 1), 1, var, na.rm = TRUE)
tokeep <- which(gene.var > 1)

exprs <- data.frame(
	ensembl = data$gene_id[tokeep],
	log2(tpm + 1)[tokeep,]
	)

write.table(
	x = exprs,
	file = generate.filename('outlier', 'GDC_TCGA-UCS_processed_TPM', 'tsv'),
	quote = FALSE,
	sep = '\t',
	row.names = FALSE
	)
