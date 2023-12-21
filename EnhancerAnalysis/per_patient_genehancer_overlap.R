### per_patient_genehancer_overlap.R ##########################################
# Description
# per patient that has an outlier gene, overlap SV segment with GH db

### HISTORY ###################################################################
# Version	Date		Developer	Comments
# 0.01		2023-12-21	jlivingstone	initial code

### PREAMBLE ##################################################################
library(BoutrosLab.Utilities)

setwd('/hot/user/jlivingstone/outlier/enhancer_analysis')

outliers <- read.delim(
	file = '/hot/user/jlivingstone/outlier/run_method/2023-12-21_Outlier_patients_with_genes_BRCA_EU_cutoff_0.01.txt',
	as.is = TRUE
	)

gh <- read.delim(
	file = '/hot/ref/database/GeneHancer-v5.18/original/GRCh38/GeneHancer_AnnotSV_elements_v5.18.txt',
	as.is = TRUE
	)

# start with amplifications - gain in enhancer or promoter = higher abundance
gains <- read.delim(
	file = '/hot/ref/cohort/ICGC/BRCA/EU/processed/gain_unique_somatic_mutation.BRCA-EU.tsv',
	as.is = TRUE
	)

gene.annot <- read.delim(
	file = '/hot/user/jlivingstone/outlier/NikZainal_2016/original/SupplementaryTable7Transcriptomic342.txt',
	as.is = TRUE
	)
cols.to.keep <- c('UNIQID', 'Ensembl', 'Source', 'Name', 'loc')
gene.annot  <- gene.annot[, match(cols.to.keep, colnames(gene.annot))]


# ugh the WGS names aren't the same as RNA
outliers$sample.name <- sub('R', 'D', outliers$patient)

# n = 104 (only 1/3 of total RNA patients)
sample.overlap <- intersect(outliers$sample.name, gains$submitted_sample_id)

outliers.parsed <- outliers[match(sample.overlap, outliers$sample.name),]

for (i in 1:nrow(outliers.parsed)) {
	genes <- unlist(
		strsplit(
			x = outliers.parsed$outlier_genes[i],
			split = ';'
			)

	temp <- gains[which(gains$submitted_sample_id %in% outliers.parsed$sample.name[i]), ]

	# overlap patient specific gains with gene enhancer element for outlier genes
	
	# need to get elements only for these genes

	}
