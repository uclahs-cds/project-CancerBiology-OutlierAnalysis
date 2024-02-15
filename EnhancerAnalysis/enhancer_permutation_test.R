## enhancer_permutation_test.R #########################################################
# Description
# perform a permutation test to see if outlier genes are more likely to have mutations
# in enhancer regions

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-02-15	jlivingstone	initial code

### PREAMBLE ###########################################################################
setwd('/hot/user/jlivingstone/outlier/enhancer_analysis')

# pick random 'non' outlier gene and overlap sample mutations with gene enhancer regions for that gene

# how many genes to pick ? see how many true outlier genes exist in genehancer
gh <- read.delim(
	file = file.path('/hot/ref/database/GeneHancer-v5.18/processed/GRCh38', 'GeneHancer_AnnotSV_elements_v5.18.txt'),
	as.is = TRUE
	)

outliers <- read.delim(
	file = file.path('/hot/user/jlivingstone/outlier/run_method/', '2023-12-21_Outlier_patients_with_genes_BRCA_EU_cutoff_0.01.txt'),
	as.is = TRUE
	)
outlier.genes <- unlist(strsplit(outliers$outlier_genes, ';'))
gh.gene <- unique(gh$symbol)

n.outlier.genes <- length(intersect(gh.gene, outlier.genes))

data <- read.delim2(
	file = '/hot/users/jlivingstone/outlier/NikZainal_2016/original/SupplementaryTable7Transcriptomic342.txt',
	header = TRUE,
	row.names = 1
	)
dataset.genes <- unique(data$Name)
overlap.genes <- intersect(dataset.genes, gh.gene)

# which samples overlap mutations and rna
muts <- read.delim(
	file = file.path(
		'/hot/ref/cohort/ICGC/BRCA/EU/processed/wgs/',
		'deletion.BRCA-EU.tsv'
		),
	as.is = TRUE
	)
muts$submitted_sample_id <- sub('a2', 'a', muts$submitted_sample_id)
muts$submitted_sample_id <- sub('a3', 'a', muts$submitted_sample_id)

muts.sample <- unique(muts$submitted_sample_id)
samples <- intersect(muts.sample, outliers$sample.name)

# create random sample and gene combination
# https://www.biorxiv.org/content/10.1101/2024.02.05.579021v1.full.pdf
seeds <- c(51404, 366306, 423647, 838004, 50135, 628019, 97782, 253505, 659767, 13142)

for (i in 1:length(seeds)) {
	set.seed(seeds[i])

	toprint <- data.frame(
		patient = character(),
		outlier_genes = character(),
		stringsAsFactors = FALSE
		)
	picked.samples <- sample(samples, n.outlier.genes, replace = TRUE)
	picked.genes <- sample(overlap.genes, n.outlier.genes, replace = TRUE)

	s <- unique(picked.samples)
	for (j in 1:length(s)) {
		toprint[j, 'patient'] <- s[j]
		toprint[j, 'outlier_genes'] <- paste0(picked.genes[grep(s[j], picked.samples)], collapse = ';')
		}

	write.table(
		x = toprint,
		file = generate.filename(seeds[i], 'random_outlier_genes', 'tsv'),
		quote = FALSE,
		row.names = FALSE,
		sep = '\t'
		)
	}
