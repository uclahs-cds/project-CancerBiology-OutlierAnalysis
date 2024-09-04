## compare_outlier_genes_called.R ######################################################
# Description
# Compare the outliers that are called by function vs JYH scripts

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-08-27	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

matrix <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/NikZainal_2016/original',
		'SupplementaryTable7Transcriptomic342.txt'
		),
	as.is = TRUE
	)

unlog <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/NikZainal_2016/processed',
		'2024-08-14_BRCA_EU_processed_unlogged.tsv'
		),
	as.is = TRUE
	)

# load script results
outlier.genes <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/run_method',
		'2024-03-15_outlier_formatted_patients_with_genes.txt'
		),
	as.is = TRUE
	)
outlier.genes$ensembl <- matrix$Ensembl[match(outlier.genes$outlier_genes, matrix$Name)]
	
# load package results
load(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr',
		'2024-08-22_BRCA_EU_results_using_package_function.rda'
		)
	)

package.genes <- names(which(outliers$num.outliers.adjusted > 0))

# when FDR is set to 0.01
# 2/3 genes overlap
overlap.genes <- intersect(package.genes, outlier.genes$ensemb)
# n = 8923 !
# 226 / 238 

# RP11-255N4.2
plot.outlier.barplot <- function(id = 'ENSG00000241313') {
	outliers.samples <- names(which(outliers$fdr[grep(id, rownames(outliers$fdr)), ] < 0.1))

	unlog.filtered <- unlog[grep(id, rownames(unlog)),]
	#matrix.filtered <- matrix[grep('ENSG00000241313', matrix$Ensembl), grep('PR', colnames(matrix))]
	samples <- colnames(unlog.filtered)

	ord <- order(as.numeric(unlog.filtered))

	toplot <- data.frame(
		order = length(unlog.filtered):1,
		value = as.numeric(unlog.filtered)[ord],
		samples = samples[ord]
		)

	bar.colour <- rep('black', nrow(toplot))
	bar.colour[match(outliers.samples, toplot$samples)] <- 'red'

	create.barplot(
		formula = value ~ order,
		data = toplot,
		file = generate.filename('outlier', paste0(id, '_outlier_barplot'), 'png'),
		width = 12,
		xaxis.tck = 0,
		xaxis.cex = 0.5,
		xaxis.lab = toplot$samples,
		xaxis.rot = 45,
		xlimit = c(0, nrow(toplot)),
		xat = match(outliers.samples, toplot$samples),
		col = bar.colour,
		xlab.label = '',
		ylab.label = '',
		resolution = 200
		)
	}

# outliers due to very low abundance ...
plot.outlier.barplot(id = 'ENSG00000241313')
plot.outlier.barplot(id = 'ENSG00000064692')

# number of outliers based on FDR cutoff
cutoff <- c(0.1, 0.05, 0.01)

outlier.samples.per.transcript <- outlier.numbers <- matrix(data = NA, ncol = length(cutoff), length(package.genes))

for (j in 1:length(package.genes)){
	id <- package.genes[j]
	for (i in 1:length(cutoff)) {
		outlier.samples.per.transcript[j,i] <- paste(
			names(which(outliers$fdr[grep(id, rownames(outliers$fdr)), ] < cutoff[i])),
			collapse = ';'
			)
		outlier.numbers[j,i] <- length(which(outliers$fdr[grep(id, rownames(outliers$fdr)), ] < cutoff[i]))
		}
	}
colnames(outlier.samples.per.transcript) <- colnames(outlier.numbers) <- cutoff
rownames(outlier.samples.per.transcript) <- rownames(outlier.numbers) <- package.genes

# length(which(outlier.numbers[,2] > 0))
# [1] 6084 - sum 17459
# length(which(outlier.numbers[,3] > 0))
# [1] 15 - sum = 33

write.table(
	x = outlier.samples.per.transcript,
	file = generate.filename('outlier', 'package_outlier_samples_per_gene_per_cutoff', 'tsv'),
	quote = FALSE,
	sep = '\t'
	)

write.table(
	x = outlier.numbers,
	file = generate.filename('outlier', 'package_outlier_numbers_per_cutoff', 'tsv'),
	quote = FALSE,
	sep = '\t'
	)

fdr.genes <- names(outlier.numbers[which(outlier.numbers[,3] > 0), 3])
# intersect(fdr.genes, outlier.genes$ensemb)
# 9/15

# get abundance value for outlier genes (n = 6084)
get.transcript.sample.pairs <- function(index = 2, fdr.value = '0.05', cutoff = 5) {
	per.transcript <- outlier.samples.per.transcript[, index][which(outlier.samples.per.transcript[, index] != '')]

	pairs <- NULL
	for (i in 1:length(per.transcript)) {
		if (i %% 1000 == 0) { print(i) }
		tmp <- unlog[names(per.transcript[i]), unlist(strsplit(per.transcript[i], ';')), drop = FALSE]
		toadd <- data.frame(
			transcript = rownames(tmp),
			sample = names(tmp),
			values = as.numeric(tmp),
			stringsAsFactors = FALSE
			)
		pairs <- rbind(pairs, toadd)
		}

	filtered.pairs <- pairs[-which(pairs$values < cutoff), ]
	write.table(
		x = filtered.pairs,
		file = generate.filename('outlier', paste0('transcript_sample_pairs_fdr_', fdr.value, '_tpm_filtered'), 'tsv'),
		row.names = FALSE,
		sep = '\t'
		)
	return(filtered.pairs)
	}

pairs.fdr.pointone <- get.transcript.sample.pairs(
	index = 1,
	fdr.value = '0.1',
	cutoff = 5
	)

pairs.fdr.pointzerofive <- get.transcript.sample.pairs(
	index = 2,
	fdr.value = '0.05',
	cutoff = 5
	)

pairs.fdr.pointzeroone <- get.transcript.sample.pairs(
	index = 3,
	fdr.value = '0.01',
	cutoff = 5
	)
