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
outlier.genes$ensemb <- matrix$Ensembl[match(outlier.genes$outlier_genes, matrix$Name)]
	
# load package results
load(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr',
		'2024-08-22_BRCA_EU_results_using_package_function.rda'
		)
	)

package.genes <- names(which(outliers$num.outliers.adjusted > 1))

# when FDR is set to 0.01
# 2/3 genes overlap
overlap.genes <- intersect(package.genes, outlier.genes$ensemb)
# n = 5803 !
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
