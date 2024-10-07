## compare_num_null.R ################################################################
# Description
# Compare the outlier genes called based on the num.null parameter

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-10-07	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)
library(VennDiagram)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr')

previous <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/run_method',
		'2024-08-28_outlier_formatted_patients_with_ensembl.txt'
		),
	as.is = TRUE
	)
script.genes <- paste(previous$patient, previous$ensembl, sep = ':')
script.genes <- sub('D', 'R', script.genes)

create.venn <- function(script.genes, threshold = 0.001, threshold.type = 'p.value', file.name, file.core) {
	load(
		file = file.path(
			'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/test_outseekr',
			file.name
			)
		)

	number.of.rounds <- length(names(outliers$outlier.test.results.list))
	package.genes <- NULL
	for (i in 1:number.of.rounds) {
		if (threshold.type == 'p.value') {
			signif <- outliers$outlier.test.results.list[[i]][which(outliers$outlier.test.results.list[[i]]$p.value < threshold), ]
		} else {
			signif <- outliers$outlier.test.results.list[[i]][which(outliers$outlier.test.results.list[[i]]$fdr < threshold), ]
			}
	out.genes <- paste(signif$sample, rownames(signif), sep = ':')
	package.genes <- c(package.genes, out.genes)
	}

	gene.list <- list(
		manuscript.genes = script.genes,
		package.genes = package.genes
		)

	venn.diagram(
		x = gene.list,
		filename = generate.filename('outliers', file.core, 'png'),
		fill = c('green', 'purple'),
		category.names = c('Outlier genes', 'Package genes'),
		margin = 0.1,
		resolution = 200,
		imagetype = 'png',
		disable.logging	= TRUE,
		width = 6,
		height = 6,
		units = 'in'
		)
	}

create.venn(
	script.genes = script.genes,
	threshold = 0.001,
	threshold.type = 'p.value',
	file.name = '2024-10-03_pvalue_0.001_remove_one_gene_1000.rda',
	file.core = 'compare_outliers_genes_venn_diagram_pvalue_thousand'
	)

create.venn(
	script.genes = script.genes,
	threshold = 0.001,
	threshold.type = 'p.value',
	file.name = '2024-10-06_pvalue_0.001_remove_one_gene_100000.rda',
	file.core = 'compare_outliers_genes_venn_diagram_pvalue_hundred_thousand'
	)

create.venn(
	script.genes = script.genes,
	threshold = 0.01,
	threshold.type = 'fdr',
	file.name = '2024-10-03_fdr_0.01_remove_one_gene_100000.rda',
	file.core = 'compare_outliers_genes_venn_diagram_fdr_hundred_thousand'
	)

create.venn(
	script.genes = script.genes,
	threshold = 0.01,
	threshold.type = 'fdr',
	file.name = '2024-10-03_fdr_0.01_remove_one_gene_10000.rda',
	file.core = 'compare_outliers_genes_venn_diagram_fdr_ten_thousand'
	)
