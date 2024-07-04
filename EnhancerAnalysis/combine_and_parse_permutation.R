### combine_and_parse_permutation #############################################################
# Description

### HISTORY ###################################################################################
# Version	Date		Developer	Comments
# 0.01		2024-03-20	jlivingstone	initial code
# 0.02		2024-03-25	jlivingstone	add density plots

### PREAMBLE #################################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation')

seed.files <- list.files(pattern = 'outlier_seeds')
if (is.elite == TRUE) {
	files <- seed.files[grep('.*elite.*txt', seed.files, perl = TRUE)]
	placeholder <- 'elite_'
} else {
	files <- seed.files[grep('.*results_bkpt.*txt', seed.files, perl = TRUE)]
	placeholder <- ''
	}

results <- data.frame()
for (i in 1:length(files)) {
	temp <- read.delim(
		files[i],
		as.is = TRUE
		)
	results <- rbind(results, temp)
	}

if (is.elite == TRUE) {
	load(
		file = file.path(
			'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation',
			'2024-03-18_true_dataset_bkpt_plus_1000_overlap_elite.rda'
			)
		)
} else {
	load(
		file = file.path(
			'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation',
			'2024-03-26_true_dataset_bkpt_plus_1000_overlap.rda'
			)
		)
	}
overlap <- which(unlist(lapply(elements.overlap, function(x) { nrow(x) })) > 0)
overlap.n <- length(overlap)

results$n.overlap[is.na(results$n.overlap)] <- 0
p.value <- length(which(results$n.overlap >= overlap.n)) / nrow(results)

create.histogram(
	x = results$n.overlap,,
	filename = generate.filename('outlier', paste0('overlap_', placeholder, 'histogram_plot'), 'png'),
	type = 'count',
:	abline.v = overlap.n,
	abline.lty = 'dashed',
	abline.lwd = 2,
	ylab.cex = 1,
	ylab.label = 'Count',
	xlab.cex = 1,
	xaxis.cex = 1,
	yaxis.cex = 1,
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
				key = list(
					text = list(
						lab = as.expression(substitute(P == paste(base %*% 10^exponent),
						list(
							base = scientific.notation(p.value, digits = 2, type = 'list')$base,
							exponent = scientific.notation(p.value, type = 'list')$exponent)
							)
						),
					cex = 1,
					fontface = 'bold'
					)
				)
			),
			x = 0.68,
			y = 0.98
			)
		)
	)

# density plot for distance analysis

final.outliers <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation',
		'2024-03-25_outliers_final_outliers_overlapped_gh_exprs.tsv'
		),
	as.is = TRUE
	)

load(
	file = file.path(
w		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/distance_analysis',
		'2024-03-19_outlier_distance_of_mutation_to_outlier_gene_gh_element.rda'
		)
	)

outliers <- final.outliers[match(names(distance.results), final.outliers$outlier_genes),]

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/distance_analysis/histogram')

for (i in 1:length(distance.results)) {

	gene <- names(distance.results)[i]
	print(gene)

	ind.to.remove <- match(outliers$sample.name[i], rownames(distance.results[[gene]]))
	null <- log10(distance.results[[gene]]$distance[-ind.to.remove] + 1)
	true <- log10(distance.results[[gene]]$distance[ind.to.remove] + 1)

	x.max <- ceiling(max(c(null, true)) + 1)
	x.min <- floor(min(c(null, true)) - 1)

	create.histogram(
		x = null,
		filename = generate.filename('outlier', paste0(gene, '_distance_histogram_plot'), 'png'),
		type = 'count',
		abline.v = true,
		abline.lty = 'dashed',
		abline.lwd = 2,
		xlimits = c(x.min, x.max),
		breaks = 16,
		main = gene,
		main.cex = 1,
		ylab.cex = 1,
		ylab.label = 'Count',
		xlab.cex = 1,
		xlab.label = expression(bold('-log')[bold('10')] * bold('distance)')),
		xaxis.cex = 1,
		yaxis.cex = 1
		)
	}
