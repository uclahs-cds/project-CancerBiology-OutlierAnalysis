### combine_and_parse_permutation #############################################################
# Description

### HISTORY ###################################################################################
# Version	Date		Developer	Comments
# 0.01		2024-03-20	jlivingstone	initial code

### PREAMBLE #################################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

setwd('/hot/user/jlivingstone/outlier/enhancer_analysis/permutation')

seed.files <- list.files(pattern = 'outlier_seeds')
files <- seed.files[grep('txt', seed.files)]

results <- data.frame()
for (i in 1:length(files)) {
	temp <- read.delim(
		files[i],
		as.is = TRUE
		)
	results <- rbind(results, temp)
	}

load('2024-03-18_true_dataset_bkpt_plus_1000_overlap_elite.rda')
overlap <- which(unlist(lapply(elements.overlap, function(x) { nrow(x) })) == 1)
overlap.n <- length(overlap)

results$n.overlap[is.na(results$n.overlap)] <- 0
p.value <- length(which(results$n.overlap >= overlap.n)) / nrow(results)

# density plot for distance analysis

final.outliers <- read.delim(
	file = file.path(
		'/hot/user/jlivingstone/outlier/enhancer_analysis/permutation',
		'2024-03-25_outliers_final_outliers_overlapped_gh_exprs.tsv'
		),
	as.is = TRUE
	)

load(
	file = file.path(
		'/hot/user/jlivingstone/outlier/enhancer_analysis/distance_analysis',
		'2024-03-19_outlier_distance_of_mutation_to_outlier_gene_gh_element.rda'
		)
	)

outliers <- final.outliers[match(names(distance.results), final.outliers$outlier_genes),]

for (i in 1:length(distance.results)) {

	gene <- names(distance.results)[i]
	print(gene)

	ind.to.remove <- match(outliers$sample.name[i], rownames(distance.results[[gene]]))
	null <- log10(distance.results[[gene]]$distance[-ind.to.remove]) 
	true <- log10(distance.results[[gene]]$distance[ind.to.remove])

	create.histogram(
		x = null,
		filename = generate.filename('outlier', paste0(gene, '_distance_histogram_plot'), 'png'),
		type = 'count',
		abline.v = true,
		abline.lty = 'dashed',
		abline.lwd = 2,
		breaks = 16,
		main = gene,
		main.cex = 1,
		ylab.cex = 1,
		ylab.label = 'Count',
		xlab.cex = 1,
		xlab.label = expression(paste('-log'[10], '(distance)')),
		xaxis.cex = 1,
		yaxis.cex = 1
		)
	}
