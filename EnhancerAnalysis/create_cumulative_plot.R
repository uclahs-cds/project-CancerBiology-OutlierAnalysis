## create_cumulative_plot.R ############################################################
# Description
# create cumulative plot to show the permutation test results

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-09-04	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation')

files <- list.files(
	pattern = 'permutation_results_bkpt_padding_1000_overlap.txt'
	)

results <- NULL
for (i in 1:length(files)) {
	data <- read.delim(
		file = file.path(
			'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation',
			files[i]
			),
		as.is = TRUE
		)
	results <- rbind(results, data)
	}
results$n.overlap[is.na(results$n.overlap)] <- 0

ecdf.data <- ecdf(results$n.overlap)

count <- knots(ecdf.data)
ecdf.plot.data <- data.frame(
	overlap.count = count,
	cumulative.probability = ecdf.data(count)
	)

create.scatterplot(
	formula = cumulative.probability ~ overlap.count,
	data = ecdf.plot.data,
	file = generate.filename('outlier', 'enhancer_permutation_cumulative_plot', 'png'),
	xlab.label = expression('Overlap count'),
	ylab.label = expression('Cumulative Probability'),
	col = 'grey10',
	lwd = 3,
	abline.v = 4, 
	abline.lty = 2,
	abline.lwd = 1.5,
	abline.col = 'red2',
	main.cex = 0,
	xaxis.cex = 1,
	yaxis.cex = 1,
	yaxis.tck = c(0.2, 0),
	xaxis.tck = c(0.2, 0),
	xlab.cex = 1.3,
	ylab.cex = 1.5,
	xat = seq(0, 8, 1),
	yat = seq(0, 1, 0.2),
	xaxis.lab = c(0:8),
	xlimits = c(-0.2, 8.5),
	ylimits = c(0, 1.15),
	xaxis.fontface = 1,
	yaxis.fontface = 1,
	add.grid = TRUE,
	xgrid.at = c(2, 4, 6, 8),
	ygrid.at = c(0.2, 0.4, 0.6, 0.8, 1),
	grid.colour = 'grey80',
	type = 'l'
	)
