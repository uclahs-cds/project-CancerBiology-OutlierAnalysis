## compare_null_data.R ################################################################
# Description
# Compare the null distributions created by JYH scripts and package

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-09-03	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/run_method')

null.data <- NULL
gene.ids <- NULL

for (i in 1:10) {
	print(i)
	load(
		file = paste('2023-11-22_Simulated_data_generation_2_BRCA_EU.', i, '.rda', sep = '')
		)

	gene.ids <- rbind(
		gene.ids,
		rownames(negative.simulated.sum)
		)

	null.data <- rbind(
		null.data,
		negative.simulated.sum
		)
	}

save(
	null.data,
	file = generate.filename('outlier', 'scripts_null_data_negative_simulated_sum', 'rda')
	)

# how often was each gene selected to create null from ?
summary(as.numeric(table(rownames(null.data))))
#   Min. 1st Qu.  Median	Mean 3rd Qu.	Max. 
#  32.00   51.00   56.00   56.49   62.00   87.00 

load('2023-11-30_Simulated_Data_5method_combine_BRCA_EU.rda')

data <- data.frame(
	gene.zrange.fraction.negative.simulated.sum.bic.5method.1M,
	stringsAsFactors = FALSE
	)

toplot <- data.frame(
	mean = data$zrange.mean,
	median = data$zrange.median,
	trimmean = data$zrange.trimmean
	)

create.densityplot(
	x = toplot,
	file = generate.filename('outlier', 'null_data_density_plot', 'png'),
	col = default.colours(3),
	xlimits = c(0, 50),
	yaxis.cex = 1,
	xaxis.cex = 1,
	ylab.cex = 1,
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
				key = list(
					points = list(
						col = default.colours(3),
						pch = 21,
						cex = 1,
						fill = default.colours(3)
						),
					text = list(
						lab = names(toplot)
						),
					cex = 1
					)
				),
			x = 0.65,
			y = 0.97
			)
		)
	 )
