### check_random_pairs.R ######################################################
# Description
# read in simulated sets and check how often a sample and gene are paired

### HISTORY ###################################################################
# Version	Date		Developer	Comments
# 0.01		2024-08-21	jlivingstone	initial code

### PREAMBLE ##################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

type <- 'not_elite'
if (type == 'elite') {
	setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation/one_thousand_permutations/elite')
	placeholder <- 'elite_'
} else {
	setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation/one_thousand_permutations/not_elite')
	placeholder <- '_'
	}
	
files <- list.files(
	pattern = 'rda'
	)
files <- files[grep('seeds', files)]
start <- nchar('2024-XX-XX_outlier_seeds_')
seed <- substr(files, start +1, start + 7)
ord <- order(
	as.numeric(
		unlist(
			lapply(seed, function(x) { strsplit(x, '_')[[1]][1] })
			)
		)
	)
ordered.files <- files[ord]

final.list <- list()
for (i in 1:length(ordered.files)) {
	print(ordered.files[i])
	# loads simulated.set object
	load(ordered.files[i])
	ind <- which(!sapply(X = simulated.set, FUN = is.null))

	final.list <- c(final.list, simulated.set[ind])
	}

save(
	final.list,
	file = generate.filename('outlier', paste0(placeholder, 'simulated_sets'), 'rda')
	)

pairings <- unlist(
	sapply(
		X = final.list,
		FUN = function(x) { paste(x$sample.name, x$outlier_genes, sep = ':') }
		)
	)

count <- as.numeric(table(pairings))
print(table(count))
#     1      2      3      4 
#170989   3905     63      3

create.histogram(
	x = count,
	filename = generate.filename('outlier', paste0(placeholder, 'gene_sample_pairing_count'), 'png'),
	type = 'count',
	ylimits = c(0, 200000),
	yat = seq(0, 200000, 50000),
	breaks = 16,
	ylab.cex = 1,
	ylab.label = 'Count',
	xlab.cex = 1,
	xlab.label = '',
	xaxis.cex = 1,
	yaxis.cex = 1
	)
