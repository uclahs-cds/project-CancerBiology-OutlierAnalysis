### explore_genehancer.R ###########################################################
# Description
# get some stats on genehancer db, data is distributed across four files
# genehancer, elements, scores, tissues

### HISTORY ########################################################################
# Version	Date		Developer	Comments
# 0.01		2023-11-02	jlivingstone	initial code
# 0.02		2023-11-08	jlivingstone	add exploratory plots

### PREAMBLE ####################################################################
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities)

setwd('/hot/users/jlivingstone/outlier')

# contains the same information as in GeneHancer_v5.18.bed + element info (n = 409,271)
elements <- read.delim(
	file = '/hot/ref/database/GeneHancer-v5.18/original/GRCh38/GeneHancer_AnnotSV_elements_v5.18.txt',
	as.is = TRUE
	)
# broken down into Enhancer, Promoter, Promoter/Enhancer
# Enhancer  Promoter  Promoter/Enhancer 
# 371839    7370	30062

# collect only the elite GHids to begin with, then merge (n = 173,025)
elite <- elements[which(elements$is_elite == 1), ]

scores <- read.delim(
	file = '/hot/ref/database/GeneHancer-v5.18/original/GRCh38/GeneHancer_AnnotSV_gene_association_scores_v5.18.txt',
	as.is = TRUE
	)
# number of genes per GHid median = 3, mean = 9.5, min = 1, max = 413

scores.elite <- scores[which(scores$is_elite == 1), ]

tissue <- read.delim(
	file = '/hot/ref/database/GeneHancer-v5.18/original/GRCh38/GeneHancer_AnnotSV_tissues_v5.18.txt',
	as.is = TRUE,
	sep = '\t'
	)
# dbSUPER  ENCODE Ensembl FANTOM5   VISTA 
# 365378 1292029 3844417   55407	1903

# number of tissues per GHid median = 9, mean = 15, min = 1, max = 311

#plan of attack, overlap SV breakpoints with GH regions, then see if scores contains outlier gene

data <- data.frame(
	table(elements$is_elite, elements$regulatory_element_type)
	)

colnames(data) <- c('elite', 'Element', 'Frequency')

create.barplot(
	formula = Frequency ~ Element,
	data = data,
	groups = elite,
	col = c('tomato1', 'royalblue'),
	filename = generate.filename('outlier', 'genehancer_element_breakdown', 'png'),
	yaxis.cex = 1,
	xaxis.cex = 1,
	ylab.cex = 1,
	xlab.cex = 1,
	stack = FALSE,
	ylimit = c(0, 250000),
	width = 8,
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
				key = list(
					points = list(
						col = 'black',
						pch = 22,
						cex = 1.25,
						fill = c('tomato1', 'royalblue')
						),
					text = list(
						lab = c('not elite', 'elite')
						),
					padding.text = 1,
					cex = 1.25
					)
				),
			x = 0.65,
			y = 0.95
			)
		),
	resolution = 200
	)

cluster.data <- table(tissue$category, tissue$source)

# for aesthetic reasons
rownames(cluster.data)[5] <- 'induced pluripotent\nstem cell line'
rownames(cluster.data)[4] <- 'in vitro\ndifferentiated cells'

create.heatmap(
	x = cluster.data,
	filename = generate.filename('outlier', 'tissue_source_heatmap', 'png'),
	yaxis.cex = 1,
	xaxis.cex = 1,
	ylab.cex = 1,
	xlab.cex = 1,
	ylab.label = 'Tissue category',
	xlab.label = 'Database source',
	xaxis.lab = colnames(cluster.data),
	xaxis.rot = 45,
	yaxis.lab = rownames(cluster.data),
	scale.data = FALSE,
	same.as.matrix = TRUE,
	print.colour.key = FALSE,
	colour.scheme = c('white', 'darkorchid2'),
	grid.row = FALSE,
	grid.col = FALSE,
	clustering.method = 'none',
	col.pos = rep(1:ncol(cluster.data), each = nrow(cluster.data)),
	row.pos = rep(nrow(cluster.data):1, times = ncol(cluster.data)),
	cell.text = cluster.data,
	text.col = 'black',
	text.cex = 0.60,
	resolution = 200
	)

genes.all <- as.numeric(
	table(scores$symbol)
	)
genes.elite <- as.numeric(
	table(scores.elite$symbol)
	)

unique.ids <- unique(scores$GHid)
genes.per.id.elite <- genes.per.id <- matrix(
	data = NA,
	nrow = length(unique.ids),
	ncol = 1
	)
for (i in 1:length(unique.ids)) {
	if (i %% 1000 == 0) {
		print(i)
		}
	# number of genes associated with that id - is it the same as genes.all
	temp <- scores[which(scores$GHid == unique.ids[i]),]
	genes.per.id[i] <- nrow(temp)
	genes.per.id.elite[i] <- length(which(temp$is_elite == 1))
	}
rownames(genes.per.id.elite) <- rownames(genes.per.id) <- unique.ids

density.data <- list(
	all = genes.all,
	elite = genes.elite
	)

create.densityplot(
	x = density.data,
	filename = generate.filename('outlier', 'number_genes_per_', 'png'),
	col = c('tomato1', 'royalblue'),
	lty = c('solid', 'dashed'),
	yaxis.cex = 1,
	xaxis.cex = 1,
	ylab.cex = 1,
	xlab.cex = 1,
	xlimit = c(0,20),
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
				key = list(
					points = list(
						col = c('tomato1', 'royalblue'),
						pch = 21,
						cex = 1,
						fill = c('tomato1', 'royalblue')
						),
					text = list(
						lab = c("All (n = 264,047)", "Elite only (n = 179,794)")
						),
					padding.text = 1,
					cex = 1
					)
				),
			x = 0.60,
			y = 0.95,
			draw = FALSE
			)
		),
	resolution = 200
	)

