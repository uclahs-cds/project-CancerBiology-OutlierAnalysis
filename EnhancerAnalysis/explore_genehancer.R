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

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone')

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

ids.per.genes.all <- data.frame(
	table(scores$symbol)
	)
ids.per.genes.elite <- data.frame(
	table(scores.elite$symbol)
	)

genes.per.id <- data.frame(
	table(scores$GHid)
	)
genes.per.id.elite <- data.frame(
	table(scores.elite$GHid)
	)

density.data <- list(
	all = ids.per.genes.all$Freq,
	elite = ids.per.genes.elite$Freq
	)

per.gene <- create.densityplot(
	x = density.data,
	# filename = generate.filename('outlier', 'GHid_per_gene_densityplot', 'png'),
	col = c('tomato1', 'royalblue'),
	lty = c('solid', 'dashed'),
	yaxis.cex = 1,
	xaxis.cex = 1,
	ylab.cex = 1,
	xlab.cex = 1,
	xlab.label = 'GHids per gene',
	xlimit = c(0, 20),
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
						lab = c('All genes (n = 264,047)', 'Elite only (n = 179,794)')
						),
					padding.text = 1,
					cex = 1
					)
				),
			x = 0.50,
			y = 0.95,
			draw = FALSE
			)
		),
	resolution = 200
	)

genes.per.id <- data.frame(
	table(scores$GHid)
	)
genes.per.id.elite <- data.frame(
	table(scores.elite$GHid)
	)

id.density.data <- list(
	all = genes.per.id$Freq,
	elite = genes.per.id.elite$Freq
	)

per.id <- create.densityplot(
	x = id.density.data,
	# filename = generate.filename('outlier', 'genes_per_GHid_densityplot', 'png'),
	col = c('tomato1', 'royalblue'),
	lty = c('solid', 'dashed'),
	yaxis.cex = 1,
	xaxis.cex = 1,
	ylab.cex = 1,
	xlab.cex = 1,
	xlab.label = 'Genes per GHid',
	xlimit = c(0, 20),
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
						lab = c('All GHids (n = 409,271)', 'Elite only (n = 241,676)')
						),
					padding.text = 1,
					cex = 1
					)
				),
			x = 0.50,
			y = 0.95,
			draw = FALSE
			)
		),
	resolution = 200
	)

create.multipanelplot(
	plot.objects = list(per.gene, per.id),
	filename = generate.filename('outlier', 'GeneHancer_densityplots', 'png'),
	resolution = 200,
	layout.width = 1,
	layout.height = 2
	)

ids <- unique(tissue$GHid)
ghid.source <- ghid.source.len <- matrix(
        data = NA,
        nrow = length(ids),
        ncol = 1
        )
for (i in 1:length(ids)) {
        if (i %% 1000 == 0) {
                print(i)
                }
        ghid.source[i] <- paste(unique(tissue$source[which(tissue$GHid == ids[i])]), collapse = ';')
        ghid.source.len[i] <- length(unique(tissue$source[which(tissue$GHid == ids[i])]))
        }
tissue.parsed <- data.frame(
	ghid = ids,
	source = ghid.source,
	number_evidence_sources = ghid.source.len,
	stringsAsFactors = FALSE
	)

data <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone',
		'NikZainal_2016/original/SupplementaryTable7Transcriptomic342.txt'
		),
	as.is = TRUE
	)

overlap.genes <- intersect(data$Name, gh$symbol)
gh <- elite[which(elite$symbol %in% overlap.genes),]

multiple <- names(which(table(gh$symbol) > 1))
num.chr <- rep(NA, length(multiple))
for (i in 1:length(multiple)) {
	num.chr[i] <- length(unique(gh$chr[which(gh$symbol == multiple[i])]))
	}
# all elements for the same gene are on the same chromosome
