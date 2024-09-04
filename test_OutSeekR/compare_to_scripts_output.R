## compare_to_scripts_output.R #########################################################
# Description
# Compare the outliers that are called by function vs JYH scripts

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-08-30	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(BoutrosLab.utilities)

load(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/run_method',
		'2023-11-20_BRCA-EU_final_outlier_rank_bic.short.rda'
		)
	)

brca <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/NikZainal_2016/processed',	
		'2024-08-14_BRCA_EU_processed_unlogged.tsv'
		),
	as.is = TRUE
	)

all(rownames(annot.filter) == rownames(gene.zrange.fraction.cosine.last.point.bic))
gene.zrange.fraction.cosine.last.point.bic$ensembl <- annot.filter$Ensembl

ind <- match(rownames(gene.rank.order.5method.cosine.last.point.bic), rownames(annot.filter))
annot.filter.ordered <- annot.filter[ind, ]
gene.rank.order.5method.cosine.last.point.bic$ensembl <- annot.filter.ordered$Ensembl

# exprs matrices are the same except that fpkm.tumor.symbol.filter has seven duplicate ensembl ids
ensembl.ids <- intersect(annot.filter$Ensembl, rownames(brca))
annot.filter <- annot.filter[match(ensembl.ids, annot.filter$Ensembl),]
bic.trim.distribution.fit <- bic.trim.distribution.fit[match(ensembl.ids, names(bic.trim.distribution.fit))]

gene.rank.order.5method.cosine.last.point.bic <- gene.rank.order.5method.cosine.last.point.bic[match(ensembl.ids, gene.rank.order.5method.cosine.last.point.bic$ensembl),]
gene.zrange.fraction.cosine.last.point.bic <- gene.zrange.fraction.cosine.last.point.bic[match(ensembl.ids, gene.zrange.fraction.cosine.last.point.bic$ensembl),]

# table(outliers$distributions == bic.trim.distribution.fit)
# FALSE  TRUE
#   1    17695

