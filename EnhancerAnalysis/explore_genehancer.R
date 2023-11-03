### explore_genehancer.R ###########################################################
# Description
# get some stats on genehancer db, data is distributed across four files
# genehancer, elements, scores, tissues

### HISTORY ########################################################################
# Version	Date		Developer	Comments
# 0.01		2023-01-02	jlivingstone	initial code

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
# 371839      7370         30062 

# collect only the elite GHids to begin with, then merge (n = 173,025)
elite <- elements[which(elements$is_elite == 1), ]

scores <- read.delim(
	file = '/hot/ref/database/GeneHancer-v5.18/original/GRCh38/GeneHancer_AnnotSV_gene_association_scores_v5.18.txt',
	as.is = TRUE
	)
# number of genes per GHid median = 3, mean = 9.5, min = 1, max = 413

scores.elite <- scores[which(scores$is_elite == 1),]

tissue <- read.delim(
	file = 'GeneHancer_AnnotSV_tissues_v5.18.txt',
	as.is = TRUE,
	sep = '\t'
	)
# dbSUPER  ENCODE Ensembl FANTOM5   VISTA 
# 365378 1292029 3844417   55407    1903

# number of tissues per GHid median = 9, mean = 15, min = 1, max = 311

#plan of attack, overlap SV breakpoints with GH regions, then see if scores contains outlier gene

