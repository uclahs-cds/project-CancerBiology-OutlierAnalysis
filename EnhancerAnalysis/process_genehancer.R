### explore_genehancer.R ###########################################################
# Description
# Create data files under /hot/ref/database/GeneHancer-v5.18/processed/GRCh38

### HISTORY ########################################################################
# Version	Date		Developer	Comments
# 0.01		2024-02-08	jlivingstone	initial code

### PREAMBLE ####################################################################

setwd('/hot/ref/database/GeneHancer-v5.18/original/GRCh38/')

# contains the same information as in GeneHancer_v5.18.bed + element info (n = 409,271)
elements <- read.delim(
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/original/GRCh38/',
		'GeneHancer_AnnotSV_elements_v5.18.txt'
		),
	as.is = TRUE
	)

# contains scores for each GHid (n = 2,498,196)
scores <- read.delim(
        file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/original/GRCh38/',
		'GeneHancer_AnnotSV_gene_association_scores_v5.18.txt'
		),
        as.is = TRUE
        )

# need to combine scores with elements to get gene symbol and genomic co-ordinates
expanded <- list()
for (i in 1:nrow(elements)) {
	if (i %% 1000 == 0) {
		print(i)
		}
	temp <- cbind(elements[i,c('chr', 'element_start', 'element_end', 'regulatory_element_type')], scores[which(scores$GHid %in% elements$GHid[i]),])
	expanded <- rbind(expanded, temp)
	}

# n = 2,498,196
write.table(
	x = expanded,
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19', 
		'GeneHancer_AnnotSV_elements_v5_18.txt'
		),
	quote = FALSE,
	sep = '\t',
	row.names = FALSE
	)

# elite only
elite <- expanded[which(expanded$is_elite == 1), ]

write.table(
	x = elite,
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19',
		'GeneHancer_AnnotSV_elements_v5.18_elite.txt'
		),
	quote = FALSE,
	sep = '\t',
	row.names = FALSE
	)

# use GHids to create same files in hg19 instead of using liftover
hg19 <- read.delim(
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19',
		'GeneHancer_AnnotSV_hg19_v5.18.txt'
		),
	as.is = TRUE,
	)
colnames(hg19) <- c('chr', 'element_start', 'element_end', 'GHid')
# length(unique(hg19$GHid)) - n = 408,144

grch <- read.delim(
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/GRCh38',
		'GeneHancer_AnnotSV_elements_v5.18.txt'
		),
	as.is = TRUE
	)
#length(unique(grch$GHid)) - n = 409,271

hg19.expanded <- hg19[match(grch$GHid, hg19$GHid),]
liftover <- grch
liftover$element_start <- hg19.expanded$element_start
liftover$element_end <- hg19.expanded$element_end

toprint <- liftover[-which(is.na(hg19.expanded$GHid)), ]
write.table(
	x = toprint,
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19',
		'GeneHancer_AnnotSV_hg19_elements_v5.18.txt'
		),
	quote = FALSE,
	sep = '\t',
	row.names = FALSE
	)

elite <- toprint[which(toprint$is_elite == 1),]
write.table(
	x = elite,
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19',
		'GeneHancer_AnnotSV_hg19_elements_v5.18_elite.txt'
		),
	quote = FALSE,
	sep = '\t',
	row.names = FALSE
	)
