### bedr_search.R #############################################################
# Description
# overlap SV breakpoints with GeneHancer-v5.18 regions
# https://github.com/uclahs-cds/BoutrosLab/blob/master/Resources/code/R/omics.plotting.recCNA/trunk/omics.plotting.recCNA/R/build.prop.cna.segment.matrix.fast.R

### HISTORY ###################################################################
# Version	Date		Developer	Comments
# 0.01		2023-11-02	jlivingstone	initial code
# 0.02		2023-11-07	jlivingstone	add ICGC BRCA-EU segments

# bedtools has to be installed on node (conda activate r_env)
### PREAMBLE ##################################################################
library(bedr)
library(BoutrosLab.utilities)

# read in GeneHancer hg19 because OS arrays are hg19 based
gh <- read.delim(
	file = file.path('/hot/ref/database/GeneHancer-v5.18/original/hg19/', 'GeneHancer_hg19_v5.18.bed'),
	header = FALSE,
	as.is = TRUE
	)
colnames(gh) <- c('chr', 'start', 'end', 'ghid')

# but only look at elite positions first
elements <- read.delim(
	file = file.path('/hot/ref/database/GeneHancer-v5.18/original/GRCh38/', 'GeneHancer_AnnotSV_elements_v5.18.txt'),
	as.is = TRUE
	)
elite <- elements[which(elements$is_elite == TRUE), ]
gh.elite <- gh[na.omit((match(elite$GHid, gh$ghid))), ]
gh.ordered <- gh.elite[(order(sub('chr', '', gh.elite$chr), gh.elite$start, gh.elite$end)), ]

# remove elements that are too short < 20bp?
gh.filtered <- gh.ordered[(gh.ordered$end - gh.ordered$start > 20),]

if (dataset == 'cpcgene') {
	# read in sample specific CNAs
	segment <- read.delim(
		file = file.path('/hot/users/jlivingstone/', '2016-02-11_CPCG_segments_200pg.txt'),
		header = FALSE,
		as.is = TRUE
		)
	colnames(segment) <- c('sample', 'chr', 'start', 'end', 'logr')
	segment$status <- NA
	segment$status[which(segment$logr >= 0.20)] <- 1
	segment$status[which(segment$logr <= -0.20)] <- -1

	# remove segments that shoulnt be there - is stil NA
	cna <- segment[!is.na(segment$status),]
	}

if (dataset == 'icgc') {
	segment <- read.delim(
		file = file.path('/hot/users/jlivingstone/outlier/', '2023-11-06_ICGC_BRCA_EU_segments.txt'),
		as.is = TRUE
		)
	colnames(segment) <- c('sample', 'chr', 'start', 'end', 'status')
	# so format is the same as with cpcg
	cna <- data.frame(segment, logr = NA)
	}

samples <- unique(cna$sample)
cn.list <- list()
cols.to.remove <- c('sample', 'logr')
for (i in 1:length(samples)) {
	sample.ind <- which(cna$sample == samples[i])
	temp <- cna[sample.ind, -match(cols.to.remove, colnames(cna))]
	cn.list[[i]] <- temp[order(temp$chr, temp$start, temp$end),]
	}
names(cn.list) <- samples

gh.positions <- paste0(gh.filtered$chr, ':', gh.filtered$start, '-', gh.filtered$end)

overlap.matrix <- matrix(
	data = NA,
	nrow = length(gh.positions),
	ncol = length(samples)
	)

for (i in 1:length(samples)) {
	print(i)

	# only care about deletions ? because these bring enhancers closer to genes ?
	dels.ind <- which(cn.list[[i]]$status < 0)
	if (length(dels.ind) == 0) {
		next()
		}

	sample.dels <- paste0('chr', cn.list[[i]]$chr, ':', cn.list[[i]]$start, '-', cn.list[[i]]$end)
	dels.overlap <- in.region(gh.positions, sample.dels[dels.ind])
	overlap.matrix[,i] <- dels.overlap
	}
colnames(overlap.matrix) <- samples
rownames(overlap.matrix) <- gh.filtered$ghid

setwd('/hot/users/jlivingstone/outlier/')

if (dataset == 'cpcgene') {
	filename <- generate.filename('outlier', 'GeneHancer_CPCG_OS_overlap', 'txt')
	}
if (dataset == 'icgc') {
	filename <- generate.filename('outlier', 'GeneHancer_ICGC_segment_overlap', 'txt')
	}

write.table(
	x = overlap.matrix,
	file = filename,
	quote = FALSE,
	sep = '\t'
	)

# check if breakpoint is in enhancer region - included 573 promoters & 30026 promoter/enhancer
elements.filtered <- elements[match(rownames(overlap.matrix), elements$GHid),]

# each GHid has more than one score, since one enhancer can promoter more than one gene
# median = 5; mean = 6.1, min = 1, max = 119
scores <- read.delim(
	file = '/hot/ref/database/GeneHancer-v5.18/original/GRCh38/GeneHancer_AnnotSV_gene_association_scores_v5.18.txt',
	as.is = TRUE
	)
scores.expanded <- scores[scores$GHid %in% rownames(overlap.matrix),]

# n = 242,701
scores.elite <- scores.expanded[which(scores.expanded$is_elite == 1),]
scores.filtered <- scores.elite[order(scores.elite$GHid),]

# Let's look at TMPRSS2
tmprss2.ids <- scores.filtered$GHid[which(scores.filtered$symbol == 'TMPRSS2')]
tmprss2.matrix <- overlap.matrix[match(tmprss2.ids, rownames(overlap.matrix)), ]

clinical <- read.delim(
	file = '/hot/project/disease/ProstateTumor/PRAD-000050-666PG/data/clinical/2021-05-14_666pg_All_updated_AUS_clinical.txt',
	as.is = TRUE
	)
clinical.gh <- clinical[match(colnames(overlap.matrix), clinical$sample_id), ]
table(clinical.gh$ets_consensus == tmprss2.matrix[1,])
# 60% of samples agree, all disagreement is ets_consensus == TRUE but missed by gh overlap

# or close to enhancer region (+ 500kb)
