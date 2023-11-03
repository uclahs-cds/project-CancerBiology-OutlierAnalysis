### bedr_search.R #############################################################
# Description
# overlap SV breakpoints with GeneHancer-v5.18 regions
# https://github.com/uclahs-cds/BoutrosLab/blob/master/Resources/code/R/omics.plotting.recCNA/trunk/omics.plotting.recCNA/R/build.prop.cna.segment.matrix.fast.R

### HISTORY ###################################################################
# Version	Date		Developer	Comments
# 0.01		2023-01-02	jlivingstone	initial code

### PREAMBLE ##################################################################
library(bedr)
library(BoutrosLab.utilities)

# read in GeneHancer hg19 because OS arrays are hg19 based
gh <- read.delim(
	file = '/hot/ref/database/GeneHancer-v5.18/original/hg19/GeneHancer_hg19_v5.18.bed',
	header = FALSE,
	as.is = TRUE
	)
colnames(gh) <- c('chr', 'start', 'end', 'ghid')

# but only look at elite positions first
elements <- read.delim(
	file = '/hot/ref/database/GeneHancer-v5.18/original/GRCh38/GeneHancer_AnnotSV_elements_v5.18.txt',
	as.is = TRUE
	)
elite <- elements[which(elements$is_elite == TRUE), ]
gh.elite <- gh[na.omit((match(elite$GHid, gh$ghid))), ]
gh.ordered <- gh.elite[(order(sub('chr', '', gh.elite$chr), gh.elite$start, gh.elite$end)), ]

# remove elements that are too short < 20bp?
gh.filtered <- gh.ordered[(gh.ordered$end - gh.ordered$start > 20),]

# read in sample specific CNAs
segment <- read.delim(
	file = '/hot/users/jlivingstone/2016-02-11_CPCG_segments_200pg.txt',
	header = FALSE,
	as.is = TRUE
	)
colnames(segment) <- c('sample', 'chr', 'start', 'end', 'logr')
segment$status <- NA
segment$status[which(segment$logr >= 0.20)] <- 1
segment$status[which(segment$logr <= -0.20)] <- -1

# remove segments that shoulnt be there - is stil NA
cna <- segment[!is.na(segment$status),]

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
write.table(
	x = overlap.matrix,
	file = generate.filename('outlier', 'GeneHancer_CPCG_OS_overlap', 'txt'),
	quote = FALSE,
	sep = '\t'
	)

# check if breakpoint is in enhancer region
# or close to enhancer region (+ 500bp?)

# read up on enhancer hijacking methods
