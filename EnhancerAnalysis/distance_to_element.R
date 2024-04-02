## distance_to_element.R ##############################################################
# Description
# determine distance of breakpoint to GeneHancer element

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-02-28	jlivingstone	initial code
# 0.02		2024-03-12	jlivingstone	add chromosome max value; heatmap
# 0.03		2024-03-28	jlivingstone	calculate distance if gh element
#						is within breakpoints

### PREAMBLE ###########################################################################
library(BoutrosLab.utilities)
library(GenomicRanges)

source('/hot/code/jlivingstone/GitHub/uclahs-cds/project-ProstateCancer-m6A/plot.general.contingency.table.R');

setwd('/hot/user/jlivingstone/outlier/enhancer_analysis/distance_analysis')

gh <- read.delim(
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19',
		'GeneHancer_AnnotSV_hg19_elements_v5.18_elite.txt'
		),
	as.is = TRUE
	)

gh.regions <- GRanges(
	seqnames = paste0('chr', gh$chr),
	ranges = IRanges(
		start = gh$element_start,
		end = gh$element_end
		)
	)

inversion <- read.delim(
	file = file.path(
		'/hot/ref/cohort/ICGC/BRCA/EU/processed/wgs/',
		'inversion.BRCA-EU.tsv'
		),
	as.is = TRUE
	)
# to account for trailing numbers
inversion$submitted_sample_id <- sub('a2', 'a', inversion$submitted_sample_id)
inversion$submitted_sample_id <- sub('a3', 'a', inversion$submitted_sample_id)

translocation <- read.delim(
	file = file.path(
		'/hot/ref/cohort/ICGC/BRCA/EU/processed/wgs/',
		'interchromosomal_rearrangement.BRCA-EU.tsv'
		),
	as.is = TRUE
	)
# to account for trailing numbers
translocation$submitted_sample_id <- sub('a2', 'a', translocation$submitted_sample_id)
translocation$submitted_sample_id <- sub('a3', 'a', translocation$submitted_sample_id)

samples.with.muts <- union(unique(translocation$submitted_sample_id), unique(inversion$submitted_sample_id))

# read in outlier samples and genes
outliers <- read.delim(
	file = file.path(
		'/hot/user/jlivingstone/outlier/run_method/',
		'2023-12-21_Outlier_patients_with_genes_BRCA_EU_cutoff_0.01.txt'
		),
	as.is = TRUE
	)

outliers$sample.name <- sub('R', 'D', outliers$patient)
outliers$sample.name <- sub('.RNA', '', outliers$sample.name)
outliers$sample.name[grep('PD6418a.2', outliers$sample.name)] <- 'PD6418a'

index <- read.delim(
	file = '/hot/ref/reference/GRCh37-EBI-GENCODE24/GRCh37.primary_assembly.genome.fa.fai',
	as.is = TRUE,
	header = FALSE
	)
colnames(index) <- c('chr', 'length', 'offset', 'nchar', 'nbytes')

exprs <- read.delim(
        file = '/hot/user/jlivingstone/outlier/NikZainal_2016/original/SupplementaryTable7Transcriptomic342.txt',
        as.is = TRUE
        )
colnames(exprs) <- sub('R', 'D', colnames(exprs))
colnames(exprs) <- sub('.RNA', '', colnames(exprs))
colnames(exprs)[grep('PD6418a.2', colnames(exprs))] <- 'PD6418a'

sample.overlap <- intersect(colnames(exprs), samples.with.muts)
samples.to.include <- intersect(sample.overlap, outliers$sample.name)
outliers.parsed <- outliers[match(samples.to.include, outliers$sample.name), ]

# analysis on outlier level
# distance between mutation and gh element will be smaller in outlier patients

# gh elements for outlier gene
outlier.genes <- unlist(strsplit(outliers.parsed$outlier_genes, ';'))
overlap.genes <- intersect(outlier.genes, unique(gh$symbol))

get.distance <- function(gh.subset, bkpts, muts) {
	closest.element <- data.frame(
		distanceToNearest(gh.subset, bkpts),
		stringsAsFactors = FALSE
		)

	if (nrow(closest.element > 1)) {
		# add gene information & co-ordinates
		closest.element$chr.from.bkpt <- muts$chr_from[closest.element$subjectHits]
		closest.element$from.bkpt <- muts$chr_from_bkpt[closest.element$subjectHits]
		closest.element$chr.to.bkpt <- muts$chr_to[closest.element$subjectHits]
		closest.element$to.bkpt <- muts$chr_to_bkpt[closest.element$subjectHits]
		closest.element$gh.gene <- gh.gene$symbol[closest.element$queryHits]
		closest.element$gh.coord.chr <- gh.gene$chr[closest.element$queryHits]
		closest.element$gh.coord.start <- gh.gene$element_start[closest.element$queryHits]
		closest.element$gh.coord.end <- gh.gene$element_end[closest.element$queryHits]

		gh.within <- closest.element$from.bkpt < closest.element$gh.coord.start & closest.element$to.bkpt > closest.element$gh.coord.end
		closest.element <- closest.element[which(gh.within == TRUE), ]
		}

	# after filtering out distances that exist but the gh element isn't within the mutation breakpoints
	if (nrow(closest.element) == 0) {
		# if no mutations are close to element, GRanges will return an empty vector so need to account for this
		closest.element <- data.frame(	
			queryHits = 0,
			subjectHits = 0,
			distance = index$length[match(paste0('chr', unique(gh.gene$chr)), index$chr)],
			chr.from.bkpt = unique(gh.gene$chr),
			from.bkpt = NA,
			chr.to.bkpt = unique(gh.gene$chr),
			to.bkpt = NA,
			gh.gene = gh.gene$symbol[which.max(gh.gene$combined_score)],
			gh.coord.chr = gh.gene$chr[which.max(gh.gene$combined_score)],
			gh.coord.start = gh.gene$element_start[which.max(gh.gene$combined_score)],
			gh.coord.end = gh.gene$element_end[which.max(gh.gene$combined_score)],
			stringsAsFactors = FALSE
			)
		}

	outcome <- closest.element[which.min(closest.element$distance),]
	return(outcome)
	}

distance.results <- list()
for (j in 1:length(overlap.genes)) {
	print(overlap.genes[j])

	gh.gene <- gh[which(gh$symbol == overlap.genes[j]),]

	gh.subset <- GRanges(
		seqnames = paste0('chr', gh.gene$chr),
		ranges = IRanges(
			start = gh.gene$element_start,
			end = gh.gene$element_end
			)
		)

	# analysis on patient level - distance between mutation and gh elements
	# just keep the closest mutation / gh element
	results <- data.frame()
	for (i in 1:length(samples.to.include)) {
		muts <- rbind(
			inversion[inversion$submitted_sample_id %in% samples.to.include[i],],
			translocation[translocation$submitted_sample_id %in% samples.to.include[i],]
			)

		#start
		from.bkpt <- GRanges(
			seqnames = paste0('chr', muts$chr_from),
			strand = muts$chr_from_strand,
			ranges = IRanges(
				start = muts$chr_from_bkpt,
				end = muts$chr_from_bkpt + 1
				)
			)

		from.results <- get.distance(
			gh.subset = gh.subset,
			bkpts = from.bkpt,
			muts = muts
			)

		#end
		to.bkpt <- GRanges(
			seqnames = paste0('chr', muts$chr_to),
			strand = muts$chr_to_strand,
			ranges = IRanges(
				start = muts$chr_to_bkpt,
				end = muts$chr_to_bkpt + 1
				)
			)

		to.results <- get.distance(
			gh.subset = gh.subset,
			bkpts = to.bkpt,
			muts = muts
			)

		# compare start.results and end.results to get shortest distance
		outcome <- rbind(to.results, from.results)
		results <- rbind(results, outcome[which.min(outcome$distance), ])
		}
	rownames(results) <- samples.to.include
	distance.results[[overlap.genes[j]]] <- results
	}

save(
	distance.results,
	file = generate.filename('outlier', 'distance_of_mutation_to_outlier_gene_gh_element', 'rda')
	)

# different scenarios - no mutation on same chromosome so no distance returned)
# distance = -1, so remove ?

true.outliers <- data.frame()
for (i in 1:nrow(outliers.parsed)) {
        temp <- data.frame(
                sample.name = outliers.parsed$sample.name[i],
                outlier_genes = unlist(strsplit(outliers.parsed$outlier_genes[i], ';')),
                stringsAsFactors = FALSE
                )
        true.outliers <- rbind(true.outliers, temp)
        }
# remove genes that aren't in GH
true.outliers <- true.outliers[match(overlap.genes,true.outliers$outlier_genes),]

# compare outlier gene distance across all outliers (n = 179) ? against non outliers (22,375)
outlier.distance <- vector()
nonoutlier.distance <- vector()
outlier.chr.max.flag <- vector()
nonoutlier.chr.max.flag <- vector()
for (i in 1:length(distance.results)) {
	print(i)
	# get samples with outlier
	outlier.sample <- true.outliers$sample.name[match(names(distance.results)[i], true.outliers$outlier_genes)]
	outlier.distance <- append(
		x = outlier.distance,
		values = distance.results[[i]]$distance[match(outlier.sample, rownames(distance.results[[i]]))]
		)
	outlier.chr.max.flag <- append(
		x = outlier.chr.max.flag,
		values = distance.results[[i]]$queryHits[match(outlier.sample, rownames(distance.results[[i]]))]
		)
	nonoutlier.distance <- append(
		x = nonoutlier.distance,
		values = distance.results[[i]]$distance[-match(outlier.sample, rownames(distance.results[[i]]))]
		)
	nonoutlier.chr.max.flag <- append(
		x = nonoutlier.chr.max.flag,
		values = distance.results[[i]]$queryHits[-match(outlier.sample, rownames(distance.results[[i]]))]
		)
	}

n.nonoutlier.nomatch <- length(which(nonoutlier.chr.max.flag == 0))
n.nonoutlier.match <- length(which(nonoutlier.chr.max.flag != 0))

n.outlier.nomatch <- length(which(outlier.chr.max.flag == 0))
n.outlier.match <- length(which(outlier.chr.max.flag != 0))

x <- matrix(
	data = c(n.nonoutlier.nomatch, n.nonoutlier.match, n.outlier.nomatch, n.outlier.match),
	ncol = 2,
	nrow = 2,
	byrow = FALSE
	)
colnames(x) <- c('Non-outlier', 'Outlier')
rownames(x) <- c('No match', 'Match')
class(x) <- 'table'
chisq.test(x)

wilcox.test(
	x = outlier.distance,
	y = nonoutlier.distance
	)

plot.general.contingency.table(
	plot.data = x,
	filename = generate.filename('outlier', 'outlier_vs_nonoutlier_distance_to_gh_element', 'png'),
	ylabel = 'On same chromosome',
	xlabel = 'Gene Status',
	scheme = 'count',
	text.cex = 0.75
	);

p.value <- rep(NA, nrow(true.outliers))
for (i in 1:nrow(true.outliers)) {
	temp <- distance.results[[i]]
	outlier.value <- temp$distance[rownames(temp) == true.outliers$sample.name[i]]
	nonoutlier.values <- temp$distance[rownames(temp) != true.outliers$sample.name[i]]
	p.value[i] <- length(which(nonoutlier.values < outlier.value)) / nrow(temp)
	}
q.value <- p.adjust(p.value, method = 'fdr')

write.table(
	x = data.frame(p.value, q.value),
	file = generate.filename('outlier', 'distance_to_gh_element_statistics', 'tsv'),
	quote = FALSE,
	sep = '\t'
	)

#distance.results[[which(names(distance.results) == 'CLCN4')]][true.outliers$sample.name[which(names(distance.results) == 'CLCN4')],]
