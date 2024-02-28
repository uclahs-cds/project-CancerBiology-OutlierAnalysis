## distance_to_element.R ##############################################################
# Description
# determine distance of breakpoint to GeneHancer element

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-02-28	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(BoutrosLab.utilities)
library(GenomicRanges)

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

#n = 125 samples with RNA & mut calls
sample.overlap <- intersect(outliers$sample.name, samples.with.muts)
outliers.parsed <- outliers[match(sample.overlap, outliers$sample.name), ]

results <- list()
for (i in 1:length(sample.overlap)) {
	print(i)
	muts <- rbind(
		inversion[inversion$submitted_sample_id %in% sample.overlap[i],],
		translocation[translocation$submitted_sample_id %in% sample.overlap[i],]
		)

	from.bkpt <- GRanges(
		seqnames = paste0('chr', muts$chr_from),
		strand = muts$chr_from_strand,
		ranges = IRanges(
			start = muts$chr_from_bkpt,
			end = muts$chr_from_bkpt + 1
			)
		)

        to.bkpt <- GRanges(
                seqnames = paste0('chr', muts$chr_to),
		strand = muts$chr_to_strand,
		ranges = IRanges(
			start = muts$chr_to_bkpt,
			end = muts$chr_to_bkpt + 1,
			)
		)

	to.closest.element <- data.frame(
		distanceToNearest(to.bkpt, gh.regions),
		stringsAsFactors = FALSE
		)

	from.closest.element <- data.frame(
		distanceToNearest(from.bkpt, gh.regions),
		stringsAsFactors = FALSE
		)

	# which distance in the smallest
	closest.element  <- from.closest.element
	ind.to.replace <- to.closest.element$distance < from.closest.element$distance
	closest.element[which(ind.to.replace),] <- to.closest.element[which(ind.to.replace),]
	
	# add gene information & co-ordinates
	closest.element$to.bkpt <- paste0(muts$chr_to, ':', muts$chr_to_bkpt)
	closest.element$from.bkpt <- paste0(muts$chr_from, ':', muts$chr_from_bkpt)
	closest.element$gh.gene <- gh$symbol[closest.element$subjectHits]
	closest.element$gh.coord <- paste0(gh$chr, ':', gh$element_start, '-', gh$element_end)[closest.element$subjectHits]
	closest.element$bkpt <- ifelse(ind.to.replace, 'to_bkpt', 'from_bkpt')

	results[[sample.overlap[i]]] <- closest.element
	}

save(
	results,
	file = generate.filename('outlier', 'distance_to_gh_element', 'rda')
	)

# only one PD18734a & LRRC69 match
#for (i in 1:length(results)) {
#	x <- match(outliers.parsed$outlier_genes[i], results[[i]]$gh.gene)
#	print(paste(i, ' ', x))
#	}
