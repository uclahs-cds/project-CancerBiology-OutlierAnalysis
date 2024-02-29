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

# analysis on outlier level
# distance between mutation and gh element will be smaller in outlier patients

# gh elements for outlier gene
outlier.genes <- unlist(strsplit(outliers.parsed$outlier_genes, ';'))
overlap.genes <- intersect(outlier.genes, unique(gh$symbol))

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

		bkpts <- c(to.bkpt, from.bkpt)
		muts.combined <- rbind(muts, muts)
		closest.element <- data.frame(
			distanceToNearest(gh.subset, bkpts),
			stringsAsFactors = FALSE
			)

		if (nrow(closest.element > 1)) {
			# add gene information & co-ordinates
			closest.element$to.bkpt <- paste0(muts.combined$chr_to[closest.element$subjectHits], ':', muts.combined$chr_to_bkpt[closest.element$subjectHits])
			closest.element$from.bkpt <- paste0(muts.combined$chr_from[closest.element$subjectHits], ':', muts.combined$chr_from_bkpt[closest.element$subjectHits])
			closest.element$gh.gene <- gh.gene$symbol[closest.element$queryHits]
			closest.element$gh.coord <- paste0(gh.gene$chr, ':', gh.gene$element_start, '-', gh.gene$element_end)[closest.element$queryHits]
		} else {
			# if no mutations are close to element, GRanges will return an empty vector
			closest.element <- results[1,]
			closest.element$queryHits <- 0
			closest.element$subjectHits <- 0
			closest.element$distance <- -1
			closest.element$to.bkpt <- NA
			closest.element$from.bkpt <- NA
			closest.element$gh.gene <- gh.gene$symbol[1]
			closest.element$gh.coord <- paste0(gh.gene$chr, ':', gh.gene$element_start, '-', gh.gene$element_end)[1]
			}
		results <- rbind(results, closest.element[which.min(closest.element$distance),])
		}
	distance.results[[overlap.genes[j]]] <- results
	}

save(
	distance.results,
	file = generate.filename('outlier', 'distance_of_mutation_to_outlier_gene_gh_element', 'rda')
	)
