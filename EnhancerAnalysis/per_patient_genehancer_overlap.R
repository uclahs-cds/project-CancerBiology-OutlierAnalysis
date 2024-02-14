### per_patient_genehancer_overlap.R ########################################################
# Description
# per patient that has an outlier gene, overlap SV segment with GH db

# Rscript per_patient_genehancer_overlap.R --mutation.type tandem_duplication --elite.only 0 --algorithm bedr
### HISTORY ##################################################################################
# Version	Date		Developer	Comments
# 0.01		2023-12-21	jlivingstone	initial code
# 0.02		2024-01-05	jlivingstone	add code to overlap with genehancer regions
# 0.03		2024-02-07	jlivingstone	update to use WGS deletion/duplication calls
#						add conversion using liftover in rtracklayer
### PREAMBLE #################################################################################
library(bedr)
library(BoutrosLab.utilities)
library(getopt)
library(rtracklayer)

params <- matrix(
	data = c(
		'mutation.type', 'm', '0', 'character',
		'elite.only', 'e', '0', 'numeric',
		'algorithm', 'a', '0', 'character'
		),
	ncol = 4,
	byrow = TRUE
	);

opt <- getopt(params);
mutation.type <- opt$mutation.type
elite.only <- opt$elite.only
engine <- opt$algorithm

setwd('/hot/user/jlivingstone/outlier/enhancer_analysis')

convert.my.regions <- function(regions, chain) {
	results <- as.data.frame(
		liftOver(
			x = regions,
			chain = chain
			),
		stringsAsFactors = FALSE
		)
	segment.length <- regions@ranges@width

	need.to.pick <- as.numeric(names(which(table(results$group) > 1)))
	need.segment.length <- segment.length[need.to.pick]

	to.keep <- results[match(as.numeric(names(which(table(results$group) == 1))), results$group),]

	# keep region with the shortest difference in length
	for (i in 1:length(need.to.pick)) {
		temp <- results[which(results$group == need.to.pick[i]),]
		to.keep <- rbind(to.keep, temp[which.min(abs(temp$width - need.segment.length[i])),])	
		}
	to.keep.ordered <- to.keep[order(to.keep$group), ]

	regions <- GRanges(
		seqnames = to.keep.ordered$seqnames,
		ranges = IRanges(
			start = to.keep.ordered$start,
			end = to.keep.ordered$end
			)
		)
	return(regions)
	}

# read in outlier genes per patient
outliers <- read.delim(
	file = file.path('/hot/user/jlivingstone/outlier/run_method/', '2023-12-21_Outlier_patients_with_genes_BRCA_EU_cutoff_0.01.txt'),
	as.is = TRUE
	)

# read in GeneHancer elements
if (elite.only  == 1) {
	gh <- read.delim(
		file = file.path('/hot/ref/database/GeneHancer-v5.18/processed/GRCh38', 'GeneHancer_AnnotSV_elements_v5.18_elite.txt'),
		as.is = TRUE
		)
} else {
	gh <- read.delim(
		file = file.path('/hot/ref/database/GeneHancer-v5.18/processed/GRCh38', 'GeneHancer_AnnotSV_elements_v5.18.txt'),
		as.is = TRUE
		)
	}

path <- '/hot/ref/cohort/ICGC/BRCA/EU/processed/wgs/'

# interchromosome rearrangment - chr_from, chr_from_bkp & chr_to, chr_to_bkpt (each breakpoint is region ?)
# read in mutation of interest
# These mutation calls are from SNP arrays
if (mutation.type == 'gains') {
	muts <- read.delim(
		file = file.path(path, 'gain_unique_somatic_mutation.BRCA-EU.tsv'),
		as.is = TRUE
		)
	}

# These mutation calls are from SNP arrays
if (mutation.type == 'loss') {
	muts <- read.delim(
		file = file.path(path, 'loss_unique_somatic_mutation.BRCA-EU.tsv'),
		as.is = TRUE
		)
	}

if (mutation.type == 'tandem_duplication') {
	muts <- read.delim(
		file = file.path(path, 'tandem_duplication.BRCA-EU.tsv'),
		as.is = TRUE
		)
	muts$chromosome <- muts$chr_to
	muts$chromosome_start <- muts$chr_from_bkpt
	muts$chromosome_end <- muts$chr_to_bkpt
	}

if (mutation.type == 'deletion') {
	muts <- read.delim(
		file = file.path(path, 'deletion.BRCA-EU.tsv'),
		as.is = TRUE
		)
	muts$chromosome <- muts$chr_to
	muts$chromosome_start <- muts$chr_from_bkpt
	muts$chromosome_end <- muts$chr_to_bkpt
	}

# inversion - chr_from, chr_from_bkp & chr_to, chr_to_bkpt (one region, reverse co-ordinates)
if (mutation.type == 'inversion') {
	muts <- read.delim(
		file = file.path(path, 'inversion.BRCA-EU.tsv'),
		as.is = TRUE
		)
	muts$chromosome <- muts$chr_from
	muts$chromosome_start <- muts$chr_from_bkpt
	muts$chromosome_end <- muts$chr_to_bkpt
	}

# interchromosomal_rearrangement - need to figure out how to deal with breakpoints on different chromosomes
if (mutation.type == 'interchromosomal_rearrangement') {
	muts <- read.delim(
		file = file.path(path, 'interchromosomal_rearrangement.BRCA-EU.tsv'),
		as.is = TRUE
		)
	# add code here
	}

# ugh the WGS names aren't the same as RNA
outliers$sample.name <- sub('R', 'D', outliers$patient)

# n = 104 (57% of total RNA patients)
sample.overlap <- intersect(outliers$sample.name, muts$submitted_sample_id)

outliers.parsed <- outliers[match(sample.overlap, outliers$sample.name), ]

elements.overlap <- list()
gh.overlap <- list()

chain <- import.chain(con = '/hot/user/jlivingstone/outlier/enhancer_analysis/hg19ToHg38.over.chain')

# for each patient, overlap the enhancer regions with sv breakpoints
for (i in 1:nrow(outliers.parsed)) {
	print(i)
	genes <- unlist(
		strsplit(
			x = outliers.parsed$outlier_genes[i],
			split = ';'
			)
		)
	elements <- gh[gh$symbol %in% genes, ]

	if (nrow(elements) > 0) {
		# GH regions
		if (engine == 'bedr') {
			element.regions <- paste0('chr', elements$chr, ':', elements$element_start, '-', elements$element_end)
			ind <- order(elements$chr, elements$element_start, elements$element_end)
			element.regions.ordered <- element.regions[ind]
		} else {
			element.regions <- GRanges(
				seqnames = paste0('chr', elements$chr),
				ranges = IRanges(
					start = elements$element_start,
					end = elements$element_end
					)
				)
			}

		sample.muts <- muts[muts$submitted_sample_id %in% outliers.parsed$sample.name[i], ]

		regions.grch38 <- GRanges(
			seqnames = paste0('chr', sample.muts$chromosome),
			ranges = IRanges(
				start = sample.muts$chromosome_start,
				end = sample.muts$chromosome_end
				)
			)

		# liftover regions
		regions = convert.my.regions(
			regions = regions.grch38,
			chain = chain
			)

		# check if mutation regions overlap with genehancer elements (object a are in object b)

		if (engine == 'bedr') {
			sample.muts.converted <- as.data.frame(regions, stringsAsFactors = FALSE)
			regions <- paste0(sample.muts.converted$seqnames, ':', sample.muts.converted$start, '-', sample.muts.converted$end)
			s.ind <- order(sub('chr', '', sample.muts.converted$seqnames), sample.muts.converted$start, sample.muts.converted$end)
			regions.ordered <- regions[s.ind]

			overlap <- element.regions.ordered[in.region(element.regions.ordered, regions.ordered)]
			elements.overlap[[outliers.parsed$sample.name[i]]] <- elements[match(overlap, element.regions), ]
			gh.overlap[[outliers.parsed$sample.name[i]]] <- sample.muts[in.region(regions.ordered, element.regions.ordered), ]
		} else {
			regions.overlap <- as.data.frame(mergeByOverlaps(element.regions, regions))
			regions.overlap.index <- findOverlaps(element.regions, regions)
			elements.overlap[[outliers.parsed$sample.name[i]]] <- elements[unique(queryHits(regions.overlap.index)), ]
			gh.overlap[[outliers.parsed$sample.name[i]]] <- sample.muts[unique(subjectHits(regions.overlap.index)), ]
			}
	} else {
		elements.overlap[[outliers.parsed$sample.name[i]]] <- 'No overlapping genes in GH'
		gh.overlap[[outliers.parsed$sample.name[i]]] <- 'No overlapping genes in GH'
		}
	}
save(
	list = c('elements.overlap', 'gh.overlap'),
	file = generate.filename(paste0('genehancer_overlap_regions_', mutation.type, '_', engine), 'BRCA_EU', 'rda')
	)
