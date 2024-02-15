### per_patient_genehancer_overlap.R ########################################################
# Description
# per patient that has an outlier gene, overlap SV segment with GH db

# Rscript per_patient_genehancer_overlap.R --mutation.type tandem_duplication --elite.only 0 --algorithm bedr --input.file --flag ''
# '/hot/user/jlivingstone/outlier/run_method/2023-12-21_Outlier_patients_with_genes_BRCA_EU_cutoff_0.01.txt'

### HISTORY ##################################################################################
# Version	Date		Developer	Comments
# 0.01		2023-12-21	jlivingstone	initial code
# 0.02		2024-01-05	jlivingstone	add code to overlap with genehancer regions
# 0.03		2024-02-07	jlivingstone	update to use WGS deletion/duplication calls
#						add conversion using liftover in rtracklayer
# 0.04		2024-02-15	jlivingstone	add permuation datasets
### PREAMBLE #################################################################################
library(bedr)
library(BoutrosLab.utilities)
library(getopt)
library(rtracklayer)

params <- matrix(
	data = c(
		'mutation.type', 'm', '1', 'character',
		'elite.only', 'e', '1', 'numeric',
		'algorithm', 'a', '1', 'character',
		'input.file', 'i', '1', 'character',
		'flag', 'f', '1', 'character'
		),
	ncol = 4,
	byrow = TRUE
	);

opt <- getopt(params);
mutation.type <- opt$mutation.type
elite.only <- opt$elite.only
engine <- opt$algorithm
input.file <- opt$input.file
# set when running permutation test ie the seed used
flag <- opt$flag

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
	file = input.file,
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
# add pading around breakpoint?
if (mutation.type == 'inversion') {
	muts <- read.delim(
		file = file.path(path, 'inversion.BRCA-EU.tsv'),
		as.is = TRUE
		)
	muts$chromosome <- muts$chr_from
	muts$chromosome_start <- muts$chr_from_bkpt
	muts$chromosome_end <- muts$chr_to_bkpt
	}

# interchromosomal_rearrangement - duplicate and add 500bp to each breakpoint
if (mutation.type == 'interchromosomal_rearrangement') {
	temp <- read.delim(
		file = file.path(path, 'interchromosomal_rearrangement.BRCA-EU.tsv'),
		as.is = TRUE
		)
	# lets add 500kb before and after breakpoint
	break_one <- break_two <- temp
	break_one$chromosome <- temp$chr_from
	break_one$chromosome_start <- temp$chr_from_bkpt - 1000
	break_one$chromosome_end <- temp$chr_from_bkpt + 1000

	break_two$chromosome <- temp$chr_to
	break_two$chromosome_start <- temp$chr_to_bkpt - 1000
	break_two$chromosome_end <- temp$chr_to_bkpt + 1000

	muts <- rbind(break_one, break_two)
	}

# ugh the WGS names aren't the same as RNA
outliers$sample.name <- sub('R', 'D', outliers$patient)
outliers$sample.name <- sub('.RNA', '', outliers$sample.name)
outliers$sample.name[grep('PD6418a.2', outliers$sample.name)] <- 'PD6418a'

# to account for trailing numbers
muts$submitted_sample_id <- sub('a2', 'a', muts$submitted_sample_id)
muts$submitted_sample_id <- sub('a3', 'a', muts$submitted_sample_id)

# n = 120 / 182 (66% of total RNA patients)
sample.overlap <- intersect(outliers$sample.name, muts$submitted_sample_id)
outliers.parsed <- outliers[match(sample.overlap, outliers$sample.name), ]

chain <- import.chain(con = '/hot/user/jlivingstone/outlier/enhancer_analysis/hg19ToHg38.over.chain')

elements.overlap <- list()
gh.overlap <- list()
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
	file = generate.filename(paste0('genehancer_overlap_regions_', mutation.type, '_', engine, flag), 'BRCA_EU', 'rda')
	)
