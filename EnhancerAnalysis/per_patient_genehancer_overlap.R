### per_patient_genehancer_overlap.R ########################################################
# Description
# per patient that has an outlier gene, overlap SV segment with GH db

# Rscript per_patient_genehancer_overlap.R --mutation.type tandem
### HISTORY ##################################################################################
# Version	Date		Developer	Comments
# 0.01		2023-12-21	jlivingstone	initial code
# 0.02		2024-01-05	jlivingstone	add code to overlap with genehancer regions
# 0.03		2024-02-07	jlivingstone	update to use WGS deletion/duplication calls
						add conversion using liftover in rtracklayer
### PREAMBLE #################################################################################
library(bedr)
library(BoutrosLab.utilities)
library(getopt)
library(rtracklayer)

params <- matrix(
        data = c(
                'mutation.type', 'm', '0', 'character'
                ),
        ncol = 4,
        byrow = TRUE
        );

opt <- getopt(params);
mutation.type <- opt$mutation.type

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
gh <- read.delim(
	file = file.path('/hot/ref/database/GeneHancer-v5.18/processed/GRCh38', 'GeneHancer_AnnotSV_elements_v5.18_elite.txt'),
	as.is = TRUE
	)

path <- '/hot/ref/cohort/ICGC/BRCA/EU/processed/wgs/'

# interchromosome rearrangment - chr_from, chr_from_bkp & chr_to, chr_to_bkpt (each breakpoint is region ?)
# read in mutation of interest
if (mutation.type == 'gains') {
	muts <- read.delim(
		#file = file.path(path, 'gain_unique_somatic_mutation.BRCA-EU.tsv'),
		file = file.path(path, 'tandem_duplication.BRCA-EU.tsv'),
		as.is = TRUE
		)
	muts$chromosome <- muts$chr_to
	muts$chromosome_start <- muts$chr_from_bkpt
	muts$chromosome_end <- muts$chr_to_bkpt
	}

if (mutation.type == 'loss') {
	muts <- read.delim(
		#file = file.path(path, 'loss_unique_somatic_mutation.BRCA-EU.tsv'),
		file = file.path(path, 'deletion.BRCA-EU.tsv'),
		as.is = TRUE
		)
	muts$chromosome <- muts$chr_to
	muts$chromosome_start <- muts$chr_from_bkpt
	muts$chromosome_end <- muts$chr_to_bkpt
	}

# tandem duplication - chr_from, chr_from_bkp & chr_to, chr_to_bkpt (one region)
if (mutation.type == 'tandem') {
	muts <- read.delim(
		file = file.path(path, 'tandem_duplication.BRCA-EU.tsv'),
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

# ugh the WGS names aren't the same as RNA
outliers$sample.name <- sub('R', 'D', outliers$patient)

# n = 104 (57% of total RNA patients)
sample.overlap <- intersect(outliers$sample.name, muts$submitted_sample_id)

outliers.parsed <- outliers[match(sample.overlap, outliers$sample.name),]

elements.overlap <- list()
gh.overlap <- list()
engine <- 'granges' #'bedr'

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
		element.regions <- GRanges(
			seqnames = paste0('chr', elements$chr),
			ranges = IRanges(
				start = elements$element_start,
				end = elements$element_end
				)
			)

		if (engine == 'bedr') {
			element.regions <- paste0('chr', elements$chr, ':', elements$element_start, '-', elements$element_end)
			ind <- order(elements$chr, elements$element_start, elements$element_end)
			element.regions.ordered <- element.regions[ind]
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

		if (engine == 'bedr') {
			regions <- paste0('chr', sample.muts$chromosome, ':', sample.muts$chromosome_start, '-', sample.muts$chromosome_end)
			s.ind <- order(sample.muts$chromosome, sample.muts$chromosome_start, sample.muts$chromosome_end)
			regions.ordered <- regions[s.ind]

			}

		# check if mutation regions overlap with genehancer elements (object a are in object b)
		regions.overlap <- as.data.frame(mergeByOverlaps(element.regions, regions))
		regions.overlap.index <- findOverlaps(element.regions, regions)

		if (length(regions.overlap.index) > 0) {
			elements.overlap[[outliers.parsed$sample.name[i]]] <- elements[unique(queryHits(regions.overlap.index)),]
			gh.overlap[[outliers.parsed$sample.name[i]]] <- sample.muts[unique(subjectHits(regions.overlap.index)), ]
			}

		if (engine == 'bedr') {
			overlap <- element.regions.ordered[in.region(element.regions.ordered, regions.ordered)]
			elements.overlap <- elements[match(overlap, element.regions), ]
			gh.overlap <- sample.muts[in.region(regions.ordered, element.regions.ordered), ]
			}
		}
	}
save(
	list = c('elements.overlap', 'gh.overlap'),
	file = generate.filename(paste0('genehancer_overlap_regions_', mutation.type), 'BRCA_EU', 'rda')
	)
