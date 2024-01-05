### per_patient_genehancer_overlap.R ##########################################
# Description
# per patient that has an outlier gene, overlap SV segment with GH db

### HISTORY ###################################################################
# Version	Date		Developer	Comments
# 0.01		2023-12-21	jlivingstone	initial code

### PREAMBLE ##################################################################
library(bedr)
library(BoutrosLab.utilities)
library(rtracklayer)

setwd('/hot/user/jlivingstone/outlier/enhancer_analysis')

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

# start with amplifications - gain in enhancer or promoter = higher abundance
gains <- read.delim(
	file = file.path('/hot/ref/cohort/ICGC/BRCA/EU/processed/', 'gain_unique_somatic_mutation.BRCA-EU.tsv'),
	as.is = TRUE
	)

# ugh the WGS names aren't the same as RNA
outliers$sample.name <- sub('R', 'D', outliers$patient)

# n = 104 (57% of total RNA patients)
sample.overlap <- intersect(outliers$sample.name, gains$submitted_sample_id)

outliers.parsed <- outliers[match(sample.overlap, outliers$sample.name),]

elements.overlap <- list()
gh.overlap <- list()
engine <- 'granges' #'bedr'
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

		sample.gains <- gains[gains$submitted_sample_id %in% outliers.parsed$sample.name[i], ]

		regions <- GRanges(
			seqnames = paste0('chr', sample.gains$chromosome),
			ranges = IRanges(
				start = sample.gains$chromosome_start,
				end = sample.gains$chromosome_end
				)
			)

		if (engine == 'bedr') {
			regions <- paste0('chr', sample.gains$chromosome, ':', sample.gains$chromosome_start, '-', sample.gains$chromosome_end)
			s.ind <- order(sample.gains$chromosome, sample.gains$chromosome_start, sample.gains$chromosome_end)
			regions.ordered <- regions[s.ind]

			}

		# check if gain regions overlap with genehancer elements (object a are in object b)
		regions.overlap <- as.data.frame(mergeByOverlaps(element.regions, regions))
		regions.overlap.index <- findOverlaps(element.regions, regions)

		if (length(regions.overlap.index) > 0) {
			elements.overlap[[outliers.parsed$sample.name[i]]] <- elements[queryHits(regions.overlap.index),]
			gh.overlap[[outliers.parsed$sample.name[i]]] <- sample.gains[unique(subjectHits(regions.overlap.index)), -match('gene_affected', colnames(sample.gains))]
			}

		if (engine == 'bedr') {
			overlap <- element.regions.ordered[in.region(element.regions.ordered, regions.ordered)]
			elements.overlap <- elements[match(overlap, element.regions), ]
			gh.overlap <- sample.gains[in.region(regions.ordered, element.regions.ordered), ]
			}
		}
	}
save(
	list = c('elements.overlap', 'gh.overlap'),
	file = generate.filename('genehancer_overlap_regions_gains', 'BRCA_EU', 'rda')
	)
