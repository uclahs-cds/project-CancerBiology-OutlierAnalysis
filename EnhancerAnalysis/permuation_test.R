### permutation_test #########################################################################
# Description
# per patient that has an outlier gene, overlap SV segment with GH db

### HISTORY ##################################################################################
# Version	Date		Developer	Comments
# 0.01		2023-12-21	jlivingstone	initial code
# 0.02		2024-03-12	jlivingstone	update to use inversions and rearrangements

### PREAMBLE #################################################################################
library(bedr)
library(BoutrosLab.utilities)

setwd('/hot/user/jlivingstone/outlier/enhancer_analysis/permutation')

gh <- read.delim(
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19',
		'GeneHancer_AnnotSV_hg19_elements_v5.18_elite.txt'
		),
	as.is = TRUE
	)

path <- '/hot/ref/cohort/ICGC/BRCA/EU/processed/wgs/'

# get sample names & mutation data
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

inversion$chromosome <- inversion$chr_from
inversion$chromosome_start <- inversion$chr_from_bkpt
inversion$chromosome_end <- inversion$chr_to_bkpt

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

break_one <- break_two <- translocation
break_one$chromosome <- translocation$chr_from
break_one$chromosome_start <- translocation$chr_from_bkpt - 1000
break_one$chromosome_end <- translocation$chr_from_bkpt + 1000

break_two$chromosome <- translocation$chr_to
break_two$chromosome_start <- translocation$chr_to_bkpt - 1000
break_two$chromosome_end <- translocation$chr_to_bkpt + 1000

combined.muts <- rbind(
	break_one,
	break_two,
	inversion
	)
muts <- combined.muts[order(combined.muts$submitted_sample_id), ]

samples.with.muts <- union(unique(translocation$submitted_sample_id), unique(inversion$submitted_sample_id))

# read in outlier genes per patient - get number of outlier genes that overlaps with gh & genes to remove
outliers <- read.delim(
	file = '/hot/user/jlivingstone/outlier/run_method/2023-12-21_Outlier_patients_with_genes_BRCA_EU_cutoff_0.01.txt',
	as.is = TRUE
	)

# ugh the WGS names aren't the same as RNA
outliers$sample.name <- sub('R', 'D', outliers$patient)
outliers$sample.name <- sub('.RNA', '', outliers$sample.name)
outliers$sample.name[grep('PD6418a.2', outliers$sample.name)] <- 'PD6418a'

outlier.genes <- unlist(strsplit(outliers$outlier_genes, ';'))

exprs <- read.delim(
	file = '/hot/user/jlivingstone/outlier/NikZainal_2016/original/SupplementaryTable7Transcriptomic342.txt',
	as.is = TRUE
	)

overlap.genes <- intersect(exprs$Name, gh$symbol)

# not all outlier genes are in genehancer so take overlap
n.outlier.genes <- length(intersect(overlap.genes, outlier.genes))

# and remove outlier genes from pool
ind.to.remove <- match(intersect(outlier.genes, overlap.genes), overlap.genes)
overlap.genes <- overlap.genes[-ind.to.remove]

muts.sample <- unique(muts$submitted_sample_id)
samples <- intersect(samples.with.muts, outliers$sample.name)

results <- data.frame(
	seed = numeric(),
	percent = numeric(),
       	n.overlap = numeric(),
	stringsAsFactors = FALSE
       	)

seeds <- 1:1000
for (i in 1:length(seeds)) {
	print(seeds[i])
	set.seed(seeds[i])

	random.outliers <- data.frame(
		patient = character(),
		outlier_genes = character(),
		stringsAsFactors = FALSE
		)
	picked.samples <- sample(samples, n.outlier.genes, replace = TRUE)
	picked.genes <- sample(overlap.genes, n.outlier.genes, replace = TRUE)

	random.outliers <- data.frame(
		patient = picked.samples,
		outlier_genes = picked.genes,
		stringsAsFactors = FALSE
		)

	elements.overlap <- list()
	gh.overlap <- list()

	# for each patient, overlap the enhancer regions with sv breakpoints
	for (j in 1:nrow(random.outliers)) {
		elements <- gh[gh$symbol %in% random.outliers$outlier_genes[j], ]

		if (nrow(elements) > 0) {
			# GH regions
			element.regions <- paste0('chr', elements$chr, ':', elements$element_start, '-', elements$element_end)
			ind <- order(elements$chr, elements$element_start, elements$element_end)
			element.regions.ordered <- element.regions[ind]
	
			sample.muts <- muts[muts$submitted_sample_id %in% random.outliers$patient[j], ]

			# check if mutation regions overlap with genehancer elements (object a are in object b)
			regions <- paste0(sample.muts$chromosome, ':', sample.muts$chromosome_start, '-', sample.muts$chromosome_end)
			s.ind <- order(sample.muts$chromosome, sample.muts$chromosome_start, sample.muts$chromosome_end)
			regions.ordered <- paste0('chr', regions[s.ind])
			sample.muts.ordered <- sample.muts[s.ind,]

			overlap <- element.regions.ordered[in.region(element.regions.ordered, regions.ordered)]
			elements.overlap[[random.outliers$patient[j]]] <- elements[match(overlap, element.regions), ]
			gh.overlap[[random.outliers$patient[j]]] <- sample.muts.ordered[in.region(regions.ordered, element.regions.ordered), ]
			}
		}

		# save results per seed
		n.overlap <- as.numeric(table(unlist(lapply(gh.overlap, function(x) { nrow(x) })) > 0)['TRUE'])
	       	results[i, 'seed'] <- seeds[i]
		results[i, 'n.overlap'] <- n.overlap
		results[i, 'percent'] <- (n.overlap / n.outlier.genes) * 100
	}
