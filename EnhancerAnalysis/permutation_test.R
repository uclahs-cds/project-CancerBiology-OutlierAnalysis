### permutation_test #########################################################################
# Description
# shuffle gene and patient to see if there is overlap of mutation (ins, trans) with gh element

# Rscript permutation_test.R --seed.start 1 --seed.end 25 --elite 1
### HISTORY ##################################################################################
# Version	Date		Developer	Comments
# 0.01		2023-12-21	jlivingstone	initial code
# 0.02		2024-03-12	jlivingstone	update to use inversions and rearrangements

### PREAMBLE #################################################################################
library(bedr)
library(BoutrosLab.utilities)
library(getopt)

setwd('/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/enhancer_analysis/permutation')

params <- matrix(
        data = c(
                'seed.start', 's', '0', 'character',
                'seed.end', 'e', '0', 'character',
                'elite', 'v', '0', 'character'
                ),
        ncol = 4,
        byrow = TRUE
        );
opt <- getopt(params);

if (opt$elite == 1) {
	filename <- 'GeneHancer_AnnotSV_hg19_elements_v5.18_elite.txt'
	placeholder <- 'elite_'
} else {
	filename <- 'GeneHancer_AnnotSV_hg19_elements_v5.18.txt'
	placeholder <- ''
	}

gh <- read.delim(
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19',
		filename
		),
	as.is = TRUE
	)

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

# only overlap brekpoints, not whole regions, so we don't have to account for size of mutation; add some padding
ibreak_one <- ibreak_two <- inversion
ibreak_one$chromosome <- inversion$chr_from
ibreak_one$chromosome_start <- inversion$chr_from_bkpt - 1000
ibreak_one$chromosome_end <- inversion$chr_from_bkpt + 1000

ibreak_two$chromosome <- inversion$chr_to
ibreak_two$chromosome_start <- inversion$chr_to_bkpt - 1000
ibreak_two$chromosome_end <- inversion$chr_to_bkpt + 1000

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
	ibreak_one,
	ibreak_two
	)
muts <- combined.muts[order(combined.muts$submitted_sample_id), ]

samples.with.muts <- union(unique(translocation$submitted_sample_id), unique(inversion$submitted_sample_id))

exprs <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone',
		'NikZainal_2016/original/SupplementaryTable7Transcriptomic342.txt'
		),
	as.is = TRUE
	)
colnames(exprs) <- sub('R', 'D', colnames(exprs))
colnames(exprs) <- sub('.RNA', '', colnames(exprs))
colnames(exprs)[grep('PD6418a.2', colnames(exprs))] <- 'PD6418a'

# read in outlier genes per patient - get number of outlier genes that overlaps with gh & genes to remove
outliers <- read.delim(
	file = file.path(
		'/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/jlivingstone/run_method',
		'2023-12-21_Outlier_patients_with_genes_BRCA_EU_cutoff_0.01.txt'
		),
	as.is = TRUE
	)

# ugh the WGS names aren't the same as RNA
outliers$sample.name <- sub('R', 'D', outliers$patient)
#outliers$sample.name <- sub('.RNA', '', outliers$sample.name)
outliers$sample.name[grep('PD6418a.2', outliers$sample.name)] <- 'PD6418a'

# need to filter out samples/genes that dont have mutation data
outlier.overlap.samples <- intersect(samples.with.muts, outliers$sample.name)
outliers.parsed <- outliers[match(outlier.overlap.samples, outliers$sample.name), ]

# run true dataset
true.outliers <- data.frame()
for (i in 1:nrow(outliers.parsed)) {
	temp <- data.frame(
		sample.name = outliers.parsed$sample.name[i],
		outlier_genes = unlist(strsplit(outliers.parsed$outlier_genes[i], ';')),
		stringsAsFactors = FALSE
		)
	true.outliers <- rbind(true.outliers, temp)
	}
overlap.genes <- intersect(exprs$Name, gh$symbol)

outliers.final <- true.outliers[match(intersect(true.outliers$outlier_genes, gh$symbol), true.outliers$outlier_genes),]
outlier.genes <- outliers.final$outlier_genes

#write.table(
#	x = outliers.final,
#	file = generate.filename('outliers', 'final_outliers_overlapped_gh_exprs', 'tsv'),
#	quote = FALSE,
#	row.names = FALSE
#	)

# not all outlier genes are in genehancer so take overlap (n = 174)
n.outlier.genes <- length(intersect(overlap.genes, outlier.genes))

# and remove outlier genes from pool
ind.to.remove <- match(intersect(outlier.genes, overlap.genes), overlap.genes)
overlap.genes <- overlap.genes[-ind.to.remove]

# overlap with samples that have exprs data - can randomly pick samples that dont have outliers
sample.pool <- intersect(samples.with.muts, colnames(exprs))

results <- data.frame(
	seed = numeric(),
	percent = numeric(),
	n.overlap = numeric(),
	stringsAsFactors = FALSE
	)

# run true dataset
#random.outliers <- outliers.final

seeds <- opt$seed.start:opt$seed.end
simulated.set <- list()

for (i in 1:length(seeds)) {
	print(seeds[i])

	random.outliers <- data.frame(
		patient = character(),
		outlier_genes = character(),
		stringsAsFactors = FALSE
		)

	set.seed(seeds[i])
	picked.samples <- sample(sample.pool, n.outlier.genes, replace = TRUE)
	set.seed(seeds[i])
	picked.genes <- sample(overlap.genes, n.outlier.genes, replace = TRUE)

	random.outliers <- data.frame(
		sample.name = picked.samples,
		outlier_genes = picked.genes,
		stringsAsFactors = FALSE
		)

	simulated.set[[i]] <- random.outliers
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
	
			sample.muts <- muts[muts$submitted_sample_id %in% random.outliers$sample.name[j], ]

			# check if mutation regions overlap with genehancer elements (object a are in object b)
			regions <- paste0(sample.muts$chromosome, ':', sample.muts$chromosome_start, '-', sample.muts$chromosome_end)
			s.ind <- order(sample.muts$chromosome, sample.muts$chromosome_start, sample.muts$chromosome_end)
			regions.ordered <- paste0('chr', regions[s.ind])
			sample.muts.ordered <- sample.muts[s.ind,]

			overlap <- element.regions.ordered[in.region(element.regions.ordered, regions.ordered)]
			elements.overlap[[j]] <- elements[match(overlap, element.regions), ]
			gh.overlap[[j]] <- sample.muts.ordered[in.region(regions.ordered, element.regions.ordered), ]
			}
		}

		# save results per seed
		n.overlap <- as.numeric(table(unlist(lapply(gh.overlap, function(x) { nrow(x) })) > 0)['TRUE'])
		results[i, 'seed'] <- seeds[i]
		results[i, 'n.overlap'] <- n.overlap
		results[i, 'percent'] <- (n.overlap / n.outlier.genes) * 100
	}

save(
	simulated.set,
	file = generate.filename(
		'outlier',
		paste0('seeds_', opt$seed.start, '_', opt$seed.end, '_simulated_sets', placeholder, 'bkpt_padding_1000_overlap'),
		'rda'
		)
	)

write.table(
	x = results,
	file = generate.filename(
		'outlier',
		paste0('seeds_', opt$seed.start, '_', opt$seed.end, '_permutation_results_', placeholder, 'bkpt_padding_1000_overlap'),
		'txt'
		),
	row.names = FALSE,
	quote = FALSE,
	sep = '\t'
	)
