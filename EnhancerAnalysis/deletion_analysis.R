## deletion_analysis.R #################################################################
# Description
# from = start; to = end;

### HISTORY ############################################################################
# Version	Date		Developer	Comments
# 0.01		2024-04-02	jlivingstone	initial code

### PREAMBLE ###########################################################################
library(BoutrosLab.utilities)
library(GenomicRanges)

setwd('/hot/user/jlivingstone/outlier/enhancer_analysis/deletion_analysis')

bp.threshold <- 10000

annot <- read.delim(
	file = file.path(
		'/hot/ref/database/RefSeq/release_2020-01-10',
		'hg19_refGene_annotation.txt'
		),
	as.is = TRUE
	)

gh <- read.delim(
	file = file.path(
		'/hot/ref/database/GeneHancer-v5.18/processed/hg19',
		'GeneHancer_AnnotSV_hg19_elements_v5.18_elite.txt'
		),
	as.is = TRUE
	)

deletion <- read.delim(
	file = file.path(
		'/hot/ref/cohort/ICGC/BRCA/EU/processed/wgs/',
		'deletion.BRCA-EU.tsv'
		),
	as.is = TRUE
	)
# to account for trailing numbers
deletion$submitted_sample_id <- sub('a2', 'a', deletion$submitted_sample_id)
deletion$submitted_sample_id <- sub('a3', 'a', deletion$submitted_sample_id)

samples.with.muts <- unique(deletion$submitted_sample_id)

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

# gh elements for outlier gene
outlier.genes <- unlist(strsplit(outliers.parsed$outlier_genes, ';'))
overlap.genes <- intersect(annot$gene, intersect(outlier.genes, unique(gh$symbol)))

# remove genes that aren't in GH or don't have annotation (strand) for
true.outliers <- data.frame()
for (i in 1:nrow(outliers.parsed)) {
        temp <- data.frame(
                sample.name = outliers.parsed$sample.name[i],
                outlier_genes = unlist(strsplit(outliers.parsed$outlier_genes[i], ';')),
                stringsAsFactors = FALSE
                )
        true.outliers <- rbind(true.outliers, temp)
        }
true.outliers <- true.outliers[match(overlap.genes,true.outliers$outlier_genes),]

results <- list()
for (j in 1:nrow(true.outliers)) {
	gene <- true.outliers$outlier_genes[j]
	print(gene)

	gene.annot <- annot[match(gene, annot$gene),]

	# sample deletions
	muts <- deletion[deletion$submitted_sample_id %in% true.outliers$sample.name[j],]
	chr.muts <- muts[which(muts$chr_from == gene.annot$chr),]

	if (nrow(chr.muts) == 0) {
		results[[j]] <- 'No mutations on the same chr as gene'
		next
		}

	# check if there are breakpoints within X bp of gene start
	if (gene.annot$strand == '+') {
		gene.annot$overlap.start <- gene.annot$start - bp.threshold
		gene.annot$overlap.end <- gene.annot$start

		# upstream is to the left, so need 'end' bkpt of deletion to be in window
		deletion.ranges <- GRanges(
			seqnames = paste0('chr', chr.muts$chr_to),
			ranges = IRanges(
				start = chr.muts$chr_to_bkpt,
				end = chr.muts$chr_to_bkpt + 1
				)
			)
	} else {
		gene.annot$overlap.start <- gene.annot$end
		gene.annot$overlap.end <- gene.annot$end + bp.threshold

		# upstream is to the right, so need 'start' bkpt of deletion to be in window
		deletion.ranges <- GRanges(
			seqnames = paste0('chr', chr.muts$chr_from),
			ranges = IRanges(
			start = chr.muts$chr_from_bkpt,
				end = chr.muts$chr_from_bkpt + 1
				)
			)
		}

	# overlap gene annotation and mutations on same chromosome
	gene.ranges <- GRanges(
		seqnames = paste0('chr', gene.annot$chr),
		strand = gene.annot$strand,
		ranges = IRanges(
			start = gene.annot$overlap.start,
			end = gene.annot$overlap.end
			)
		)

	overlap <- findOverlaps(deletion.ranges, gene.ranges, type = 'within')
	if (length(overlap) == 0) {
		results[[j]] <- 'No overlap with deletion'
		next
		}

	# see if GH elements are upstream
	gh.filtered <- gh[which(gh$chr == gene.annot$chr),]
	gh.ranges <- GRanges(
		seqnames = paste0('chr', gh.filtered$chr),
		ranges = IRanges(
			start = gh.filtered$element_start,
			end = gh.filtered$element_end
			)
		)

	# deletion event that brings elements together
	candidate <- chr.muts[subjectHits(overlap), ]

 	if (gene.annot$strand == '+') {
		candidate$overlap.start <- candidate$chr_from_bkpt - bp.threshold
		candidate$overlap.end <- candidate$chr_from_bkpt
	} else {
		candidate$overlap.start <- candidate$chr_to_bkpt
		candidate$overlap.end <- candidate$chr_to_bkpt + bp.threshold
		}

	candidate.ranges <- GRanges(
		seqnames = paste0('chr', candidate$chr_from),
		ranges = IRanges(
			start = candidate$overlap.start,
			end = candidate$overlap.end
			)
		)

	#gh.overlap <- findOverlaps(gh.ranges, candidate.ranges)
	closest <- data.frame(
		distanceToNearest(gh.ranges, candidate.ranges),
		stringsAsFactors = FALSE
		)
	gh.overlap <- closest$subjectHits[which(closest$distance < bp.threshold)]
	if (length(gh.overlap) == 0) {
		results[[j]] <- 'GH element not in range'
	} else {
		results[[j]] <- cbind(candidate[gh.overlap$subjectHits, ], gh.filtered[gh.overlap$queryHits,])
		}
	}
names(results) <- paste(true.outliers$sample.name, true.outliers$outlier_genes, sep = '-')

save(
	results,
	file = generate.filename('outlier', paste0(bp.threshold, 'bp_deletion_analysis'), 'rda')
	)
