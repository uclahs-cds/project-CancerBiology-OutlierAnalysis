### parse_brca_eu.R #############################################################
# Description
# copy_number_somatic_mutation.BRCA-EU.tsv.gz has duplicate lines per segment

### HISTORY ###################################################################
# Version	Date		Developer	Comments
# 0.01		2023-01-07	jlivingstone	initial code

### PREAMBLE ##################################################################

gain <- read.delim(
	file = '/hot/ref/cohort/ICGC/BRCA/EU/processed/gain_somatic_mutation.BRCA-EU.tsv',
	as.is = TRUE
	)

gain$segment <- paste0(gain$submitted_sample_id, '-', gain$chromosome, ':', gain$chromosome_start, '-', gain$chromosome_end)
ids <- unique(gain$segment)

genes.affected <- matrix(
	data = NA,
	nrow = length(ids),
	ncol = 1
	)
for (i in 1:length(ids)) {
	if (i %% 1000 == 0) {
		print(i)
		}
	genes.affected[i] <- paste(gain$gene_affected[which(gain$segment == ids[i])], collapse = ';')
	}

gain.filtered <- gain[match(ids, gain$segment), ]
gain.filtered$gene_affected <- genes.affected

write.table(
	x = gain.filtered,
	file = 'gain_unique_somatic_mutation.BRCA-EU.tsv',
	quote = FALSE,
	row.names = FALSE
	)

loss <- read.delim(
	file = '/hot/ref/cohort/ICGC/BRCA/EU/processed/loss_somatic_mutation.BRCA-EU.tsv',
	as.is = TRUE
	)

loss$segment <- paste0(loss$submitted_sample_id, '-', loss$chromosome, ':', loss$chromosome_start, '-', loss$chromosome_end)
loss.ids <- unique(loss$segment)

loss.genes.affected <- matrix(
	data = NA,
	nrow = length(loss.ids),
	ncol = 1
	)
for (i in 1:length(loss.ids)) {
	if (i %% 1000 == 0) {
		print(i)
		}
	loss.genes.affected[i] <- paste(loss$gene_affected[which(loss$segment == loss.ids[i])], collapse = ';')
	}

loss.filtered <- loss[match(loss.ids, loss$segment), ]
loss.filtered$gene_affected <- loss.genes.affected

write.table(
	x = loss.filtered,
	file = '/hot/ref/cohort/ICGC/BRCA/EU/processed/loss_unique_somatic_mutation.BRCA-EU.tsv',
	quote = FALSE,
	row.names = FALSE
	)
