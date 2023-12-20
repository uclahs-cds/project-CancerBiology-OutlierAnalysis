### Usage ####
# Rscript 1.Outlier_Detection.R --dataset.name BRCA_EU --working.directory /hot/users/jlivingstone/outlier/run_method \
# --outlier.rank.file '/hot/users/jlivingstone/outlier/run_method/2023-11-20_BRCA-EU_final_outlier_rank_bic.short.rda' \
# --qvalue.cutoff 0.1

### 9.Patient_Outlier_Identification.R ####################################################
library(BoutrosLab.utilities)
library(getopt)

# add cutoff pvalue as parameter
params <- matrix(
        data = c(
                'dataset.name', 'd', '0', 'character',
                'working.directory', 'w', '0', 'character',
                'outlier.rank.file', 'r', '0', 'character',
                'qvalue.cutoff', 'p', '0', 'numeric'
                ),
        ncol = 4,
        byrow = TRUE
        );

opt <- getopt(params);
dataset.name <- opt$dataset.name
working.directory <- opt$working.directory
outlier.rank.file <- opt$outlier.rank.file
qvalue.cutoff <- opt$pvalue.cutoff

setwd(working.directory)

files <- list.files(
	pattern = 'Significant_Outlier_Pvalue_Calculation',
	recursive = TRUE,
	full.names = TRUE
	)

# order files based on patients.to.remove indicator
chars.to.remove <- nchar('.rda')
patients.removed <- rep(x = NA, length = length(files))
for (i in 1:length(files)) {
	patients.removed[i] <- as.numeric(substr(files[i], nchar(files[i]) - chars.to.remove, nchar(files[i]) - chars.to.remove))
	}
files.ordered <- files[order(patients.removed)]

load(
	file = outlier.rank.file
	)
fpkm.tumor.symbol.filter <- data.matrix(fpkm.tumor.symbol.filter)

outlier.patient.list <- list()
for (i in 1:length(files.ordered)) {
	print(files.ordered[i])
	variable <- load(files.ordered[i])

	pvalue <- get(variable)
	pvalue$qvalue <- p.adjust(
		p = pvalue$obs.p.value,
		method = 'fdr'
		)

	patients.to.remove <- as.numeric(
		x = substr(
			files.ordered[i],
			start = nchar(files.ordered[i]) - chars.to.remove,
			stop = nchar(files.ordered[i]) - chars.to.remove
			)
		)

	outlier.genes <- pvalue[which(pvalue$qvalue < qvalue.cutoff),]
	annot.filtered <- annot.filter[match(outlier.genes$gene, rownames(annot.filter)),]
	outlier.matrix <- fpkm.tumor.symbol.filter[outlier.genes$gene,]

	#for each gene get the patient with the highest (outlier) value
	outlier.patient <- data.frame(
		gene.index = character(),
		gene.name = character(),
		ensembl = character(),
		patient = character(),
		value = numeric(),
		median = numeric(),
		mean = numeric(),
		stringsAsFactors = FALSE
		)

	for (j in 1:nrow(outlier.matrix)) {
		# need to remove outlier patients based on method iteration
		values <- sort(outlier.matrix[j, ])[1:(ncol(outlier.matrix) - patients.to.remove)]

		outlier.patient[j, 'gene.index'] <- rownames(outlier.matrix)[j]
		outlier.patient[j, 'gene.name'] <- annot.filtered$Name[j]
		outlier.patient[j, 'ensembl'] <- annot.filtered$Ensembl[j]
		outlier.patient[j, 'patient'] <- names(which.max(values))
		outlier.patient[j, 'value'] <- values[which.max(values)]
		outlier.patient[j, 'median'] <- median(values[-which.max(values)])
		outlier.patient[j, 'mean'] <- mean(values[-which.max(values)])
		}
	outlier.patient.list[[i]] <- outlier.patient
	}

save(
	x = outlier.patient.list,
	file = generate.filename('Outlier_patient_identification', dataset.name, 'rda')
	)

# create file per gene
gene.list <- outlier.patient.list[[1]][,c('gene.index', 'gene.name', 'patient')]
for (j in 2:length(outlier.patient.list)) {
	for (i in 1:nrow(gene.list)) {
		matched.index <- match(gene.list$gene.index[i], outlier.patient.list[[j]]$gene.index)
		if (!is.na(matched.index)) {
			outlier.patient.to.add <- outlier.patient.list[[j]][matched.index, 'patient']
			gene.list$patient[i] <- paste(gene.list$patient[i], outlier.patient.to.add, sep = ';')
			}
		}
	}

# and file per sample
outlier.patient.name <- unique(
	unlist(
		lapply(
			X = outlier.patient.list,
			FUN = function(x) {
				x$patient
				}
			)
		)
	)
outlier.patient.genes <- data.frame(
	patient = character(),
	outlier_genes = character(),
	stringsAsFactors = FALSE
	)
for (j in 1:length(outlier.patient.name)) {
	patient.list <- list()
	for (i in 1:length(outlier.patient.list)) {
		patient.list[[i]] <- outlier.patient.list[[i]]$gene.name[which(outlier.patient.list[[i]]$patient %in% outlier.patient.name[j])]
		outlier.patient.genes[j, 'patient'] <- outlier.patient.name[j]
		outlier.patient.genes[j, 'outlier_genes'] <- paste(unique(unlist(patient.list)), collapse = ';')
		}
	}

write.table(
	x = outlier.patient.genes,
	file = generate.filename('Outlier_patient_genes', dataset.name, 'txt'),
	quote = FALSE,
	sep = '\t',
	row.names = FALSE
	)
