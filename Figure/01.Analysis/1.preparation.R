# Load necessary libraries
library(BoutrosLab.utilities)
library(outlierAnalysisSupport) # Load custom helper library

### DATA PREPARATION ############################################################
# Attach the file path for accessing outlier data
attach(get.outlier.data.path())

# Filter outlier genes with FDR < 0.01 from each dataset (BRCA, METABRIC, etc.)
outlier.gene.fdr.01 <- list(
    brca = outlier.gene.fdr.all.brca[outlier.gene.fdr.all.brca$fdr < 0.01, ],
    meta = outlier.gene.fdr.all.meta[outlier.gene.fdr.all.meta$fdr < 0.01, ],
    matador = outlier.gene.fdr.all.matador[outlier.gene.fdr.all.matador$fdr < 0.01, ],
    ispy = outlier.gene.fdr.all.ispy[outlier.gene.fdr.all.ispy$fdr < 0.01, ],
    icgc = outlier.gene.fdr.all.icgc[outlier.gene.fdr.all.icgc$fdr < 0.01, ]
    )

# Identify the position of underscores in row names of the matador dataset
pos <- which(strsplit(rownames(outlier.gene.fdr.01$matador), '')[[1]] == '_')

# Map outlier gene symbols to their respective datasets
outlier.symbol <- list(
    metabric = fpkm.tumor.symbol.filter.meta.symbol[rownames(outlier.gene.fdr.01$meta), ]$Symbol,
    brca = fpkm.tumor.symbol.filter.brca[rownames(outlier.gene.fdr.01$brca), ]$Symbol,
    matador = substring(rownames(outlier.gene.fdr.01$matador), pos + 1),
    ispy = rownames(outlier.gene.fdr.01$ispy),
    icgc = fpkm.tumor.symbol.filter.symbol.icgc[rownames(outlier.patient.tag.01.icgc), ]$Symbol
    )

# Combine unique symbols across all datasets
outlier.symbol$unique <- na.omit(unique(c(
    outlier.symbol$metabric,
    outlier.symbol$brca,
    outlier.symbol$matador,
    outlier.symbol$ispy,
    outlier.symbol$icgc
    )))

### MATCH OUTLIER PATIENTS ACROSS DATASETS #####################################

# 1. MATADOR
# Extract the symbol part after the underscore for the Metador dataset
outlier.patient.tag.01.metador.pos <- which(strsplit(rownames(outlier.patient.tag.01.metador), '')[[1]] == '_')
outlier.patient.tag.01.metador.symbol <- substring(rownames(outlier.patient.tag.01.metador), outlier.patient.tag.01.metador.pos + 1)

# Match unique symbols to Metador outlier patients
outlier.patient.tag.01.metador.match.five <- outlier.patient.tag.01.metador[
    match(outlier.symbol$unique, outlier.patient.tag.01.metador.symbol),
    ]

# 2. TCGA-BRCA
# Match unique symbols to BRCA outlier patients
outlier.patient.tag.01.brca.symbol <- fpkm.tumor.symbol.filter.brca[rownames(outlier.patient.tag.01.brca), ]$Symbol
outlier.patient.tag.01.brca.match.five <- outlier.patient.tag.01.brca[
    match(outlier.symbol$unique, outlier.patient.tag.01.brca.symbol),
    ]

# 3. METABRIC
# Match unique symbols to METABRIC outlier patients
outlier.patient.tag.01.meta.symbol <- fpkm.tumor.symbol.filter.meta.symbol[rownames(outlier.patient.tag.01.meta), ]$Symbol
outlier.patient.tag.01.meta.match.five <- outlier.patient.tag.01.meta[
    match(outlier.symbol$unique, outlier.patient.tag.01.meta.symbol),
    ]

# 4. ISPY-2
# Match unique symbols to ISPY-2 outlier patients
outlier.patient.tag.01.ispy.symbol <- rownames(outlier.patient.tag.01.ispy)
outlier.patient.tag.01.ispy.match.five <- outlier.patient.tag.01.ispy[
    match(outlier.symbol$unique, outlier.patient.tag.01.ispy.symbol),
    ]

# 5. ICGC BRCA-EU
# Match unique symbols to ICGC outlier patients
outlier.patient.tag.01.icgc.symbol <- fpkm.tumor.symbol.filter.symbol.icgc[rownames(outlier.patient.tag.01.icgc), ]$Symbol
outlier.patient.tag.01.icgc.match.five <- outlier.patient.tag.01.icgc[
    match(outlier.symbol$unique, outlier.patient.tag.01.icgc.symbol),
    ]

### COMBINE DATA ###############################################################
# Combine matched patient data from all sources into a single data frame
outlier.patient.all.five.01 <- data.frame(
    cbind(
        outlier.patient.tag.01.brca.match.five,
        outlier.patient.tag.01.meta.match.five,
        outlier.patient.tag.01.ispy.match.five,
        outlier.patient.tag.01.metador.match.five,
        outlier.patient.tag.01.icgc.match.five
        )
    )

# Set row names of the combined data frame to unique outlier symbols
rownames(outlier.patient.all.five.01) <- outlier.symbol$unique

### SAVE VARIABLES #############################################################
# Cache the important variables for later use
cache.multiple.computed.variables(c(
    'outlier.symbol',
    'outlier.gene.fdr.01',
    'outlier.patient.all.five.01'
    ))

# Save the session profile for reproducibility
save.session.profile(file.path('output', '1.preparation.txt'))
