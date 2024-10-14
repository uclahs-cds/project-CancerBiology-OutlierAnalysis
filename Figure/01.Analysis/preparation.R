library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());


# Get the FDR information of only outlier genes
outlier.gene.fdr.01 <- list(
    brca = outlier.gene.fdr.all.brca[outlier.gene.fdr.all.brca$fdr < 0.01, ],
    meta = outlier.gene.fdr.all.meta[outlier.gene.fdr.all.meta$fdr < 0.01, ],
    matador = outlier.gene.fdr.all.matador[outlier.gene.fdr.all.matador$fdr < 0.01, ],
    ispy = outlier.gene.fdr.all.ispy[outlier.gene.fdr.all.ispy$fdr < 0.01, ],
    icgc = outlier.gene.fdr.all.icgc[outlier.gene.fdr.all.icgc$fdr < 0.01, ]
    );

pos <- which(strsplit(rownames(outlier.gene.fdr.01$matador), '')[[1]] == '_');

outlier.symbol <- list(
    metabric = fpkm.tumor.symbol.filter.meta.symbol[rownames(outlier.gene.fdr.01$meta), ]$Symbol,
    brca = fpkm.tumor.symbol.filter.brca[rownames(outlier.gene.fdr.01$brca), ]$Symbol,
    matador = substring(rownames(outlier.gene.fdr.01$matador), pos + 1),
    ispy = rownames(outlier.gene.fdr.01$ispy),
    icgc = fpkm.tumor.symbol.filter.symbol.icgc[rownames(outlier.patient.tag.01.icgc), ]$Symbol
    );

outlier.symbol$unique <- na.omit(unique(c(
    outlier.symbol$metabric,
    outlier.symbol$brca,
    outlier.symbol$matador,
    outlier.symbol$ispy,
    outlier.symbol$icgc
    )));

### 1. MATADOR
# Find the position of "_" in Metador data row names
outlier.patient.tag.01.metador.pos <- which(
    strsplit(rownames(outlier.patient.tag.01.metador), '')[[1]] == '_'
    );
# Extract the substring after "_"
outlier.patient.tag.01.metador.symbol <- substring(
    rownames(outlier.patient.tag.01.metador),
    outlier.patient.tag.01.metador.pos + 1
    );
# Match and extract rows from Metador data based on unique symbols
outlier.patient.tag.01.metador.match.five <- outlier.patient.tag.01.metador[
    match(outlier.symbol$unique, outlier.patient.tag.01.metador.symbol),
    ];

### 2. TCGA-BRCA
outlier.patient.tag.01.brca.symbol <- fpkm.tumor.symbol.filter.brca[
    rownames(outlier.patient.tag.01.brca),
    ]$Symbol;
outlier.patient.tag.01.brca.match.five <- outlier.patient.tag.01.brca[
    match(outlier.symbol$unique, outlier.patient.tag.01.brca.symbol),
    ];


### 3. METABRIC
outlier.patient.tag.01.meta.symbol <- fpkm.tumor.symbol.filter.meta.symbol[
    rownames(outlier.patient.tag.01.meta),
    ]$Symbol;
outlier.patient.tag.01.meta.match.five <- outlier.patient.tag.01.meta[
    match(outlier.symbol$unique, outlier.patient.tag.01.meta.symbol),
    ];


### 4. ISPY-2
outlier.patient.tag.01.ispy.symbol <- rownames(outlier.patient.tag.01.ispy);
outlier.patient.tag.01.ispy.match.five <- outlier.patient.tag.01.ispy[
    match(outlier.symbol$unique, outlier.patient.tag.01.ispy.symbol),
    ];


### 5. ICGC BRCA-EU
outlier.patient.tag.01.icgc.symbol <- fpkm.tumor.symbol.filter.symbol.icgc[
    rownames(outlier.patient.tag.01.icgc),
    ]$Symbol;
# Match and extract rows from ICGC data based on unique symbols
outlier.patient.tag.01.icgc.match.five <- outlier.patient.tag.01.icgc[
    match(outlier.symbol$unique, outlier.patient.tag.01.icgc.symbol),
    ];

# Combine matched data from all sources into a single data frame
outlier.patient.all.five.01 <- data.frame(
    cbind(
        outlier.patient.tag.01.brca.match.five,
        outlier.patient.tag.01.meta.match.five,
        outlier.patient.tag.01.ispy.match.five,
        outlier.patient.tag.01.metador.match.five,
        outlier.patient.tag.01.icgc.match.five
        )
    );

# Set row names of the combined data frame
rownames(outlier.patient.all.five.01) <- outlier.symbol$unique;

# Save these variables for later scripts
cache.multiple.computed.variables(c(
    'outlier.symbol',
    'outlier.gene.fdr.01',
    'outlier.patient.all.five.01'
    ));

save.session.profile(file.path('output', 'preparation.txt'));
