### HISTORY #####################################################################
# This script generates a scatter plot to visualize the correlation between 
# protein abundance (CPTAC/RPPA) and mRNA abundance (FPKM) in 
# TCGA-BRCA outlier patients. 
# This code is linked with the analysis in Figure 3a and b.
# Date: 2024-08-14

library(BoutrosLab.plotting.general);


### 1. CPTAC example

# Two inputs: IGF1R, CLU - run separately
do.plot.3d <- function(i, filename) {

# Extract matching protein data for a specific gene (i)
brca.protein.cptac.outlier.match <- brca.protein.cptac[,colnames(brca.protein.cptac) %in% substr(colnames(outlier.patient.tag.01.brca), 1, 15)];
brca.protein.cptac.outlier.match.i <- brca.protein.cptac.outlier.match[i, ];

# Match corresponding FPKM data
fpkm.protein.cptac.match.i <- fpkm.tumor.symbol.filter.brca[
    fpkm.tumor.symbol.filter.brca$Symbol == i, 
    match(colnames(brca.protein.cptac.outlier.match.i), substr(colnames(fpkm.tumor.symbol.filter.brca), 1, 15))
    ];

# Identify outlier patients for the specific gene (i)
outlier.patient.tag.01.brca.protein.cptac.zscore.match.i <- colnames(
    outlier.patient.tag.01.brca.protein.cptac.zscore.match
    )[outlier.patient.tag.01.brca.protein.cptac.zscore.match[
        rownames(fpkm.protein.cptac.match.i), ] == 1
    ];

# Extract FPKM and protein CPTAC data for outlier patients
i.fpkm <- fpkm.protein.cptac.match.i[outlier.patient.tag.01.brca.protein.cptac.zscore.match.i];
i.protein.cptac <- brca.protein.cptac.outlier.match.i[
    substr(outlier.patient.tag.01.brca.protein.cptac.zscore.match.i, 1, 15)
    ];

# Combine protein and mRNA data into a comparison data frame
protein.cptac.rna.i.comparison <- data.frame(
    cbind(
        as.numeric(brca.protein.cptac.outlier.match.i[1, ]),
        as.numeric(fpkm.protein.cptac.match.i)
        )
    );

# Rename columns for clarity
colnames(protein.cptac.rna.i.comparison) <- c('non', 'out');

# Define colors for scatter plot based on outlier status
dot.colours <- rep('black', nrow(protein.cptac.rna.i.comparison));
dot.colours[
    which(colnames(fpkm.protein.cptac.match.i) == outlier.patient.tag.01.brca.protein.cptac.zscore.match.i)
    ] <- 'red2';

cptac.gene.scatter <- create.scatterplot(
    formula = non ~ log2(out),
    data = protein.cptac.rna.i.comparison,
    col = dot.colours,
    alpha = .55,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    add.grid = TRUE,
    grid.colour = 'grey80',
    cex = 0.75,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.5,
    main = i,
    xlab.label = expression('mRNA abundance (log'[2]*'(FPKM))'),
    ylab.label = expression('Protein abundance (z-score)'),
    text.x = log2(as.numeric(i.fpkm)),
    text.y = as.numeric(i.protein.cptac),
    text.labels = '*outlier patient',
    text.guess.labels = TRUE,
    text.guess.label.position = 180,
    text.guess.radius.factor = 1.5,
    text.fontface = 1,
    text.col = 'red2',
    add.text = TRUE,
    type = c('p', 'r', 'g'), 
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = protein.cptac.rna.i.comparison$non,
                    y = protein.cptac.rna.i.comparison$out,
                    label.items = c('spearman'),
                    alpha.background = 0,
                    key.cex = 1.1
                    ) 
                ),
            x = 0.03,
            y = 0.95,
            corner = c(0,1)
            )
        ),
    );

# Save the plot as a PNG
write.plot(
    trellis.object = cptac.gene.scatter,
    filename = filename,
    width = 6,
    height = 6,
    size.units = 'in',
    resolution = 1200
);

}

output.directory <- get0('output.directory', ifnotfound = 'figures');

do.plot.3d('IGF1R', file.path(output.directory, 'Figure_3_d_IGR1R.png'));
do.plot.3d('CLU', file.path(output.directory, 'Figure_3_d_CLU.png'));

### 2. RPPA example
i <- 'NRAS';
# Extract i protein data from the matched dataset
brca.protein.outlier.match.i <- brca.protein.outlier.match[
    rownames(protein.antibody.outlier)[protein.antibody.outlier$gene.name == i], 
    ];

# Match corresponding FPKM data for i
fpkm.protein.match.i <- fpkm.tumor.symbol.filter.brca[
    fpkm.tumor.symbol.filter.brca$Symbol == i, 
    colnames(brca.protein.outlier.match.i)
    ];

# Combine protein and mRNA data into a comparison data frame for i
protein.rna.i.comparison <- data.frame(
    cbind(
        as.numeric(brca.protein.outlier.match.i[1, ]),
        as.numeric(fpkm.protein.match.i)
        )
    );

# Rename columns for clarity
colnames(protein.rna.i.comparison) <- c('non', 'out');

# Identify outlier patients for i
outlier.patient.tag.01.brca.protein.match.no.p.i <- colnames(
    outlier.patient.tag.01.brca.protein.match.no.p
    )[outlier.patient.tag.01.brca.protein.match.no.p[
        rownames(fpkm.protein.match.i), ] == 1
    ];

# Extract FPKM and protein data for outlier patients
i.fpkm <- fpkm.protein.match.i[
    outlier.patient.tag.01.brca.protein.match.no.p.i
    ];

i.protein <- brca.protein.outlier.match.i[
    1, outlier.patient.tag.01.brca.protein.match.no.p.i
    ];

# Define colors for scatter plot based on outlier status
dot.colours <- rep('black', nrow(protein.rna.i.comparison));
dot.colours[
    which(colnames(fpkm.protein.match.i) %in% outlier.patient.tag.01.brca.protein.match.no.p.i)
    ] <- 'red2';


rppa.gene.scatter <- BoutrosLab.plotting.general::create.scatterplot(
    formula = non ~ log2(out),
    data = protein.rna.i.comparison[nrow(protein.rna.i.comparison):1, ],
    col = rev(dot.colours),
    alpha = .55,
    grid.colour = 'grey80',
    cex = 0.75,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.5,
    main = i,
    xlab.label = expression('mRNA abundance (log'[2]*'(FPKM))'),
    ylab.label = expression('Protein abundance (RPPA)'),
    text.x = log2(as.numeric(i.fpkm)),
    text.y = as.numeric(i.protein),
    text.labels = rep('*outlier patient', length(i.protein)),
    text.guess.labels = TRUE,
    text.guess.label.position = 160,
    text.guess.radius.factor = 1.5,
    text.fontface = 1,
    text.col = 'red2',
    add.text = TRUE,
    type = c('p', 'r', 'g'),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = protein.rna.i.comparison$non,
                    y = protein.rna.i.comparison$out,
                    label.items = c('spearman'),
                    alpha.background = 0,
                    key.cex = 1.1
                )
            ),
            x = 0.03,
            y = 0.95,
            corner = c(0, 1)
        )
    )
);

output.directory <- get0('output.directory', ifnotfound = 'figures');

# Save the plot as a PNG
write.plot(
    trellis.object = rppa.gene.scatter,
    filename = file.path(output.directory, 'Figure_3_d_NRAS.png'),
    width = 6,
    height = 6,
    size.units = 'in',
    resolution = 1200
);
