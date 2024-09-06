### HISTORY #####################################################################
# This script analyzes the protein abundance of outlier genes identified in CCLE 
# overlapped with tissue datasets.
# Date: 2024-08-16
#################################################################################

library(BoutrosLab.plotting.general);

# Haven't upload these files on cluster. These are included as variables
# Get protein info
# protein.info.pancancer <- read.delim2(
#     file = '1-s2.0-S1535610822002744-mmc3_name.csv',
#     header = T, row.names = 1, sep = ','
#     );

protein.info.pancancer.breast <- protein.info.pancancer[
    protein.info.pancancer$BROAD_ID %in% colnames(fpkm.tumor.symbol.filter.ccle),
    ];

protein.info.pancancer.breast.match <- protein.info.pancancer.breast[
    , 3:ncol(protein.info.pancancer.breast)
    ];


library(stringr);

colnames(protein.info.pancancer.breast.match) <- str_extract(
    colnames(protein.info.pancancer.breast.match), 
    "(?<=\\.).+(?=_)"
    );
rownames(protein.info.pancancer.breast.match) <- protein.info.pancancer.breast$BROAD_ID;

protein.info.pancancer.breast.match.t <- t(protein.info.pancancer.breast.match);

ccle.sample.outlier.status.only <- ccle.sample.outlier.status[rownames(ccle.outlier.rank.fdr.05),];
ccle.sample.outlier.status.only.five <- ccle.sample.outlier.status.only[sub("\\..*", "", rownames(ccle.sample.outlier.status.only)) %in% five.data.outlier.symbol,];

# Tissue overall outlier status
ccle.sample.outlier.status.only.five.symbol <- sub("\\..*", "", rownames(ccle.sample.outlier.status.only.five));

protein.info.pancancer.breast.match.05 <- protein.info.pancancer.breast.match.t[
    rownames(protein.info.pancancer.breast.match.t) %in% ccle.sample.outlier.status.only.five.symbol, 
    ];

# Only outlier gene's FPKM
fpkm.tumor.symbol.filter.ccle.outlier.pancancer <- fpkm.tumor.symbol.filter.ccle[
    sub("\\..*", "", rownames(fpkm.tumor.symbol.filter.ccle)) %in% rownames(protein.info.pancancer.breast.match.05),
    colnames(protein.info.pancancer.breast.match.05)
    ];

# Only outlier gene's status
ccle.sample.outlier.status.protein.match.pancancer <- ccle.sample.outlier.status[
    rownames(fpkm.tumor.symbol.filter.ccle.outlier.pancancer), 
    colnames(protein.info.pancancer.breast.match.05)
    ];

rownames(fpkm.tumor.symbol.filter.ccle.outlier.pancancer) <- sub(
    "\\..*", "", 
    rownames(fpkm.tumor.symbol.filter.ccle.outlier.pancancer)
    );
rownames(ccle.sample.outlier.status.protein.match.pancancer) <- sub(
    "\\..*", "", 
    rownames(ccle.sample.outlier.status.protein.match.pancancer)
    );

# Only overlapped outlier genes with tissue
outlier.protein.ccle.zscore.list.pancancer <- list();
non.outlier.protein.ccle.zscore.list.pancancer <- list();
target.gene.ccle.zscore.list.pancancer <- NULL;

for (i in 1:nrow(protein.info.pancancer.breast.match.05)) {
    # Protein target name
    target.gene.name.protein <- rownames(protein.info.pancancer.breast.match.05)[i];
    
    # Outlier patient column
    target.col <- colnames(ccle.sample.outlier.status.protein.match.pancancer)[
        ccle.sample.outlier.status.protein.match.pancancer[target.gene.name.protein,] == 1
        ];
    non.target.col <- colnames(ccle.sample.outlier.status.protein.match.pancancer)[
        ccle.sample.outlier.status.protein.match.pancancer[target.gene.name.protein,] == 0
        ];
    
    target.gene.ccle.zscore.list.pancancer <- c(
        target.gene.ccle.zscore.list.pancancer, 
        target.gene.name.protein
        );
    
    outlier.protein.ccle.zscore.list.pancancer[[i]] <- protein.info.pancancer.breast.match.05[
        i, target.col
        ];
    
    non.outlier.protein.ccle.zscore.list.pancancer[[i]] <- protein.info.pancancer.breast.match.05[
        i, non.target.col
        ];
    }

outlier.protein.ccle.zscore.value.pancancer <- as.numeric(
    unlist(outlier.protein.ccle.zscore.list.pancancer)
    );
non.outlier.protein.ccle.zscore.value.pancancer <- as.numeric(
    unlist(non.outlier.protein.ccle.zscore.list.pancancer)
    );

# Box plot - compare the values between patients
# - excluding the genes having no outlier patient info

names(outlier.protein.ccle.zscore.list.pancancer) <- target.gene.ccle.zscore.list.pancancer;
outlier.protein.ccle.list.no.p.na.pancancer <- na.omit(
    unlist(outlier.protein.ccle.zscore.list.pancancer)
    );
names(non.outlier.protein.ccle.zscore.list.pancancer) <- target.gene.ccle.zscore.list.pancancer;
non.outlier.protein.ccle.list.no.p.na.pancancer <- non.outlier.protein.ccle.zscore.list.pancancer[
    names(outlier.protein.ccle.list.no.p.na.pancancer)
    ];

protein.ccle.na.value.pancancer <- data.frame(
    protein.ccle.na.value = c(
        as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na.pancancer)), 
        as.numeric(unlist(outlier.protein.ccle.list.no.p.na.pancancer))
        )
    );
protein.ccle.na.value.box.pancancer <- data.frame(
    cbind(
        protein.ccle.na.value.pancancer$protein.ccle.na.value, 
        c(
            rep('non', length(as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na.pancancer)))), 
            rep('out', length(as.numeric(unlist(outlier.protein.ccle.list.no.p.na.pancancer))))
            )
        )
    );
colnames(protein.ccle.na.value.box.pancancer) <- c('protein.ccle.value', 'status');
protein.ccle.na.value.box.pancancer[,1] <- as.numeric(protein.ccle.na.value.box.pancancer[,1]);

wilcox.result.protien.na.pancancer <- wilcox.test(
    as.numeric(unlist(outlier.protein.ccle.list.no.p.na.pancancer)),
    as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na.pancancer)),
    alternative = "two.sided", conf.int = TRUE
    );

text.pvalue.protien.na.pancancer <- display.statistical.result(
    x = wilcox.result.protien.na.pancancer$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

key.protien.na.pancancer <- list(
    text = list(
        lab = text.pvalue.protien.na.pancancer, 
        cex = 1
        ),
    x = 0.25,
    y = 0.95
    );


ccle.protein.box <- BoutrosLab.plotting.general::create.boxplot(
    formula = protein.ccle.value ~ status,
    data = protein.ccle.na.value.box.pancancer,
    main = expression('Protein abundance of outlier genes'),
    main.cex = 1.3,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Protein abundance'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.lab = c('Non-outlier Samples', 'Outlier Samples'),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.rot = 90,
    outliers = FALSE,
    key = key.protien.na.pancancer,
    ylimits = c(-3, 12.5),
    add.stripplot = TRUE,
    points.pch = 1,
    points.cex = 0.8,
    points.col = 'grey60',
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 4),
    xright.rectangle =c(4, 5),
    ybottom.rectangle = -6,
    ytop.rectangle = 15,
    col.rectangle = "grey",
    alpha.rectangle = 0.25,
    lwd = 1.2,
    col = c('red2', 'dodgerblue3'),
    alpha = 0.3
    );


# Save the plot as a PNG
output.directory <- get0('output.directory', ifnotfound = 'figures');

# Save the box plot as a PNG
write.plot(
    trellis.object = ccle.protein.box,
    filename = file.path(output.directory, 'Figure_4_c.png'),
    width = 3.5,
    height = 6.5,
    size.units = 'in',
    resolution = 1200
);
