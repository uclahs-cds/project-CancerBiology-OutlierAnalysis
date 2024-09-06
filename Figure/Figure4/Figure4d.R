### HISTORY #####################################################################
# This script analyzes the protein abundance of all outlier genes identified in CCLE. 
# Date: 2024-08-16

library(BoutrosLab.plotting.general);

# All CCLE outlier genes
protein.info.pancancer.breast.match.all <- protein.info.pancancer.breast.match.t[
    rownames(protein.info.pancancer.breast.match.t) %in% sub("\\..*", "", rownames(ccle.sample.outlier.status.only)), 
    ];


fpkm.tumor.symbol.filter.ccle.outlier.pancancer.all <- fpkm.tumor.symbol.filter.ccle[
    sub("\\..*", "", rownames(fpkm.tumor.symbol.filter.ccle)) %in% rownames(protein.info.pancancer.breast.match.all), 
    colnames(protein.info.pancancer.breast.match.all)
    ];


ccle.sample.outlier.status.protein.match.pancancer.all <- ccle.sample.outlier.status[
    rownames(fpkm.tumor.symbol.filter.ccle.outlier.pancancer.all), 
    colnames(protein.info.pancancer.breast.match.all)
    ];

rownames(fpkm.tumor.symbol.filter.ccle.outlier.pancancer.all) <- sub(
    "\\..*", "", 
    rownames(fpkm.tumor.symbol.filter.ccle.outlier.pancancer.all)
    );
rownames(ccle.sample.outlier.status.protein.match.pancancer.all) <- sub(
    "\\..*", "", 
    rownames(ccle.sample.outlier.status.protein.match.pancancer.all)
    );

outlier.protein.ccle.zscore.list.pancancer.all <- list();
non.outlier.protein.ccle.zscore.list.pancancer.all <- list();
target.gene.ccle.zscore.list.pancancer.all <- NULL;

for (i in 1:nrow(protein.info.pancancer.breast.match.all)) {
    # Protein target name
    target.gene.name.protein <- rownames(protein.info.pancancer.breast.match.all)[i];
    
    # Outlier patient column
    target.col <- colnames(ccle.sample.outlier.status.protein.match.pancancer.all)[
        ccle.sample.outlier.status.protein.match.pancancer.all[target.gene.name.protein,] == 1
        ];
    non.target.col <- colnames(ccle.sample.outlier.status.protein.match.pancancer.all)[
        ccle.sample.outlier.status.protein.match.pancancer.all[target.gene.name.protein,] == 0
        ];
    
    target.gene.ccle.zscore.list.pancancer.all <- c(
        target.gene.ccle.zscore.list.pancancer.all, 
        target.gene.name.protein
        );
    outlier.protein.ccle.zscore.list.pancancer.all[[i]] <- protein.info.pancancer.breast.match.all[
        i, target.col
        ];
    non.outlier.protein.ccle.zscore.list.pancancer.all[[i]] <- protein.info.pancancer.breast.match.all[
        i, non.target.col
        ];
    }

outlier.protein.ccle.zscore.value.pancancer.all <- as.numeric(
    unlist(outlier.protein.ccle.zscore.list.pancancer.all)
    );
non.outlier.protein.ccle.zscore.value.pancancer.all <- as.numeric(
    unlist(non.outlier.protein.ccle.zscore.list.pancancer.all)
    );

# Box plot - compare the values between patients
# - excluding the genes having no outlier patient info

names(outlier.protein.ccle.zscore.list.pancancer.all) <- target.gene.ccle.zscore.list.pancancer.all;
outlier.protein.ccle.list.no.p.na.pancancer.all <- na.omit(
    unlist(outlier.protein.ccle.zscore.list.pancancer.all)
    );
names(non.outlier.protein.ccle.zscore.list.pancancer.all) <- target.gene.ccle.zscore.list.pancancer.all;
non.outlier.protein.ccle.list.no.p.na.pancancer.all <- non.outlier.protein.ccle.zscore.list.pancancer.all[
    names(outlier.protein.ccle.list.no.p.na.pancancer.all)
    ];

protein.ccle.na.value.pancancer.all <- data.frame(
    protein.ccle.na.value = c(
        as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na.pancancer.all)), 
        as.numeric(unlist(outlier.protein.ccle.list.no.p.na.pancancer.all))
        )
    );
protein.ccle.na.value.box.pancancer.all <- data.frame(
    cbind(
        protein.ccle.na.value.pancancer.all$protein.ccle.na.value, 
        c(
            rep('non', length(as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na.pancancer.all)))), 
            rep('out', length(as.numeric(unlist(outlier.protein.ccle.list.no.p.na.pancancer.all))))
            )
        )
    );
colnames(protein.ccle.na.value.box.pancancer.all) <- c('protein.ccle.value', 'status');
protein.ccle.na.value.box.pancancer.all[,1] <- as.numeric(protein.ccle.na.value.box.pancancer.all[,1]);

wilcox.result.protien.na.pancancer.all <- wilcox.test(
    as.numeric(unlist(outlier.protein.ccle.list.no.p.na.pancancer.all)),
    as.numeric(unlist(non.outlier.protein.ccle.list.no.p.na.pancancer.all)),
    alternative = "two.sided", conf.int = TRUE
    );

text.pvalue.protien.na.pancancer.all <- display.statistical.result(
    x = wilcox.result.protien.na.pancancer.all$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

key.protien.na.pancancer.all <- list(
    text = list(
        lab = text.pvalue.protien.na.pancancer.all, 
        cex = 1
        ),
    x = 0.25,
    y = 0.95
    );


ccle.protein.box.all <- BoutrosLab.plotting.general::create.boxplot(
    formula = protein.ccle.value ~ status,
    data = protein.ccle.na.value.box.pancancer.all,
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
    key = key.protien.na.pancancer.all,
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
    trellis.object = ccle.protein.box.all,
    filename = file.path(output.directory, 'Figure_4_d.png'),
    width = 3.5,
    height = 6.5,
    size.units = 'in',
    resolution = 1200
);
