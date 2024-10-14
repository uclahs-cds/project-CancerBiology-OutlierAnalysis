library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'outlier.symbol'
    ));

gene.dependency.breast.t <- t(gene.dependency.breast);
gene.dependency.breast.t.num.match <- as.data.frame(apply(gene.dependency.breast.t, 2, as.numeric));
rownames(gene.dependency.breast.t.num.match) <- rownames(gene.dependency.breast.t);
colnames(gene.dependency.breast.t.num.match) <- colnames(gene.dependency.breast.t);
gene.dependency.breast.t.num.match <- gene.dependency.breast.t.num.match[
    , colnames(fpkm.tumor.symbol.filter.ccle)
    ];

gene.dependency.breast.t.num.match.05 <- gene.dependency.breast.t.num.match[rownames(ccle.outlier.rank.fdr.05), ];
gene.dependency.breast.t.num.match.05.na <- na.omit(gene.dependency.breast.t.num.match.05);
ccle.sample.outlier.status.overlap <- ccle.sample.outlier.status[rownames(ccle.outlier.rank.fdr.05), ];
ccle.sample.outlier.status.overlap.na <- ccle.sample.outlier.status.overlap[rownames(gene.dependency.breast.t.num.match.05.na), ];

# Filter for FDR < 0.05 and match gene names
gene.rnai.breast.t.num.match.05 <- rnai.effect.breast[
    rownames(rnai.effect.breast) %in% gsub('\\..*$', '', rownames(ccle.outlier.rank.fdr.05)),
    ];
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05;
sample.outlier.05.overlap <- ccle.sample.outlier.status.overlap
sample.outlier.05.overlap.na <- sample.outlier.05.overlap[match(rownames(gene.rnai.breast.t.num.match.05.na), gsub('\\..*$', '', rownames(sample.outlier.05.overlap))), ];
sample.outlier.05.overlap.na <- sample.outlier.05.overlap.na[, colnames(gene.rnai.breast.t.num.match.05.na)]; # should chnage the name
sample.outlier.05.overlap.na.rnai <- sample.outlier.05.overlap.na[, colnames(gene.rnai.breast.t.num.match.05.na)];
rownames(gene.rnai.breast.t.num.match.05.na) <- rownames(sample.outlier.05.overlap.na);
sample.outlier.05.overlap.na.sum <- apply(sample.outlier.05.overlap.na, 1, sum);

sample.outlier.05.overlap.na <- sample.outlier.05.overlap.na[sample.outlier.05.overlap.na.sum > 0, ];
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05.na[sample.outlier.05.overlap.na.sum > 0, ];


rnai.quantile.05 <- list();
outlier.gene.rnai.score.05 <- list();
nonoutlier.gene.rnai.score.05 <- list();
for (i in 1:nrow(sample.outlier.05.overlap.na)) {
    outlier.gene <- gene.rnai.breast.t.num.match.05.na[i, which(sample.outlier.05.overlap.na[i, ] == 1), drop = FALSE];
    non.outlier.gene <- gene.rnai.breast.t.num.match.05.na[i, -(which(sample.outlier.05.overlap.na[i, ] == 1))];

    ecdf.obj <- ecdf(as.numeric(non.outlier.gene));
    quantile.value <- ecdf.obj(outlier.gene);
    quantile.value <- t(data.frame(quantile.value));
    colnames(quantile.value) <- colnames(outlier.gene);
    rownames(quantile.value) <- rownames(outlier.gene);
    rnai.quantile.05[[i]] <- quantile.value;

    outlier.gene.rnai.score.05[[i]] <- outlier.gene;
    nonoutlier.gene.rnai.score.05[[i]] <- non.outlier.gene;
    }



# Calculate mean RNAi scores and differences
outlier.gene.rnai.score.05.mean <- sapply(outlier.gene.rnai.score.05, function(x) {
    mean(na.omit(unlist(x)));
    });
outlier.gene.rnai.score.05.mean <- data.frame(rnai = outlier.gene.rnai.score.05.mean);
rownames(outlier.gene.rnai.score.05.mean) <- rownames(sample.outlier.05.overlap.na);

nonoutlier.gene.rnai.score.05.mean <- sapply(nonoutlier.gene.rnai.score.05, function(x) {
    mean(na.omit(unlist(x)));
    });
nonoutlier.gene.rnai.score.05.mean <- data.frame(rnai = nonoutlier.gene.rnai.score.05.mean);
rownames(nonoutlier.gene.rnai.score.05.mean) <- rownames(sample.outlier.05.overlap.na);

# Create difference matrix for RNAi scores
gene.rnai.diff.matrix.05 <- data.frame(
    out = outlier.gene.rnai.score.05.mean$rnai,
    non = nonoutlier.gene.rnai.score.05.mean$rnai,
    diff = outlier.gene.rnai.score.05.mean$rnai - nonoutlier.gene.rnai.score.05.mean$rnai,
    symbol = sub('\\..*', '', rownames(outlier.gene.rnai.score.05.mean))
    );
rownames(gene.rnai.diff.matrix.05) <- rownames(outlier.gene.rnai.score.05.mean);

gene.rnai.diff.matrix.05.overlap <- gene.rnai.diff.matrix.05[gene.rnai.diff.matrix.05$symbol %in% outlier.symbol$unique, ];


ccle.sample.outlier.status.fdr.05 <- ccle.sample.outlier.status[rownames(ccle.outlier.rank.fdr.05), ];
ccle.sample.outlier.status.fdr.05.five <- ccle.sample.outlier.status.fdr.05[sub('\\..*', '', rownames(ccle.sample.outlier.status.fdr.05)) %in% outlier.symbol$unique, ];
ccle.sample.outlier.status.fdr.05.five.symbol <- sub('\\..*', '', rownames(ccle.sample.outlier.status.fdr.05.five));





# Score of the overlap outliers
gene.rnai.diff.matrix.05.overlap.minus.05 <- gene.rnai.diff.matrix.05.overlap[gene.rnai.diff.matrix.05.overlap$diff < -0.4, ];
gene.rnai.diff.matrix.05.overlap.minus.05 <- na.omit(gene.rnai.diff.matrix.05.overlap.minus.05);
gene.rnai.diff.matrix.05.overlap.minus.05.order <- gene.rnai.diff.matrix.05.overlap.minus.05[order(gene.rnai.diff.matrix.05.overlap.minus.05$diff), ];


rnai.score.05.overlap.minus.05 <- gene.rnai.breast.t.num.match.05.na[rownames(gene.rnai.diff.matrix.05.overlap.minus.05.order), ];
outlier.status.05.overlap.minus.05 <- ccle.sample.outlier.status.overlap.na[rownames(gene.rnai.diff.matrix.05.overlap.minus.05.order), colnames(rnai.score.05.overlap.minus.05)];

gene.name.box.05 <- NULL;
for (i in 1:nrow(outlier.status.05.overlap.minus.05)) {
    gene.name <- rep(sub('\\..*', '', rownames(rnai.score.05.overlap.minus.05))[i], ncol(rnai.score.05.overlap.minus.05));
    gene.name.box.05 <- c(gene.name.box.05, gene.name);
    }

rnai.05.box <- data.frame(cbind(
    score = as.numeric(unlist(t(rnai.score.05.overlap.minus.05))),
    gene = gene.name.box.05,
    status = as.numeric(unlist(t(outlier.status.05.overlap.minus.05)))
    ));
rnai.05.box$score <- as.numeric(rnai.05.box$score);
rnai.05.box$status <- as.numeric(rnai.05.box$status);



cache.multiple.computed.variables(c(
    'ccle.sample.outlier.status.fdr.05.five',
    'ccle.sample.outlier.status.fdr.05.five.symbol',
    'ccle.sample.outlier.status.overlap',
    'ccle.sample.outlier.status.overlap.na',
    'gene.dependency.breast.t.num.match.05.na',
    'gene.rnai.breast.t.num.match.05.na',
    'gene.rnai.diff.matrix.05.overlap',
    'rnai.05.box',
    'rnai.score.05.overlap.minus.05',
    'sample.outlier.05.overlap.na'
    ));


save.session.profile(file.path('output', 'cell.line.analysis.txt'));
