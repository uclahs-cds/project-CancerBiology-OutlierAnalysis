library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'outlier.symbol',
    'ccle.sample.outlier.status.overlap.na',
    'gene.rnai.diff.matrix.05.overlap'
    ));

# Get gene effect score from Cas-CRISPR dataset
cas.effect.breast.05 <- cas.effect.breast[rownames(ccle.outlier.rank.fdr.05), ];
cas.effect.breast.05.na <- na.omit(cas.effect.breast.05);

effect.quantile.05 <- list();
outlier.gene.effect.score.05 <- list();
nonoutlier.gene.effect.score.05 <- list();
for (i in 1:nrow(ccle.sample.outlier.status.overlap.na)) {
    outlier.gene <- cas.effect.breast.05.na[i, which(ccle.sample.outlier.status.overlap.na[i, ] == 1), drop = FALSE];
    non.outlier.gene <- cas.effect.breast.05.na[i, -(which(ccle.sample.outlier.status.overlap.na[i, ] == 1))];
    ecdf.obj <- ecdf(as.numeric(non.outlier.gene));
    quantile.value <- ecdf.obj(outlier.gene);
    ecdf.obj <- ecdf(as.numeric(non.outlier.gene));
    quantile.value <- ecdf.obj(outlier.gene);
    quantile.value <- t(data.frame(quantile.value));
    colnames(quantile.value) <- colnames(outlier.gene);
    rownames(quantile.value) <- rownames(outlier.gene);
    effect.quantile.05[[i]] <- quantile.value;

    outlier.gene.effect.score.05[[i]] <- outlier.gene;
    nonoutlier.gene.effect.score.05[[i]] <- non.outlier.gene;
    }

#   - Use absolute difference of gene effect
outlier.gene.effect.score.05.mean <- sapply(outlier.gene.effect.score.05, function(x) {
    mean(unlist(x)) # Convert inner list to a numeric vector and calculate mean
    });
outlier.gene.effect.score.05.mean <- data.frame(effect = outlier.gene.effect.score.05.mean);
rownames(outlier.gene.effect.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na);

nonoutlier.gene.effect.score.05.mean <- sapply(nonoutlier.gene.effect.score.05, function(x) {
    mean(unlist(x)) # Convert inner list to a numeric vector and calculate mean
    });
nonoutlier.gene.effect.score.05.mean <- data.frame(effect = nonoutlier.gene.effect.score.05.mean);
rownames(nonoutlier.gene.effect.score.05.mean) <- rownames(ccle.sample.outlier.status.overlap.na);

gene.effect.diff.matrix.05 <- data.frame(
    out = outlier.gene.effect.score.05.mean$effect,
    non = nonoutlier.gene.effect.score.05.mean$effect,
    diff = outlier.gene.effect.score.05.mean$effect - nonoutlier.gene.effect.score.05.mean$effect,
    symbol = sub('\\..*', '', rownames(outlier.gene.effect.score.05.mean))
    );
rownames(gene.effect.diff.matrix.05) <- rownames(outlier.gene.effect.score.05.mean);

gene.effect.diff.matrix.05.overlap <- gene.effect.diff.matrix.05[gene.effect.diff.matrix.05$symbol %in% outlier.symbol$unique, ];

cache.multiple.computed.variables(c(
    'cas.effect.breast.05.na',
    'gene.effect.diff.matrix.05.overlap'
    ));

save.session.profile(file.path('output', '5.crispr.rnai.analysis.txt'));
