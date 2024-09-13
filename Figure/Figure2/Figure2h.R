### HISTORY #####################################################################
# This script processes and merges DNA methylation data. 
# Date: 2024-08-14

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));





# Use promoter region TSS ~ +500bp
me.out.symbol.two.500 <- unique(c(
    rownames(meta.me.outlier.match),
    rownames(brca.me.outlier.match)
    ));


# 2. Merge outlier status data
two.outlier.patient.status.merge.filter.list.500 <- list();

for (i in 1:length(me.out.symbol.two.500)) {
    if (me.out.symbol.two.500[i] %in% rownames(brca.me.outlier.match)) {
        row.brca <- rownames(fpkm.tumor.symbol.filter.brca)[
            fpkm.tumor.symbol.filter.brca$Symbol == me.out.symbol.two.500[i]
            ];
        row.brca <- row.brca[row.brca %in% rownames(outlier.patient.tag.01.brca.me.match)];
        target.gene.brca <- as.numeric(
            outlier.patient.tag.01.brca.me.match[row.brca, ]
            );
        } 
    else {
        target.gene.brca <- rep('NA', ncol(brca.me.outlier.match));
        }
    
    if (me.out.symbol.two.500[i] %in% rownames(meta.me.outlier.match)) {
        target.gene.meta <- as.numeric(
            outlier.patient.tag.01.meta.me.match[
                rownames(fpkm.tumor.symbol.filter.meta.symbol)[
                    fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% me.out.symbol.two.500[i]
                    ], 
                ]
            );
        } 
    else {
        target.gene.meta <- rep('NA', ncol(outlier.patient.tag.01.meta.me.match));
        }
    
    both.target.gene <- c(target.gene.brca, target.gene.meta);
    
    two.outlier.patient.status.merge.filter.list.500[[i]] <- both.target.gene;
    }

outlier.sample.me.two.unlist.mean.500 <- lapply(
    outlier.sample.me.two.500, 
    function(x) { mean(na.omit(as.numeric(x))); }
    );

non.outlier.sample.me.two.unlist.mean.500 <- lapply(
    non.outlier.sample.me.two.500, 
    function(x) { mean(na.omit(as.numeric(x))); }
    );

# Change y-axis: mean(outliers) - mean(non-outliers)
mean.beta.merge.two.500 <- apply(
    two.outlier.promoter.symbol.sample.match.merge.filter.500, 
    1, 
    function(x) { mean(na.omit(as.numeric(x))); }
    );

minus.beta.merge.two.500 <- as.numeric(
    outlier.sample.me.two.unlist.mean.500
    ) - as.numeric(non.outlier.sample.me.two.unlist.mean.500);

# Delta beta > 0.2
dot.colours <- vector(
    length = length(mean.minus.ma.merge.two.500$minus.beta)
    );

dot.colours <- rep(
    'grey60', 
    length(mean.minus.ma.merge.two.500$minus.beta)
    );

dot.colours[mean.minus.ma.merge.two.500$minus.beta < -0.2] <- 'dodgerblue3';
dot.colours[mean.minus.ma.merge.two.500$minus.beta > 0.2] <- 'red2';

p.me <- wilcox.test(
    outlier.sample.me.two.unlist.500, 
    non.outlier.sample.me.two.unlist.500, 
    alternative = "two.sided", 
    conf.int = TRUE
    );

p_value.com <- sprintf("%.1e", p.me$p.value);
p_value_parts.com <- strsplit(p_value.com, split = "e")[[1]];
base.com <- as.numeric(p_value_parts.com[1]);
exponent.com <- as.numeric(p_value_parts.com[2]);

text.pvalue <- display.statistical.result(
    x = p.me$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

key.minus.beta <- list(
    text = list(
        lab = text.pvalue, 
        cex = 1
        ), 
    x = 0.25,
    y = 0.95
    );

me.merge.scatter <- create.scatterplot(
    formula = minus.beta ~ mean.beta,
    data = mean.minus.ma.merge.two.500,
    col = dot.colours,
    alpha = .80,
    ylimits = c(-1.01, 1.01),
    xlimits = c(-0.01, 1.01),
    yat = seq(-1.0, 1.0, 0.5),
    xat = seq(0.0, 1.0, 0.2),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    add.grid = TRUE,
    grid.colour = 'grey80',
    cex = 0.9,
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.2,
    ylab.cex = 1.1,
    main.cex = 1.5,
    left.padding = 0,
    main = expression('DNA methylation on promoter of outlier patients'),
    xlab.label = expression(paste('Mean of ', beta, ' value')),
    ylab.label = expression(paste(Delta, ' Mean of ', beta, ' value')),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = key.minus.beta
                ),
            x = 0.71,
            y = 0.95,
            corner = c(0, 1)
            )
        ),
    abline.h = 0,
    abline.col = 'black',
    abline.lwd = 2,
    abline.lty = 3
    );

save.outlier.figure(
    me.merge.scatter,
    c('methylation', 'diff', 'merge', 'scatter'),
    width = 6,
    height = 4.8
    );
