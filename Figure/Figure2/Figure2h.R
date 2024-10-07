### HISTORY ######################################################################
# This script processes and merges DNA methylation data.
# Date: 2024-08-14

### DESCRIPTION ##################################################################
# This script analyzes DNA methylation data for breast cancer. It merges data
# from multiple sources, processes outlier statuses, and creates visualizations
# of DNA methylation differences between outlier and non-outlier samples. The script
# focuses on promoter regions (TSS ~ +500bp) and includes statistical analysis.

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-10-03_Figure1_2_3_4_min_input.rda'));

load.multiple.computed.variables(c(
    'brca.outlier.symbol',
    'metabric.outlier.symbol'
    ));

# Get DNA methylation data
# 1. TCGA-BRCA
brca.me.outlier.match <- brca.me.data[rownames(brca.me.data) %in% brca.outlier.symbol,];
brca.outlier.non.promoter.symbol.sample.match.merge.500 <- brca.me.data[!(rownames(brca.me.data) %in% brca.outlier.symbol),];
# 2. METABRIC
meta.me.outlier.match <- meta.me.data[rownames(meta.me.data) %in% metabric.outlier.symbol,];
meta.outlier.non.promoter.symbol.sample.match.merge.500 <- meta.me.data[!(rownames(meta.me.data) %in% metabric.outlier.symbol),];


# Use promoter region TSS ~ +500bp
me.out.symbol.two.500 <- unique(c(
    rownames(meta.me.outlier.match),
    rownames(brca.me.outlier.match)
    ));



# merge two dataset
#   - 1. merge methylation data
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- list();
for (i in 1:length(me.out.symbol.two.500)){
    if(me.out.symbol.two.500[i] %in% rownames(brca.me.outlier.match)) {
        target.gene.brca <- as.numeric(brca.me.outlier.match[me.out.symbol.two.500[i],]);
        }
    else {
        target.gene.brca <- rep('NA', ncol(brca.me.outlier.match))
        }
    
    if(me.out.symbol.two.500[i] %in% rownames(meta.me.outlier.match)) {
        target.gene.meta <- as.numeric(meta.me.outlier.match[me.out.symbol.two.500[i],]);
        }
    else {
        target.gene.meta <- rep('NA', ncol(meta.me.outlier.match))
        }
    
    both.target.gene <- c(target.gene.brca, target.gene.meta)
    
    two.outlier.promoter.symbol.sample.match.merge.filter.500[[i]] <- both.target.gene;
    }

two.outlier.promoter.symbol.sample.match.merge.filter.500 <- do.call(rbind, two.outlier.promoter.symbol.sample.match.merge.filter.500);
rownames(two.outlier.promoter.symbol.sample.match.merge.filter.500) <- me.out.symbol.two.500;
colnames(two.outlier.promoter.symbol.sample.match.merge.filter.500) <- c(colnames(brca.me.outlier.match), colnames(meta.me.outlier.match));
two.outlier.promoter.symbol.sample.match.merge.filter.500 <- data.frame(two.outlier.promoter.symbol.sample.match.merge.filter.500);


# 2. Merge outlier status data
outlier.patient.tag.01.brca.me.match <- outlier.patient.tag.01.brca[,colnames(brca.me.outlier.match)];
outlier.patient.tag.01.brca.me.match <- outlier.patient.tag.01.brca.me.match[
    rownames(outlier.patient.tag.01.brca.me.match) %in% rownames(fpkm.tumor.symbol.filter.brca)[
        fpkm.tumor.symbol.filter.brca$Symbol %in% 
            rownames(brca.me.outlier.match)],
    ];
outlier.patient.tag.01.meta.me.match <- outlier.patient.tag.01.meta[,colnames(meta.me.outlier.match)];
outlier.patient.tag.01.meta.me.match <- outlier.patient.tag.01.meta.me.match[
    rownames(outlier.patient.tag.01.meta.me.match) %in% rownames(fpkm.tumor.symbol.filter.meta.symbol)[
        fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% 
            rownames(meta.me.outlier.match)],
    ];

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
        } else {
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
        } else {
        target.gene.meta <- rep('NA', ncol(outlier.patient.tag.01.meta.me.match));
        }

    both.target.gene <- c(target.gene.brca, target.gene.meta);

    two.outlier.patient.status.merge.filter.list.500[[i]] <- both.target.gene;
    }

two.outlier.patient.status.merge.filter.500 <- do.call(rbind, two.outlier.patient.status.merge.filter.list.500);
two.outlier.patient.status.merge.filter.500 <- data.frame(two.outlier.patient.status.merge.filter.500);
rownames(two.outlier.patient.status.merge.filter.500) <- me.out.symbol.two.500;
colnames(two.outlier.patient.status.merge.filter.500) <- c(
    colnames(brca.me.data), 
    colnames(meta.me.data));




# divide into outlier and non-outlier
outlier.sample.me.two.500 <- list();
non.outlier.sample.me.two.500 <- list();
for (i in 1:nrow(two.outlier.promoter.symbol.sample.match.merge.filter.500)) {
    
    outlier.patient.gene <- two.outlier.patient.status.merge.filter.500[i,];
    outlier.me <- two.outlier.promoter.symbol.sample.match.merge.filter.500[i,][outlier.patient.gene == 1];
    non.outlier.me <- two.outlier.promoter.symbol.sample.match.merge.filter.500[i,][outlier.patient.gene == 0];
    
    outlier.sample.me.two.500[[i]] <- outlier.me;
    non.outlier.sample.me.two.500[[i]] <- non.outlier.me;
    
    }

outlier.sample.me.two.unlist.500 <- as.numeric(unlist(outlier.sample.me.two.500));
non.outlier.sample.me.two.unlist.500 <- as.numeric(unlist(non.outlier.sample.me.two.500));

outlier.sample.me.two.unlist.mean.500 <- lapply(outlier.sample.me.two.500, function(x) {mean(na.omit(as.numeric(x)))});
non.outlier.sample.me.two.unlist.mean.500 <- lapply(non.outlier.sample.me.two.500, function(x) {mean(na.omit(as.numeric(x)))});

mean.beta.merge.two.500 <- apply(
    two.outlier.promoter.symbol.sample.match.merge.filter.500, 
    1, 
    function(x) {mean(na.omit(as.numeric(x)))}
    );
minus.beta.merge.two.500 <- as.numeric(outlier.sample.me.two.unlist.mean.500) - as.numeric(non.outlier.sample.me.two.unlist.mean.500);
mean.minus.ma.merge.two.500 <- data.frame(cbind(
    as.numeric(mean.beta.merge.two.500), 
    as.numeric(minus.beta.merge.two.500), 
    rownames(two.outlier.promoter.symbol.sample.match.merge.filter.500)
    ));
mean.minus.ma.merge.two.500[,1] <- as.numeric(mean.minus.ma.merge.two.500[,1]);
mean.minus.ma.merge.two.500[,2] <- as.numeric(mean.minus.ma.merge.two.500[,2]);
colnames(mean.minus.ma.merge.two.500) <- c('mean.beta', 'minus.beta', 'Symbol');

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
    alternative = 'two.sided',
    conf.int = TRUE
    );

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
    c('Figure2h', 'methylation', 'diff', 'merge', 'scatter'),
    width = 6,
    height = 4.8
    );

cache.multiple.computed.variables(c(
    'p.me',
    'brca.outlier.non.promoter.symbol.sample.match.merge.500',
    'me.out.symbol.two.500',
    'meta.me.outlier.match',
    'meta.outlier.non.promoter.symbol.sample.match.merge.500',
    'non.outlier.sample.me.two.500',
    'outlier.patient.tag.01.brca.me.match',
    'outlier.patient.tag.01.meta.me.match',
    'outlier.sample.me.two.500',
    'two.outlier.patient.status.merge.filter.500',
    'two.outlier.promoter.symbol.sample.match.merge.filter.500'
    ));

save.session.profile(file.path('output', 'Figure2h.txt'));
