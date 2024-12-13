### HISTORY ######################################################################
# This script performs driver gene analysis for breast cancer subtypes.
# Date: 2024-08-13

### DESCRIPTION ##################################################################
# This script analyzes driver gene mutations in breast cancer subtypes (Luminal A and B)
# across multiple datasets (METABRIC, TCGA-BRCA, ICGC BRCA-EU). It calculates
# odds ratios, performs meta-analysis, and creates visualizations of the results.

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
library(metafor);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

prepare.data <- function(row, outlier) {
    df <- data.frame(mutation = row, outlier = outlier);
    df$outlier <- as.numeric(df$outlier);
    table <- table(df)

    if (nrow(table) == 1) {
        odd.table <- cbind(sum(table[, 2:ncol(table)]), table[, 1]);
        odd.table <- rbind(
            c(0, 0),
            odd.table
            )
        } else {
        odd.table <- t(rbind(apply(table[, 2:ncol(table)], 1, sum), table[, 1]));
        }

    odd.table <- odd.table + 1
    colnames(odd.table) <- c('out', 'non');
    return(odd.table);
    };


perform.fisher.test <- function(table) {
    test.result <- fisher.test(table, alternative = 'two.sided');
    return(list(
        odds.ratio = test.result$estimate,
        ci = test.result$conf.int,
        p.value = test.result$p.value
        ));
    };


calculate.odds.ratios <- function(mutation.data, outlier.data) {
    results <- apply(mutation.data, 1, function(row) {
        table <- prepare.data(row, outlier.data);
        return(perform.fisher.test(table));
        });

    odds.ratio.unlist <- sapply(results, function(x) x$odds.ratio);
    ci.list <- lapply(results, function(x) x$ci);
    p.value.unlist <- sapply(results, function(x) x$p.value);
    ci.df <- t(sapply(ci.list, c));
    colnames(ci.df) <- c('lower', 'upper');

    return(list(
        odds.ratio = odds.ratio.unlist,
        ci = ci.df,
        p.value = p.value.unlist
        ));
    };


outlier.patient.tag.01.sum.patient.brca <- apply(outlier.patient.tag.01.brca, 2, sum)
outlier.patient.tag.01.sum.patient.meta <- apply(outlier.patient.tag.01.meta, 2, sum);
outlier.patient.tag.01.sum.patient.icgc <- apply(outlier.patient.tag.01.icgc, 2, sum);

outlier.patient.tag.01.t.p.order.sum <- c(
    outlier.patient.tag.01.sum.patient.brca,
    outlier.patient.tag.01.sum.patient.meta,
    outlier.patient.tag.01.sum.patient.icgc
    );

### Get the mutation data
# 1. TCGA-BRCA
brca.mutation.driver.list.gene.vector.data.convert <- brca.mutation.driver.list.gene.vector.data;
brca.mutation.driver.list.gene.vector.data.convert[!(brca.mutation.driver.list.gene.vector.data.convert == '0')] <- 'mutation';
brca.mutation.driver.list.gene.vector.data.convert[brca.mutation.driver.list.gene.vector.data.convert == '0'] <- 'normal';

# 2. METABRIC
meta.mutation.driver.list.gene.vector.data.convert <- meta.mutation.driver.list.gene.vector.data;
meta.mutation.driver.list.gene.vector.data.convert[!(meta.mutation.driver.list.gene.vector.data.convert == '0')] <- 'mutation';
meta.mutation.driver.list.gene.vector.data.convert[meta.mutation.driver.list.gene.vector.data.convert == '0'] <- 'normal';

# 3. ICGC BRCA-EU
icgc.mutation.driver.list.gene.vector.data.convert <- icgc.mutation.driver.list.gene.vector.data.mis;
icgc.mutation.driver.list.gene.vector.data.convert[!(icgc.mutation.driver.list.gene.vector.data.convert == '0')] <- 'mutation';
icgc.mutation.driver.list.gene.vector.data.convert[icgc.mutation.driver.list.gene.vector.data.convert == '0'] <- 'normal';
icgc.mutation.driver.list.gene.vector.data.convert <- icgc.mutation.driver.list.gene.vector.data.convert[
    match(
        rownames(brca.mutation.driver.list.gene.vector.data),
        rownames(icgc.mutation.driver.list.gene.vector.data.convert)
        ),
    ];
rownames(icgc.mutation.driver.list.gene.vector.data.convert) <- rownames(brca.mutation.driver.list.gene.vector.data);
rownames(icgc.mutation.driver.list.gene.vector.data.convert) <- rownames(brca.mutation.driver.list.gene.vector.data);

meta.mutation.driver.list.gene.vector.data.convert.na <- cbind(
    brca.mutation.driver.list.gene.vector.data.convert,
    meta.mutation.driver.list.gene.vector.data.convert,
    icgc.mutation.driver.list.gene.vector.data.convert
    );
meta.mutation.driver.list.gene.vector.data.convert.na <- apply(meta.mutation.driver.list.gene.vector.data.convert.na, 1, function(x) {
    x[is.na(x)] <- 'normal'
    return(x)
    })
meta.mutation.driver.list.gene.vector.data.convert.na <- t(meta.mutation.driver.list.gene.vector.data.convert.na);

# Luma subtype patient
subtype.total.outlier.num.luma <- subtype.total.outlier.num[subtype.total.outlier.num$subtype == 3, ];
outlier.patient.tag.01.t.p.order.sum.luma <- outlier.patient.tag.01.t.p.order.sum[rownames(na.omit(subtype.total.outlier.num.luma))];
meta.mutation.driver.list.gene.vector.data.convert.na.luma <- meta.mutation.driver.list.gene.vector.data.convert.na[
    ,
    colnames(meta.mutation.driver.list.gene.vector.data.convert.na) %in% rownames(na.omit(subtype.total.outlier.num.luma))
    ];

outlier.patient.tag.01.t.p.order.sum.luma.brca <-
    outlier.patient.tag.01.t.p.order.sum.luma[substr(names(outlier.patient.tag.01.t.p.order.sum.luma), 1, 4) == 'TCGA'];
outlier.patient.tag.01.t.p.order.sum.luma.meta <-
    outlier.patient.tag.01.t.p.order.sum.luma[substr(names(outlier.patient.tag.01.t.p.order.sum.luma), 1, 2) == 'MB'];
outlier.patient.tag.01.t.p.order.sum.luma.icgc <-
    outlier.patient.tag.01.t.p.order.sum.luma[substr(names(outlier.patient.tag.01.t.p.order.sum.luma), 1, 2) == 'PR'];

meta.mutation.driver.list.gene.vector.data.convert.na.luma.brca <-
    meta.mutation.driver.list.gene.vector.data.convert.na.luma[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma), 1, 4) == 'TCGA'];
meta.mutation.driver.list.gene.vector.data.convert.na.luma.meta <-
    meta.mutation.driver.list.gene.vector.data.convert.na.luma[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma), 1, 2) == 'MB'];
meta.mutation.driver.list.gene.vector.data.convert.na.luma.icgc <-
    meta.mutation.driver.list.gene.vector.data.convert.na.luma[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma), 1, 2) == 'PR'];

# Lumb subtype patient
subtype.total.outlier.num.lumb <- subtype.total.outlier.num[subtype.total.outlier.num$subtype == 4, ];
outlier.patient.tag.01.t.p.order.sum.lumb <- outlier.patient.tag.01.t.p.order.sum[rownames(na.omit(subtype.total.outlier.num.lumb))];
meta.mutation.driver.list.gene.vector.data.convert.na.lumb <- meta.mutation.driver.list.gene.vector.data.convert.na[
    ,
    colnames(meta.mutation.driver.list.gene.vector.data.convert.na) %in% rownames(na.omit(subtype.total.outlier.num.lumb))
    ];

outlier.patient.tag.01.t.p.order.sum.lumb.brca <-
    outlier.patient.tag.01.t.p.order.sum.lumb[substr(names(outlier.patient.tag.01.t.p.order.sum.lumb), 1, 4) == 'TCGA'];
outlier.patient.tag.01.t.p.order.sum.lumb.meta <-
    outlier.patient.tag.01.t.p.order.sum.lumb[substr(names(outlier.patient.tag.01.t.p.order.sum.lumb), 1, 2) == 'MB'];
outlier.patient.tag.01.t.p.order.sum.lumb.icgc <-
    outlier.patient.tag.01.t.p.order.sum.lumb[substr(names(outlier.patient.tag.01.t.p.order.sum.lumb), 1, 2) == 'PR'];

meta.mutation.driver.list.gene.vector.data.convert.na.lumb.brca <-
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb), 1, 4) == 'TCGA'];
meta.mutation.driver.list.gene.vector.data.convert.na.lumb.meta <-
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb), 1, 2) == 'MB'];
meta.mutation.driver.list.gene.vector.data.convert.na.lumb.icgc <-
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb), 1, 2) == 'PR'];

# 1. TCGA-BRCA
brca.mutation.silent.0 <- c(as.character(brca.mutation.silent), '0');

brca.mutation.driver.list.gene.vector.data.convert.mis <- brca.mutation.driver.list.gene.vector.data;
brca.mutation.driver.list.gene.vector.data.convert.mis[!(brca.mutation.driver.list.gene.vector.data.convert.mis %in% brca.mutation.silent.0)] <- 'mutation';
brca.mutation.driver.list.gene.vector.data.convert.mis[brca.mutation.driver.list.gene.vector.data.convert.mis %in% brca.mutation.silent.0] <- 'normal';

brca.mutation.driver.list.gene.vector.data.convert.mis.luma.brca <-
    brca.mutation.driver.list.gene.vector.data.convert.mis[, colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma.brca)];

# 2. METABRIC
meta.mutation.silent.0 <- c(as.character(meta.mutation.silent), '0');

meta.mutation.driver.list.gene.vector.data.convert.mis <- meta.mutation.driver.list.gene.vector.data;
meta.mutation.driver.list.gene.vector.data.convert.mis[!(meta.mutation.driver.list.gene.vector.data.convert.mis %in% meta.mutation.silent.0)] <- 'mutation';
meta.mutation.driver.list.gene.vector.data.convert.mis[meta.mutation.driver.list.gene.vector.data.convert.mis %in% meta.mutation.silent.0] <- 'normal';

meta.mutation.driver.list.gene.vector.data.convert.mis.luma.meta <-
    meta.mutation.driver.list.gene.vector.data.convert.mis[, colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma.meta)];

# 3. ICGC BRCA-EU
icgc.mutation.silent.0 <- c(as.character(icgc.mutation.silent), '0');

icgc.mutation.driver.list.gene.vector.data.convert.mis <- icgc.mutation.driver.list.gene.vector.data.mis;
icgc.mutation.driver.list.gene.vector.data.convert.mis[!(icgc.mutation.driver.list.gene.vector.data.convert.mis %in% icgc.mutation.silent.0)] <- 'mutation';
icgc.mutation.driver.list.gene.vector.data.convert.mis[icgc.mutation.driver.list.gene.vector.data.convert.mis %in% icgc.mutation.silent.0] <- 'normal';

colnames(icgc.mutation.driver.list.gene.vector.data.convert.mis) <- colnames(icgc.mutation.driver.list.gene.vector.data.mis);

icgc.mutation.driver.list.gene.vector.data.convert.mis.luma.icgc <-
    icgc.mutation.driver.list.gene.vector.data.convert.mis[, colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma.icgc)];

# Ensure the column names are preserved for ICGC BRCA-EU data


### Luminal A

# 1. TCGA-BRCA
brca.results.luma <- calculate.odds.ratios(
    brca.mutation.driver.list.gene.vector.data.convert.mis.luma.brca,
    outlier.patient.tag.01.t.p.order.sum.luma.brca
    );

# 2. ICGC BRCA-EU
icgc.results.luma <- calculate.odds.ratios(
    icgc.mutation.driver.list.gene.vector.data.convert.mis.luma.icgc,
    outlier.patient.tag.01.t.p.order.sum.luma.icgc
    );

# 3. METABRIC
meta.results.luma <- calculate.odds.ratios(
    meta.mutation.driver.list.gene.vector.data.convert.mis.luma.meta,
    outlier.patient.tag.01.t.p.order.sum.luma.meta
    );


brca.luma.odds.ratio.unlist <- brca.results.luma$odds.ratio;
brca.luma.ci.df <- brca.results.luma$ci;
brca.luma.p.unlist <- brca.results.luma$p.value;
names(brca.luma.odds.ratio.unlist) <- mutation.driver.list.gene;

# ICGC
icgc.luma.odds.ratio.unlist <- icgc.results.luma$odds.ratio;
icgc.luma.ci.df <- icgc.results.luma$ci;
icgc.luma.p.unlist <- icgc.results.luma$p.value;
names(icgc.luma.odds.ratio.unlist) <- rownames(icgc.mutation.driver.list.gene.vector.data.mis);

# METABRIC
meta.luma.odds.ratio.unlist <- meta.results.luma$odds.ratio;
meta.luma.ci.df <- meta.results.luma$ci;
meta.luma.p.unlist <- meta.results.luma$p.value;
names(meta.luma.odds.ratio.unlist) <- mutation.driver.list.gene;








convert.to.binary <- function(data) {
    return(ifelse(data == 'mutation', 1, 0));
    }

calculate.sum <- function(data) {
    return(apply(data, 1, sum));
    }

filter.genes <- function(data, threshold) {
    return(data[data > threshold]);
    }

calculate.ln.odds.and.se <- function(odds.ratio, ci) {
    ln.odds <- log(odds.ratio);
    se <- (log(ci[, 2]) - log(ci[, 1])) / 3.92;
    return(list(ln.odds = ln.odds, se = se));
    }



### 1. Luminal A

datasets.luma <- list(
    icgc = icgc.mutation.driver.list.gene.vector.data.convert.mis.luma.icgc,
    meta = meta.mutation.driver.list.gene.vector.data.convert.mis.luma.meta,
    brca = brca.mutation.driver.list.gene.vector.data.convert.mis.luma.brca
    );

results.luma <- list();

for (name in names(datasets.luma)) {
    binary.data <- convert.to.binary(datasets.luma[[name]]);
    sum.data <- calculate.sum(binary.data);
    threshold <- length(get(paste0('outlier.patient.tag.01.t.p.order.sum.luma.', name))) * 0.05;
    results.luma[[name]] <- list(
        sum = sum.data,
        filtered = filter.genes(sum.data, threshold)
        );
    }




names.list <- lapply(results.luma, function(x) names(x$filtered))
all.names <- unlist(names.list)
name.counts <- table(all.names)
common.genes.luma <- names(name.counts[name.counts >= 2])




odds.ratios.luma <- list(
    meta = meta.luma.odds.ratio.unlist[names(results.luma$meta$filtered)][common.genes.luma],
    brca = brca.luma.odds.ratio.unlist[names(results.luma$brca$filtered)][common.genes.luma],
    icgc = icgc.luma.odds.ratio.unlist[names(results.luma$icgc$filtered)][common.genes.luma]
    );
odds.ratios.luma.df <- do.call(rbind, odds.ratios.luma);


p.values.luma <- list(
    meta = meta.luma.p.unlist[names(results.luma$meta$filtered)][common.genes.luma],
    brca = brca.luma.p.unlist[names(results.luma$brca$filtered)][common.genes.luma],
    icgc = icgc.luma.p.unlist[names(results.luma$icgc$filtered)][common.genes.luma]
    );
p.values.luma.df <- do.call(rbind, p.values.luma);

confidence.intervals.luma <- list(
    meta = meta.luma.ci.df[names(results.luma$meta$filtered), ][common.genes.luma, ],
    brca = brca.luma.ci.df[names(results.luma$brca$filtered), ][common.genes.luma, ],
    icgc = icgc.luma.ci.df[names(results.luma$icgc$filtered), ][match(common.genes.luma, rownames(icgc.luma.ci.df[names(results.luma$icgc$filtered), ])), ]
    );


ln.odds.se.luma <- lapply(names(odds.ratios.luma), function(name) {
    calculate.ln.odds.and.se(odds.ratios.luma[[name]], confidence.intervals.luma[[name]]);
    });
names(ln.odds.se.luma) <- names(odds.ratios.luma);


# Meta-analysis
meta.analysis.results.luma <- lapply(common.genes.luma, function(gene) {
    data <- data.frame(
        yi = sapply(ln.odds.se.luma, function(x) x$ln.odds[gene]),
        sei = sapply(ln.odds.se.luma, function(x) x$se[gene])
        );
    data <- data[!is.infinite(data$yi) & !is.infinite(data$sei), ];

    if (nrow(data) > 1) {
        meta.result <- rma.uni(yi = yi, sei = sei, data = data, method = 'DL');
        return(c(exp(meta.result$beta), exp(meta.result$ci.lb), exp(meta.result$ci.ub), meta.result$pval));
        } else {
        return(c(NA, NA, NA, NA));
        };
    });

# Prepare final results
final.results.luma <- do.call(rbind, meta.analysis.results.luma);
colnames(final.results.luma) <- c('odds.ratio', 'ci.lower', 'ci.upper', 'p.value');
rownames(final.results.luma) <- common.genes.luma;

# Filter significant results
final.results.luma <- data.frame(final.results.luma);





# use multiple steps
background.cutoff <- 2;
colourkey.labels.at <- seq(0, background.cutoff, by = 2);
colourkey.labels <- sapply(
    X = colourkey.labels.at,
    FUN = function(x) {
        if (x == 0) {
            return(expression('10'^'0'));
            } else if (x != background.cutoff) {
            return(as.expression(bquote('10'^-.(as.character(x)))));
            } else {
            return(as.expression(bquote('<10'^-.(as.character(x)))));
            }
        }
    );

legend <- legend.grob(
    list(
        legend = list(
            title = expression(underline('P-value')),
            continuous = TRUE,
            colours = c('white', 'black'),
            total.colours = 100,
            labels = colourkey.labels,
            cex = 0.9, # label text size
            at = seq(0, 100, length.out = length(colourkey.labels)),
            height = 3
            )
        ),
    label.cex = 1,
    title.cex = 1,
    title.just = 'left',
    title.fontface = 'plain',
    between.row = 4
    );



# functions to decide the spot size and colour can be included # if not, default functions will be used
spot.size.function <- function(x) {
    0.1 + (2 * abs(x));
    }
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[sign(x) == -1] <- default.colours(2, palette.type = 'dotmap')[1];
    colours[sign(x) == 1] <- default.colours(2, palette.type = 'dotmap')[2];
    return(colours);
    }

dot.each <- create.dotmap(
    x = log2(odds.ratios.luma.df),
    # top.padding = 4,
    main = expression('Association with the driver gene mutation'),
    main.cex = 1.4,
    yaxis.cex = 1.2,
    # yaxis.rot = 90,
    xaxis.rot = 90,
    xaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.lab = c('METABRIC', 'TCGA-BRCA', 'ICGC BRCA-EU'),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    text = list(
                        lab = expression('Odds Ratio'),
                        cex = 1
                        ),
                    padding.text = 4.5
                    )
                ),
            x = 1,
            y = 1
            ),
        inside = list(
            fun = legend,
            x = 1.07,
            y = 0
            )
        ),
    key = list(
        space = 'right',
        points = list(
            cex = spot.size.function(seq(-2, 2, 1)),
            col = spot.colour.function(seq(-2, 2, 1)),
            pch = 19
            ),
        text = list(
            lab = c(
                '0.25', ' 0.5', '1',
                '2', '4'
                ),
            cex = 1,
            adj = 1,
            fontface = 'bold'
            ),
        padding.text = 8
        ),
    key.top = 1,
    right.padding = 2,
    # add borders to points
    pch = 21,
    pch.border.col = 'white',
    # add the background
    bg.data = -log10(p.values.luma.df),
    # add a colourkey
    colourkey = FALSE,
    colourkey.cex = 0.95,
    bg.alpha = 1,
    colour.scheme = c('white', 'black'),
    at = seq(0, 2, 0.01),
    row.colour = 'white',
    col.colour = 'white',
    row.lwd = 1.3,
    col.lwd = 1.3
    );



dot.all <- create.dotmap(
    x = t(log2(final.results.luma$odds.ratio)),
    xaxis.lab = rownames(final.results.luma),
    # main.cex = 1.4,
    yaxis.cex = 1.2,
    # yaxis.rot = 90,
    top.padding = 10,
    xaxis.rot = 90,
    xaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.lab = c('All patients'),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    # control spacing at top of key
    key.top = 1,
    right.padding = 2,
    # add borders to points
    pch = 21,
    pch.border.col = 'white',
    # add the background
    bg.data = t(-log10(final.results.luma$p.value)),
    # add a colourkey
    colourkey = FALSE,
    # colourkey.cex = 0.95,
    bg.alpha = 1,
    # set colour scheme for background data
    colour.scheme = c('white', 'black'),
    at = seq(0, 2, 0.01),
    row.colour = 'white',
    col.colour = 'white',
    row.lwd = 1.3,
    col.lwd = 1.3
    );





dot.multi.luma <- create.multipanelplot(
    list(dot.each, dot.all),
    layout.height = 2,
    layout.width = 1,
    layout.skip = c(FALSE, FALSE),
    plot.objects.heights = c(10.5, 5),
    # xlab.label = xlabel,
    # ylab.label = ylabel,
    x.spacing = -1,
    y.spacing = -10,
    bottom.padding = 0,
    top.padding = 2,
    right.padding = 0
    )
dot.multi.luma

save.outlier.figure(
    dot.multi.luma,
    c('Figure2ef', 'drivergene', 'luma', 'multipanel'),
    width = 6.9,
    height = 5.7
    );

### 2. Luminal B
# 1. TCGA-BRCA
brca.mutation.silent.0 <- c(as.character(brca.mutation.silent), '0');

brca.mutation.driver.list.gene.vector.data.convert.mis <- brca.mutation.driver.list.gene.vector.data;
brca.mutation.driver.list.gene.vector.data.convert.mis[!(brca.mutation.driver.list.gene.vector.data.convert.mis %in% brca.mutation.silent.0)] <- 'mutation';
brca.mutation.driver.list.gene.vector.data.convert.mis[brca.mutation.driver.list.gene.vector.data.convert.mis %in% brca.mutation.silent.0] <- 'normal';

brca.mutation.driver.list.gene.vector.data.convert.mis.lumb.brca <-
    brca.mutation.driver.list.gene.vector.data.convert.mis[, colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb.brca)];

# 2. METABRIC
meta.mutation.silent.0 <- c(as.character(meta.mutation.silent), '0');

meta.mutation.driver.list.gene.vector.data.convert.mis <- meta.mutation.driver.list.gene.vector.data;
meta.mutation.driver.list.gene.vector.data.convert.mis[!(meta.mutation.driver.list.gene.vector.data.convert.mis %in% meta.mutation.silent.0)] <- 'mutation';
meta.mutation.driver.list.gene.vector.data.convert.mis[meta.mutation.driver.list.gene.vector.data.convert.mis %in% meta.mutation.silent.0] <- 'normal';

meta.mutation.driver.list.gene.vector.data.convert.mis.lumb.meta <-
    meta.mutation.driver.list.gene.vector.data.convert.mis[, colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb.meta)];

# 3. ICGC BRCA-EU
icgc.mutation.silent.0 <- c(as.character(icgc.mutation.silent), '0');

icgc.mutation.driver.list.gene.vector.data.convert.mis <- icgc.mutation.driver.list.gene.vector.data.mis;
icgc.mutation.driver.list.gene.vector.data.convert.mis[!(icgc.mutation.driver.list.gene.vector.data.convert.mis %in% icgc.mutation.silent.0)] <- 'mutation';
icgc.mutation.driver.list.gene.vector.data.convert.mis[icgc.mutation.driver.list.gene.vector.data.convert.mis %in% icgc.mutation.silent.0] <- 'normal';

icgc.mutation.driver.list.gene.vector.data.convert.mis.lumb.icgc <-
    icgc.mutation.driver.list.gene.vector.data.convert.mis[, colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb.icgc)];





# 1. TCGA-BRCA
brca.results.lumb <- calculate.odds.ratios(
    brca.mutation.driver.list.gene.vector.data.convert.mis.lumb.brca,
    outlier.patient.tag.01.t.p.order.sum.lumb.brca
    );

# 2. ICGC BRCA-EU
icgc.results.lumb <- calculate.odds.ratios(
    icgc.mutation.driver.list.gene.vector.data.convert.mis.lumb.icgc,
    outlier.patient.tag.01.t.p.order.sum.lumb.icgc
    );

# 3. METABRIC
meta.results.lumb <- calculate.odds.ratios(
    meta.mutation.driver.list.gene.vector.data.convert.mis.lumb.meta,
    outlier.patient.tag.01.t.p.order.sum.lumb.meta
    );


brca.lumb.odds.ratio.unlist <- brca.results.lumb$odds.ratio;
brca.lumb.ci.df <- brca.results.lumb$ci;
brca.lumb.p.unlist <- brca.results.lumb$p.value;
names(brca.lumb.odds.ratio.unlist) <- mutation.driver.list.gene;

# ICGC
icgc.lumb.odds.ratio.unlist <- icgc.results.lumb$odds.ratio;
icgc.lumb.ci.df <- icgc.results.lumb$ci;
icgc.lumb.p.unlist <- icgc.results.lumb$p.value;
names(icgc.lumb.odds.ratio.unlist) <- rownames(icgc.mutation.driver.list.gene.vector.data.mis);

# METABRIC
meta.lumb.odds.ratio.unlist <- meta.results.lumb$odds.ratio;
meta.lumb.ci.df <- meta.results.lumb$ci;
meta.lumb.p.unlist <- meta.results.lumb$p.value;
names(meta.lumb.odds.ratio.unlist) <- mutation.driver.list.gene;








convert.to.binary <- function(data) {
    return(ifelse(data == 'mutation', 1, 0));
    }

calculate.sum <- function(data) {
    return(apply(data, 1, sum));
    }

filter.genes <- function(data, threshold) {
    return(data[data > threshold]);
    }

calculate.ln.odds.and.se <- function(odds.ratio, ci) {
    ln.odds <- log(odds.ratio);
    se <- (log(ci[, 2]) - log(ci[, 1])) / 3.92;
    return(list(ln.odds = ln.odds, se = se));
    }



datasets.lumb <- list(
    icgc = icgc.mutation.driver.list.gene.vector.data.convert.mis.lumb.icgc,
    meta = meta.mutation.driver.list.gene.vector.data.convert.mis.lumb.meta,
    brca = brca.mutation.driver.list.gene.vector.data.convert.mis.lumb.brca
    );

results.lumb <- list();

for (name in names(datasets.lumb)) {
    binary.data <- convert.to.binary(datasets.lumb[[name]]);
    sum.data <- calculate.sum(binary.data);
    threshold <- length(get(paste0('outlier.patient.tag.01.t.p.order.sum.lumb.', name))) * 0.05;
    results.lumb[[name]] <- list(
        sum = sum.data,
        filtered = filter.genes(sum.data, threshold)
        );
    }




names.list <- lapply(results.lumb, function(x) names(x$filtered))
all.names <- unlist(names.list)
name.counts <- table(all.names)
# common.genes.lumb <- names(name.counts[name.counts >= 2])



common.genes.lumb <- union(
    intersect(names(results.lumb$meta$filtered), names(results.lumb$brca$filtered)),
    union(
        intersect(names(results.lumb$meta$filtered), names(results.lumb$icgc$filtered)),
        intersect(names(results.lumb$brca$filtered), names(results.lumb$icgc$filtered))
        )
    );






odds.ratios.lumb <- list(
    meta = meta.lumb.odds.ratio.unlist[names(results.lumb$meta$filtered)][common.genes.lumb],
    brca = brca.lumb.odds.ratio.unlist[names(results.lumb$brca$filtered)][common.genes.lumb],
    icgc = icgc.lumb.odds.ratio.unlist[names(results.lumb$icgc$filtered)][common.genes.lumb]
    );
odds.ratios.lumb.df <- do.call(rbind, odds.ratios.lumb);

p.values.lumb <- list(
    meta = meta.lumb.p.unlist[names(results.lumb$meta$filtered)][common.genes.lumb],
    brca = brca.lumb.p.unlist[names(results.lumb$brca$filtered)][common.genes.lumb],
    icgc = icgc.lumb.p.unlist[names(results.lumb$icgc$filtered)][common.genes.lumb]
    );
p.values.lumb.df <- do.call(rbind, p.values.lumb);

confidence.intervals.lumb <- list(
    meta = meta.lumb.ci.df[names(results.lumb$meta$filtered), ][match(common.genes.lumb, rownames(meta.lumb.ci.df[names(results.lumb$meta$filtered), ])), ],
    brca = brca.lumb.ci.df[names(results.lumb$brca$filtered), ][match(common.genes.lumb, rownames(brca.lumb.ci.df[names(results.lumb$brca$filtered), ])), ],
    icgc = icgc.lumb.ci.df[names(results.lumb$icgc$filtered), ][match(common.genes.lumb, rownames(icgc.lumb.ci.df[names(results.lumb$icgc$filtered), ])), ]
    );


ln.odds.se.lumb <- lapply(names(odds.ratios.lumb), function(name) {
    calculate.ln.odds.and.se(odds.ratios.lumb[[name]], confidence.intervals.lumb[[name]]);
    });
names(ln.odds.se.lumb) <- names(odds.ratios.lumb);


# Meta-analysis
meta.analysis.results.lumb <- lapply(common.genes.lumb, function(gene) {
    data <- data.frame(
        yi = sapply(ln.odds.se.lumb, function(x) x$ln.odds[gene]),
        sei = sapply(ln.odds.se.lumb, function(x) x$se[gene])
        );
    data <- data[!is.infinite(data$yi) & !is.infinite(data$sei), ];

    if (nrow(data) > 1) {
        meta.result <- rma.uni(yi = yi, sei = sei, data = data, method = 'DL');
        return(c(exp(meta.result$beta), exp(meta.result$ci.lb), exp(meta.result$ci.ub), meta.result$pval));
        } else {
        return(c(NA, NA, NA, NA));
        };
    });

# Prepare final results
final.results.lumb <- do.call(rbind, meta.analysis.results.lumb);
colnames(final.results.lumb) <- c('odds.ratio', 'ci.lower', 'ci.upper', 'p.value');
rownames(final.results.lumb) <- common.genes.lumb;

# Filter significant results
final.results.lumb <- data.frame(final.results.lumb);


dot.each <- create.dotmap(
    x = log2(odds.ratios.lumb.df),
    # top.padding = 4,
    main = expression('Association with the driver gene mutation'),
    main.cex = 1.4,
    yaxis.cex = 1.2,
    # yaxis.rot = 90,
    xaxis.rot = 90,
    xaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.lab = c('METABRIC', 'TCGA-BRCA', 'ICGC BRCA-EU'),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    text = list(
                        lab = expression('Odds Ratio'),
                        cex = 1
                        ),
                    padding.text = 4.5
                    )
                ),
            x = 1,
            y = 1
            ),
        inside = list(
            fun = legend,
            x = 1.07,
            y = 0
            )
        ),
    key = list(
        space = 'right',
        points = list(
            cex = spot.size.function(seq(-2, 2, 1)),
            col = spot.colour.function(seq(-2, 2, 1)),
            pch = 19
            ),
        text = list(
            lab = c(
                '0.25', ' 0.5', '1',
                '2', '4'
                ),
            cex = 1,
            adj = 1,
            fontface = 'bold'
            ),
        padding.text = 8
        ),
    key.top = 1,
    right.padding = 2,
    # add borders to points
    pch = 21,
    pch.border.col = 'white',
    # add the background
    bg.data = -log10(p.values.lumb.df),
    # add a colourkey
    colourkey = FALSE,
    colourkey.cex = 0.95,
    bg.alpha = 1,
    colour.scheme = c('white', 'black'),
    at = seq(0, 2, 0.01),
    row.colour = 'white',
    col.colour = 'white',
    row.lwd = 1.3,
    col.lwd = 1.3
    );



dot.all <- create.dotmap(
    x = t(log2(final.results.lumb$odds.ratio)),
    xaxis.lab = rownames(final.results.lumb),
    # main.cex = 1.4,
    yaxis.cex = 1.2,
    # yaxis.rot = 90,
    top.padding = 10,
    xaxis.rot = 90,
    xaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.lab = c('All patients'),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    # control spacing at top of key
    key.top = 1,
    right.padding = 2,
    # add borders to points
    pch = 21,
    pch.border.col = 'white',
    # add the background
    bg.data = t(-log10(final.results.lumb$p.value)),
    # add a colourkey
    colourkey = FALSE,
    # colourkey.cex = 0.95,
    bg.alpha = 1,
    # set colour scheme for background data
    colour.scheme = c('white', 'black'),
    at = seq(0, 2, 0.01),
    row.colour = 'white',
    col.colour = 'white',
    row.lwd = 1.3,
    col.lwd = 1.3
    );

dot.multi.lumb <- create.multipanelplot(
    list(dot.each, dot.all),
    layout.height = 2,
    layout.width = 1,
    layout.skip = c(FALSE, FALSE),
    plot.objects.heights = c(10.5, 5),
    # xlab.label = xlabel,
    # ylab.label = ylabel,
    x.spacing = -1,
    y.spacing = -10,
    bottom.padding = 0,
    top.padding = 2,
    right.padding = 0
    );
dot.multi.lumb

save.outlier.figure(
    dot.multi.lumb,
    c('Figure2ef', 'drivergene', 'lumb', 'multipanel'),
    width = 6.9,
    height = 5.7
    );

# Helper function to subset and organize data by subtype
process.subtype.data <- function(subtype.id) {
    subtype.data <- subtype.total.outlier.num[subtype.total.outlier.num$subtype == subtype.id, ];
    outlier.patients <- outlier.patient.tag.01.t.p.order.sum[rownames(na.omit(subtype.data))]
    mutation.data <- meta.mutation.driver.list.gene.vector.data.convert.na[, rownames(na.omit(subtype.data))]

    outlier.brca <- outlier.patients[substr(names(outlier.patients), 1, 4) == 'TCGA']
    outlier.meta <- outlier.patients[substr(names(outlier.patients), 1, 2) == 'MB']
    outlier.icgc <- outlier.patients[substr(names(outlier.patients), 1, 2) == 'PR']

    mutation.brca <- mutation.data[, substr(colnames(mutation.data), 1, 4) == 'TCGA']
    mutation.meta <- mutation.data[, substr(colnames(mutation.data), 1, 2) == 'MB']
    mutation.icgc <- mutation.data[, substr(colnames(mutation.data), 1, 2) == 'PR']

    list(
        outlier.patients.all = c(outlier.meta, outlier.brca, outlier.icgc),
        mutation.data_all = data.frame(mutation.meta, mutation.brca, mutation.icgc)
        )
    }

# Process data for Luminal A (luma) and Luminal B (lumb)
luma.data <- process.subtype.data(3)
lumb.data <- process.subtype.data(4)

# Define mutated and non-mutated PIK3CA cases
gene.alias <- 'PIK3CA'
luma.mut <- luma.data$outlier.patients.all[luma.data$mutation.data_all[gene.alias, ] == 'mutation']
luma.non <- luma.data$outlier.patients.all[luma.data$mutation.data_all[gene.alias, ] == 'normal']
lumb.mut <- lumb.data$outlier.patients.all[lumb.data$mutation.data_all[gene.alias, ] == 'mutation']
lumb.non <- lumb.data$outlier.patients.all[lumb.data$mutation.data_all[gene.alias, ] == 'normal']

# Perform statistical tests
wilcox.test(luma.mut, luma.non, alternative = 'two.sided', conf.int = TRUE)
wilcox.test(lumb.mut, lumb.non, alternative = 'two.sided', conf.int = TRUE)

# Create density plots
create.density.data <- function(values, bw.value) {
    density.data <- density(na.omit(values), bw = bw.value, from = 0, to = max(as.numeric(values)))
    return(as.data.frame(density.data[c('x', 'y')]))
    }

bw.value <- 0.6
dens.luma.mut.df <- create.density.data(luma.mut, bw.value)
dens.luma.non.df <- create.density.data(luma.non, bw.value)
dens.lumb.mut.df <- create.density.data(lumb.mut, bw.value)
dens.lumb.non.df <- create.density.data(lumb.non, bw.value)

# Combine density data
combine.density <- rbind(
    dens.luma.mut.df, dens.luma.non.df, dens.lumb.mut.df, dens.lumb.non.df
    )

combine.density.group <- cbind(
    combine.density,
    group = c(
        rep('a', nrow(dens.luma.mut.df)),
        rep('b', nrow(dens.luma.non.df)),
        rep('c', nrow(dens.lumb.mut.df)),
        rep('d', nrow(dens.lumb.non.df))
        )
    )

# Create scatter plot
mutation.density <- BoutrosLab.plotting.general::create.scatterplot(
    y ~ log2(x + 1),
    data = combine.density.group,
    type = 'l',
    xlimits = c(-0.05, 6.4),
    groups = combine.density.group$group,
    ylab.label = expression('Density'),
    xlab.label = expression('Number of outlier genes per patient'),
    ylimits = c(-0.002, 0.45),
    yat = seq(0, 3, 0.2),
    xat = seq(0, 10, 2),
    xaxis.lab = expression(2^0, 2^2, 2^4, 2^6),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.cex = 1,
    yaxis.cex = 1,
    cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.4,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    lwd = 3,
    lty = c('dashed', 'solid', 'dashed', 'solid'),
    col = c('red3', 'red3', 'dodgerblue2', 'dodgerblue2'),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = c('red3', 'red3', 'dodgerblue2', 'dodgerblue2'),
                        pch = 21,
                        cex = 2,
                        fill = c('red3', 'red3', 'dodgerblue2', 'dodgerblue2')
                        ),
                    text = list(
                        lab = c(
                            'Luminal A - PIK3CA mutated patients',
                            'Luminal A - PIK3CA non-mutated patients',
                            'Luminal B - PIK3CA mutated patients',
                            'Luminal B - PIK3CA non-mutated patients'
                            )
                        ),
                    padding.text = c(0, 5, 0, 5),
                    cex = 1
                    )
                ),
            x = 0.2,
            y = 0.97,
            draw = FALSE
            )
        ),
    main = expression('Number of outlier genes per patient')
    )

# Save plot and session profile
save.outlier.figure(
    mutation.density,
    c('Figure2ef', 'PIK3CA', 'mutation', 'luma', 'lumb', 'scatter', 'density'),
    width = 6,
    height = 5.5
    )

save.session.profile(file.path('output', 'Figure2ef.txt'));
