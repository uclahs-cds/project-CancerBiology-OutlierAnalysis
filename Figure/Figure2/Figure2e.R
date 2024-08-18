### DRIVER GENE ANALYSIS ######################################################
# Date: 2024-08-13

### PREAMBLE ###################################################################
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
library(metafor);


load(file = '2023-11-18_driver_gene_mutation_meta.rda');
load(file = '2023-11-18_driver_gene_mutation_brca.rda');
load(file = '2024-01-01_driver_gene_mutation_icgc.rda');





prepare.data <- function(row, outlier) {
    df <- data.frame(mutation = row, outlier = outlier);
    df$mutation <- as.numeric(df$mutation);
    df$outlier <- as.numeric(df$outlier);
    table <- table(df) + 1;  # Add pseudocount
    
    if (nrow(table) == 1) {
        table <- rbind(c(1, 1), table);
        };
    
    colnames(table) <- c('out', 'non');
    return(table);
    };


perform.fisher.test <- function(table) {
    test.result <- fisher.test(table, alternative = "two.sided");
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
    colnames(ci.df) <- c("lower", "upper");
    
    return(list(
        odds.ratio = odds.ratio.unlist,
        ci = ci.df,
        p.value = p.value.unlist
        ));
    };


### Luminal A

# 1. TCGA-BRCA
brca.results.luma <- calculate.odds.ratios(
    brca.mutation.driver.list.gene.vector.data.convert.mis.luma.brca,
    outlier.patient.tag.01.t.p.order.sum.luma.brca
    );

# 2. ICGC BRCA-EU
icgc.results.luma <- calculate.odds.ratios(
    icgc.mutation.driver.list.gene.vector.data.convert.mis[
        , colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma.icgc)
        ],
    outlier.patient.tag.01.t.p.order.sum.luma.icgc
    );

# 3. METABRIC
meta.results.luma <- calculate.odds.ratios(
    meta.mutation.driver.list.gene.vector.data.convert.mis[
        , colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma.meta)
        ],
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
names(icgc.luma.odds.ratio.unlist) <- rownames(icgc.mutation.driver.list.gene.vector.data);

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
    se <- (log(ci[,2]) - log(ci[,1])) / 3.92;
    return(list(ln.odds = ln.odds, se = se));
    }



### 1. Luminal A

datasets.luma <- list(
    icgc = icgc.mutation.driver.list.gene.vector.data.convert.mis.luma.icgc,
    meta = meta.mutation.driver.list.gene.vector.data.convert.mis.luma.meta,
    brca = brca.mutation.driver.list.gene.vector.data.convert.mis.luma.brca
    );

results.luma <- list();

for (name in names(datasets)) {
    binary.data <- convert.to.binary(datasets[[name]]);
    sum.data <- calculate.sum(binary.data);
    threshold <- length(get(paste0("outlier.patient.tag.01.t.p.order.sum.luma.", name))) * 0.05;
    results[[name]] <- list(
        sum = sum.data,
        filtered = filter.genes(sum.data, threshold)
        );
    }


common.genes.luma <- Reduce(intersect, lapply(results.luma, function(x) names(x$filtered)));


odds.ratios.luma <- list(
    meta = meta.luma.odds.ratio.unlist[common.genes.luma],
    brca = brca.luma.odds.ratio.unlist[common.genes.luma],
    icgc = icgc.luma.odds.ratio.unlist[common.genes.luma]
    );

p.values.luma <- list(
    meta = meta.luma.p.unlist[common.genes.luma],
    brca = brca.luma.p.unlist[common.genes.luma],
    icgc = icgc.luma.p.unlist[common.genes.luma]
    );

confidence.intervals.luma <- list(
    meta = meta.luma.ci.df[common.genes.luma,],
    brca = brca.luma.ci.df[common.genes.luma,],
    icgc = icgc.luma.ci.df[common.genes.luma,]
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
    data <- data[!is.infinite(data$yi) & !is.infinite(data$sei),];
    
    if (nrow(data) > 1) {
        meta.result <- rma.uni(yi = yi, sei = sei, data = data, method = 'DL');
        return(c(exp(meta.result$beta), exp(meta.result$ci.lb), exp(meta.result$ci.ub), meta.result$pval));
    } else {
        return(c(NA, NA, NA, NA));
    };
});

# Prepare final results
final.results.luma <- do.call(rbind, meta.analysis.results.luma);
colnames(final.results.luma) <- c("odds.ratio", "ci.lower", "ci.upper", "p.value");
rownames(final.results.luma) <- common.genes.luma;

# Filter significant results
significant.genes.luma <- rownames(final.results.luma)[final.results.luma[,"p.value"] < 0.05];






# use multiple steps
background.cutoff <- 2;
colourkey.labels.at <- seq(0, background.cutoff, by = 2);
colourkey.labels <- sapply(
    X = colourkey.labels.at,
    FUN = function(x) {
        if (x == 0) {
            return(expression('10'^'0'));
            } 
        else if ( x != background.cutoff) {
            return(as.expression(bquote('10'^-.(as.character(x)))));
            } 
        else {
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
spot.size.function <- function(x) { 0.1 + (2 * abs(x)); }
spot.colour.function <- function(x) {
    colours <- rep("white", length(x));
    colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1]; 
    colours[sign(x) == 1] <- default.colours(2, palette.type = "dotmap")[2]; 
    return(colours);
    }

dot.each <- create.dotmap(
    x = t(log2(odds.ratios.luma)),
    # top.padding = 4,
    main = expression('Association with the driver.pseudo gene mutation'),
    main.cex = 1.4,
    yaxis.cex = 1.2,
    # yaxis.rot = 90,
    xaxis.rot = 90,
    xaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.lab = c('METABRIC', 'TCGA-BRCA', 'ICGC BRCA-EU'),
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    
    legend = list(
        
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    text = list(
                        lab =expression('Odds Ratio'),
                        cex = 1
                    ),
                    padding.text = 4.5
                )
            ),
            x = 1,
            y = 1

            ),
        
        inside = list(fun = legend,
                      x = 1.07,
                      y = 0
                )
        
            ),
    key = list(
        space = "right",
        points = list(
            cex = spot.size.function(seq(-2, 2, 1)),
            col = spot.colour.function(seq(-2, 2, 1)),
            pch = 19
            ),
        text = list(
                lab = c( "0.25", " 0.5", "1",
                     "2", "4"),
            cex = 1,
            adj = 1,
            fontface = "bold"
            ),
        padding.text = 8
        ),
    key.top = 1,
    right.padding = 2,
    # add borders to points
    pch = 21,
    pch.border.col = "white",
    # add the background
    bg.data = t(-log10(p.values.luma)),
    # add a colourkey
    colourkey = FALSE,
    colourkey.cex = 0.95,
    bg.alpha = 1,
    colour.scheme = c("white", "black"),
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
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    # control spacing at top of key
    key.top = 1,
    right.padding = 2,
    # add borders to points
    pch = 21,
    pch.border.col = "white",
    # add the background
    bg.data = t(-log10(final.results.luma$p.value)),
    # add a colourkey
    colourkey = FALSE,
    # colourkey.cex = 0.95,
    bg.alpha = 1,
    # set colour scheme for background data
    colour.scheme = c("white", "black"),
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

# Save the multi plot as a PDF
pdf(
    file = generate.filename(
        'drivergene_luma', 
        'multipanel', 
        'pdf'
        ), 
    width = 6.9, 
    height = 5.7
    );
dot.multi.luma;
dev.off();

# Save the multi plot as a PNG
png(
    file = generate.filename(
        'drivergene_luma', 
        'multipanel', 
        'png'
        ), 
    width = 6.9, 
    height = 5.7,
    unit = 'in', 
    res = 1200
    );
dot.multi.luma;
dev.off();





### 2. Luminal B

datasets.lumb <- list(
    icgc = icgc.mutation.driver.list.gene.vector.data.convert.mis.lumb.icgc,
    meta = meta.mutation.driver.list.gene.vector.data.convert.mis.lumb.meta,
    brca = brca.mutation.driver.list.gene.vector.data.convert.mis.lumb.brca
    );

results.lumb <- list();

for (name in names(datasets)) {
    binary.data <- convert.to.binary(datasets[[name]]);
    sum.data <- calculate.sum(binary.data);
    threshold <- length(get(paste0("outlier.patient.tag.01.t.p.order.sum.lumb.", name))) * 0.05;
    results[[name]] <- list(
        sum = sum.data,
        filtered = filter.genes(sum.data, threshold)
        );
    }


common.genes.lumb <- Reduce(intersect, lapply(results.lumb, function(x) names(x$filtered)));


odds.ratios.lumb <- list(
    meta = meta.lumb.odds.ratio.unlist[common.genes.lumb],
    brca = brca.lumb.odds.ratio.unlist[common.genes.lumb],
    icgc = icgc.lumb.odds.ratio.unlist[common.genes.lumb]
    );

p.values.lumb <- list(
    meta = meta.lumb.p.unlist[common.genes.lumb],
    brca = brca.lumb.p.unlist[common.genes.lumb],
    icgc = icgc.lumb.p.unlist[common.genes.lumb]
    );

confidence.intervals.lumb <- list(
    meta = meta.lumb.ci.df[common.genes.lumb,],
    brca = brca.lumb.ci.df[common.genes.lumb,],
    icgc = icgc.lumb.ci.df[common.genes.lumb,]
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
    data <- data[!is.infinite(data$yi) & !is.infinite(data$sei),];
    
    if (nrow(data) > 1) {
        meta.result <- rma.uni(yi = yi, sei = sei, data = data, method = 'DL');
        return(c(exp(meta.result$beta), exp(meta.result$ci.lb), exp(meta.result$ci.ub), meta.result$pval));
        } 
    else {
        return(c(NA, NA, NA, NA));
        };
    });

# Prepare final results
final.results.lumb <- do.call(rbind, meta.analysis.results.lumb);
colnames(final.results.lumb) <- c("odds.ratio", "ci.lower", "ci.upper", "p.value");
rownames(final.results.lumb) <- common.genes.lumb;

# Filter significant results
significant.genes.lumb <- rownames(final.results.lumb)[final.results.lumb[,"p.value"] < 0.05];




dot.each <- create.dotmap(
    x = t(log2(odds.ratios.lumb)),
    # top.padding = 4,
    main = expression('Association with the driver.pseudo gene mutation'),
    main.cex = 1.4,
    yaxis.cex = 1.2,
    # yaxis.rot = 90,
    xaxis.rot = 90,
    xaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.lab = c('METABRIC', 'TCGA-BRCA', 'ICGC BRCA-EU'),
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    
    legend = list(
        
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    text = list(
                        lab =expression('Odds Ratio'),
                        cex = 1
                    ),
                    padding.text = 4.5
                    )
                ),
            x = 1,
            y = 1

            ),
        
        inside = list(fun = legend,
                      x = 1.07,
                      y = 0
                    )
        
                ),
    key = list(
        space = "right",
        points = list(
            cex = spot.size.function(seq(-2, 2, 1)),
            col = spot.colour.function(seq(-2, 2, 1)),
            pch = 19
            ),
        text = list(
                lab = c( "0.25", " 0.5", "1",
                     "2", "4"),
            cex = 1,
            adj = 1,
            fontface = "bold"
            ),
        padding.text = 8
        ),
    key.top = 1,
    right.padding = 2,
    # add borders to points
    pch = 21,
    pch.border.col = "white",
    # add the background
    bg.data = t(-log10(p.values.lumb)),
    # add a colourkey
    colourkey = FALSE,
    colourkey.cex = 0.95,
    bg.alpha = 1,
    colour.scheme = c("white", "black"),
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
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    spot.size.function = spot.size.function,
    spot.colour.function = spot.colour.function,
    # control spacing at top of key
    key.top = 1,
    right.padding = 2,
    # add borders to points
    pch = 21,
    pch.border.col = "white",
    # add the background
    bg.data = t(-log10(final.results.lumb$p.value)),
    # add a colourkey
    colourkey = FALSE,
    # colourkey.cex = 0.95,
    bg.alpha = 1,
    # set colour scheme for background data
    colour.scheme = c("white", "black"),
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

# Save the multi plot as a PDF
pdf(
    file = generate.filename(
        'drivergene_lumb', 
        'multipanel', 
        'pdf'
        ), 
    width = 6.9, 
    height = 5.7
    );
dot.multi.lumb;
dev.off();

# Save the multi plot as a PNG
png(
    file = generate.filename(
        'drivergene_lumb', 
        'multipanel', 
        'png'
        ), 
    width = 6.9, 
    height = 5.7,
    unit = 'in', 
    res = 1200
    );
dot.multi.lumb;
dev.off();



