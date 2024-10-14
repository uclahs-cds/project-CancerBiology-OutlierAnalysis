library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'outlier.symbol'
    ));

### 1. METABRIC

# Filter meta.outlier.symbol and subset the meta.cnv.chr.new.gis dataset
meta.outlier.symbol <- fpkm.tumor.symbol.filter.meta.symbol[rownames(outlier.patient.tag.01.meta), 'Symbol']
meta.cnv.chr.new.gis.match <- meta.cnv.chr.new.gis[na.omit(match(meta.outlier.symbol, meta.cnv.chr.new.gis$Hugo_Symbol)), ]
rownames(meta.cnv.chr.new.gis.match) <- meta.cnv.chr.new.gis.match$Hugo_Symbol

# Subset columns based on common columns with outlier.patient.tag.01.meta
cnv_columns_match <- colnames(meta.cnv.chr.new.gis.match) %in% colnames(outlier.patient.tag.01.meta)
meta.cnv.chr.new.gis.match <- meta.cnv.chr.new.gis.match[, cnv_columns_match]

# Match outlier symbols with row names of outlier.patient.tag.01.meta and remove NAs
meta.outlier.symbol.match <- na.omit(rownames(outlier.patient.tag.01.meta)[match(rownames(meta.cnv.chr.new.gis.match), meta.outlier.symbol)])

# Ensure we remove any `NA`s in column matches to avoid undefined columns
valid_colnames <- colnames(outlier.patient.tag.01.meta) %in% colnames(meta.cnv.chr.new.gis.match)
outlier.patient.tag.01.meta.cnv.match <- outlier.patient.tag.01.meta[meta.outlier.symbol.match, valid_colnames, drop = FALSE]

# Use lapply instead of a for loop for better performance
outlier.nonoutlier_samples <- lapply(seq_len(nrow(meta.cnv.chr.new.gis.match)), function(i) {
    out.patient <- colnames(outlier.patient.tag.01.meta.cnv.match)[outlier.patient.tag.01.meta.cnv.match[i, ] == 1]
    non.out.patient <- colnames(meta.cnv.chr.new.gis.match)[!(colnames(meta.cnv.chr.new.gis.match) %in% out.patient)]

    list(
        outlier = meta.cnv.chr.new.gis.match[i, out.patient, drop = FALSE],
        non_outlier = meta.cnv.chr.new.gis.match[i, non.out.patient, drop = FALSE]
        )
    })

# Split results into separate lists for outlier and non-outlier samples
meta.outlier.sample.cnv.new.gis <- lapply(outlier.nonoutlier_samples, `[[`, 'outlier')
meta.non.outlier.sample.cnv.new.gis <- lapply(outlier.nonoutlier_samples, `[[`, 'non_outlier')

# Set the names of the lists
names(meta.outlier.sample.cnv.new.gis) <- rownames(meta.cnv.chr.new.gis.match)
names(meta.non.outlier.sample.cnv.new.gis) <- rownames(meta.cnv.chr.new.gis.match)

# Convert lists to numeric vectors and calculate medians for meta data
meta.outlier.cnv <- as.numeric(
    unlist(meta.outlier.sample.cnv.new.gis)
    ); # Outlier sample CNA values
meta.non.outlier.cnv <- as.numeric(
    unlist(meta.non.outlier.sample.cnv.new.gis)
    ); # Non-outlier sample CNA values

# Combine outlier and non-outlier values into a single data frame
meta.gis.value <- data.frame(new.gis.value = c(meta.non.outlier.cnv, meta.outlier.cnv));
meta.gis.box <- data.frame(
    cbind(
        meta.gis.value$new.gis.value,
        c(
            rep('non', length(meta.non.outlier.cnv)),
            rep('out', length(meta.outlier.cnv))
            )
        )
    );
colnames(meta.gis.box) <- c('new.gis.value', 'status');
meta.gis.box[, 1] <- as.numeric(meta.gis.box[, 1]);

# Create frequency tables for the CNV values
meta.gis.table.df <- data.frame(table(meta.gis.box));
meta.gis.table <- data.frame(table(meta.gis.box));

# Calculate the relative frequencies
meta.gis.table$Freq[meta.gis.table$status == 'non'] <- meta.gis.table$Freq[meta.gis.table$status == 'non'] / sum(meta.gis.table$Freq[meta.gis.table$status == 'non']);
meta.gis.table$Freq[meta.gis.table$status == 'out'] <- meta.gis.table$Freq[meta.gis.table$status == 'out'] / sum(meta.gis.table$Freq[meta.gis.table$status == 'out']);

# Add sample sizes and calculate standard errors
meta.gis.table.num <- cbind(meta.gis.table, sample = meta.gis.table.df$Freq);
meta.gis.table.num$se <- sqrt(meta.gis.table.num$Freq * (1 - meta.gis.table.num$Freq) / meta.gis.table.num$sample);



### 2. TCGA-BRCA
brca.cnv.chr.new.gis.match <- brca.cnv.chr.new.gis[na.omit(match(outlier.symbol$brca, brca.cnv.chr.new.gis$Hugo_Symbol)), ];
rownames(brca.cnv.chr.new.gis.match) <- brca.cnv.chr.new.gis.match$Hugo_Symbol;
brca.cnv.chr.new.gis.match <- brca.cnv.chr.new.gis.match[, 3:ncol(brca.cnv.chr.new.gis.match)];
brca.cnv.chr.new.gis.match <- brca.cnv.chr.new.gis.match[, colnames(brca.cnv.chr.new.gis.match) %in% substr(colnames(outlier.patient.tag.01.brca), 1, 15)];
outlier.patient.tag.01.brca.symbol <- fpkm.tumor.symbol.filter.brca[rownames(outlier.patient.tag.01.brca), ]$Symbol;
brca.outlier.symbol.match <- rownames(outlier.patient.tag.01.brca)[match(rownames(brca.cnv.chr.new.gis.match), outlier.patient.tag.01.brca.symbol)]
outlier.patient.tag.01.brca.cnv.match <- outlier.patient.tag.01.brca[brca.outlier.symbol.match, match(colnames(brca.cnv.chr.new.gis.match), substr(colnames(outlier.patient.tag.01.brca), 1, 15))];

brca.outlier.sample.cnv.new.gis <- list();
brca.non.outlier.sample.cnv.new.gis <- list();
for (i in 1:nrow(brca.cnv.chr.new.gis.match)) {
    out.patient <- substr(colnames(outlier.patient.tag.01.brca.cnv.match)[outlier.patient.tag.01.brca.cnv.match[i, ] == 1], 1, 15);
    non.out.patient <- colnames(brca.cnv.chr.new.gis.match)[!(colnames(brca.cnv.chr.new.gis.match) %in% out.patient)];
    brca.outlier.sample.cnv.new.gis[[i]] <- brca.cnv.chr.new.gis.match[i, out.patient];
    brca.non.outlier.sample.cnv.new.gis[[i]] <- brca.cnv.chr.new.gis.match[i, non.out.patient];
    }

names(brca.outlier.sample.cnv.new.gis) <- rownames(brca.cnv.chr.new.gis.match);
names(brca.non.outlier.sample.cnv.new.gis) <- rownames(brca.cnv.chr.new.gis.match);

brca.outlier.cnv <- as.numeric(
    unlist(brca.outlier.sample.cnv.new.gis)
    );
brca.non.outlier.cnv <- as.numeric(
    unlist(brca.non.outlier.sample.cnv.new.gis)
    );

brca.gis.value <- data.frame(new.gis.value = c(brca.non.outlier.cnv, brca.outlier.cnv));
brca.gis.box <- data.frame(
    cbind(
        brca.gis.value$new.gis.value,
        c(
            rep('non', length(brca.non.outlier.cnv)),
            rep('out', length(brca.outlier.cnv))
            )
        )
    );
colnames(brca.gis.box) <- c('new.gis.value', 'status');
brca.gis.box[, 1] <- as.numeric(brca.gis.box[, 1]);

brca.gis.table.df <- data.frame(table(brca.gis.box));
brca.gis.table <- data.frame(table(brca.gis.box));

brca.gis.table$Freq[brca.gis.table$status == 'non'] <- brca.gis.table$Freq[brca.gis.table$status == 'non'] / sum(brca.gis.table$Freq[brca.gis.table$status == 'non']);
brca.gis.table$Freq[brca.gis.table$status == 'out'] <- brca.gis.table$Freq[brca.gis.table$status == 'out'] / sum(brca.gis.table$Freq[brca.gis.table$status == 'out']);

brca.gis.table.num <- cbind(brca.gis.table, sample = brca.gis.table.df$Freq);
brca.gis.table.num$se <- sqrt(brca.gis.table.num$Freq * (1 - brca.gis.table.num$Freq) / brca.gis.table.num$sample);



### 3. ICGC BRCA-EU
# classify copy number using cut-off
icgc.cnv.chr.new.gis.raw.patient.only <- icgc.cnv.chr.new.gis.raw[, 4:ncol(icgc.cnv.chr.new.gis.raw)];
icgc.cnv.chr.new.gis.new <- NULL;
for (i in 1:ncol(icgc.cnv.chr.new.gis.raw.patient.only)) {
    cnv.class <- as.numeric(icgc.cnv.chr.new.gis.raw.patient.only[, i]);
    low.cut <- as.numeric(icgc.cnv.chr.new.gis.cutoff[i, ]$Low);
    high.cut <- as.numeric(icgc.cnv.chr.new.gis.cutoff[i, ]$High);

    condition1 <- cnv.class >= -0.15 & cnv.class <= 0.15
    condition2 <- cnv.class < -0.15 & cnv.class >= low.cut
    condition3 <- cnv.class < low.cut
    condition4 <- cnv.class > 0.15 & cnv.class <= high.cut
    condition5 <- cnv.class > high.cut

    cnv.class[condition1] <- 0
    cnv.class[condition2] <- -1
    cnv.class[condition3] <- -2
    cnv.class[condition4] <- 1
    cnv.class[condition5] <- 2

    icgc.cnv.chr.new.gis.new <- cbind(icgc.cnv.chr.new.gis.new, cnv.class);
    }

icgc.cnv.chr.new.gis.new <- data.frame(icgc.cnv.chr.new.gis.new);
colnames(icgc.cnv.chr.new.gis.new) <- colnames(icgc.cnv.chr.new.gis.raw.patient.only);
icgc.cnv.chr.new.gis <- icgc.cnv.chr.new.gis.new;
icgc.outlier.symbol <- fpkm.data.icgc$Name[as.numeric(rownames(outlier.patient.tag.01.icgc))];

icgc.cnv.chr.new.gis.match <- icgc.cnv.chr.new.gis[na.omit(match(icgc.outlier.symbol, sub('\\|.*', '', icgc.cnv.chr.new.gis.raw$Gene.Symbol))), ];
rownames(icgc.cnv.chr.new.gis.match) <- sub('\\|.*', '', icgc.cnv.chr.new.gis.raw$Gene.Symbol)[na.omit(match(icgc.outlier.symbol, sub('\\|.*', '', icgc.cnv.chr.new.gis.raw$Gene.Symbol)))];

col_names <- colnames(icgc.cnv.chr.new.gis.match)
col_names <- gsub('^PD', '', col_names)
col_names <- gsub('_2$', '', col_names)
col_names <- gsub('a2?$', '', col_names)
col_names <- gsub('b2?$', '', col_names)
col_names <- gsub('c2?$', '', col_names)
col_names.1 <- colnames(outlier.patient.tag.01.icgc)
col_names.1 <- gsub('^PR', '', col_names.1)
col_names.1 <- gsub('.RNA$', '', col_names.1)
col_names.1 <- gsub('.2$', '', col_names.1)
col_names.1 <- gsub('a2?$', '', col_names.1)
col_names.1 <- gsub('b2?$', '', col_names.1)
col_names.1 <- gsub('c2?$', '', col_names.1)
col_names.1[col_names.1 %in% '7201a3'] <- '7201';
icgc.cnv.chr.new.gis.match <- icgc.cnv.chr.new.gis.match[, col_names %in% col_names.1];

col_names.2 <- colnames(icgc.cnv.chr.new.gis.match)
col_names.2 <- gsub('^PD', '', col_names.2)
col_names.2 <- gsub('_2$', '', col_names.2)
col_names.2 <- gsub('a2?$', '', col_names.2)
col_names.2 <- gsub('b2?$', '', col_names.2)
col_names.2 <- gsub('c2?$', '', col_names.2)
icgc.outlier.symbol.match <- rownames(outlier.patient.tag.01.icgc)[match(rownames(icgc.cnv.chr.new.gis.match), icgc.outlier.symbol)]

outlier.patient.tag.01.icgc.cnv.match <- outlier.patient.tag.01.icgc[icgc.outlier.symbol.match, match(col_names.2, col_names.1)];
colnames(icgc.cnv.chr.new.gis.match) <- colnames(outlier.patient.tag.01.icgc.cnv.match);

icgc.cnv.chr.new.gis.col <- icgc.cnv.chr.new.gis[, col_names %in% col_names.1];
colnames(icgc.cnv.chr.new.gis.col) <- colnames(outlier.patient.tag.01.icgc.cnv.match);

icgc.cnv.chr.new.gis.fpkm.order.match <- icgc.cnv.chr.new.gis.col;

icgc.cnv.chr.new.gis.fpkm.order.match.chr <- icgc.cnv.chr.new.gis.raw$Cytoband;
icgc.cnv.chr.new.gis.fpkm.order.match.chr <- gsub('[pq].*', '', icgc.cnv.chr.new.gis.fpkm.order.match.chr);

icgc.outlier.sample.cnv.new.gis <- list();
icgc.non.outlier.sample.cnv.new.gis <- list();
for (i in 1:nrow(icgc.cnv.chr.new.gis.match)) {
    out.patient <- substr(colnames(outlier.patient.tag.01.icgc.cnv.match)[outlier.patient.tag.01.icgc.cnv.match[i, ] == 1], 1, 15);
    non.out.patient <- colnames(icgc.cnv.chr.new.gis.match)[!(colnames(icgc.cnv.chr.new.gis.match) %in% out.patient)];
    icgc.outlier.sample.cnv.new.gis[[i]] <- icgc.cnv.chr.new.gis.match[i, out.patient];
    icgc.non.outlier.sample.cnv.new.gis[[i]] <- icgc.cnv.chr.new.gis.match[i, non.out.patient];
    }
names(icgc.outlier.sample.cnv.new.gis) <- rownames(icgc.cnv.chr.new.gis.match);
names(icgc.non.outlier.sample.cnv.new.gis) <- rownames(icgc.cnv.chr.new.gis.match);

icgc.outlier.cnv <- as.numeric(
    unlist(icgc.outlier.sample.cnv.new.gis)
    );
icgc.non.outlier.cnv <- as.numeric(
    unlist(icgc.non.outlier.sample.cnv.new.gis)
    );

icgc.gis.value <- data.frame(new.gis.value = c(icgc.non.outlier.cnv, icgc.outlier.cnv));
icgc.gis.box <- data.frame(
    cbind(
        icgc.gis.value$new.gis.value,
        c(
            rep('non', length(icgc.non.outlier.cnv)),
            rep('out', length(icgc.outlier.cnv))
            )
        )
    );

colnames(icgc.gis.box) <- c('new.gis.value', 'status');
icgc.gis.box[, 1] <- as.numeric(icgc.gis.box[, 1]);
icgc.gis.table.df <- data.frame(table(icgc.gis.box));
icgc.gis.table <- data.frame(table(icgc.gis.box));

icgc.gis.table$Freq[icgc.gis.table$status == 'non'] <- icgc.gis.table$Freq[icgc.gis.table$status == 'non'] / sum(icgc.gis.table$Freq[icgc.gis.table$status == 'non']);
icgc.gis.table$Freq[icgc.gis.table$status == 'out'] <- icgc.gis.table$Freq[icgc.gis.table$status == 'out'] / sum(icgc.gis.table$Freq[icgc.gis.table$status == 'out']);

icgc.gis.table.num <- cbind(icgc.gis.table, sample = icgc.gis.table.df$Freq);
icgc.gis.table.num$se <- sqrt(icgc.gis.table.num$Freq * (1 - icgc.gis.table.num$Freq) / icgc.gis.table.num$sample);



### META ANALYSIS

library(metafor);

# Initialize vectors to store meta-analysis results
metafor.cnv.odd.all <- NULL;
metafor.cnv.lower.all <- NULL;
metafor.cnv.upper.all <- NULL;
metafor.cnv.p.all <- NULL;

# Perform meta-analysis for each of the five CNV categories
for (i in 1:5) {
    # Prepare the data for meta-analysis
    data <- data.frame(
        measure = c('OR', 'OR', 'OR'),
        ai = c(
            meta.gis.table.num$sample[meta.gis.table.num$status == 'out'][i],
            brca.gis.table.num$sample[brca.gis.table.num$status == 'out'][i],
            icgc.gis.table.num$sample[icgc.gis.table.num$status == 'out'][i]
            ),
        n1i = c(
            sum(meta.gis.table.num$sample[meta.gis.table.num$status == 'out']),
            sum(brca.gis.table.num$sample[brca.gis.table.num$status == 'out']),
            sum(icgc.gis.table.num$sample[icgc.gis.table.num$status == 'out'])
            ),
        ci = c(
            meta.gis.table.num$sample[meta.gis.table.num$status == 'non'][i],
            brca.gis.table.num$sample[brca.gis.table.num$status == 'non'][i],
            icgc.gis.table.num$sample[icgc.gis.table.num$status == 'non'][i]
            ),
        n2i = c(
            sum(meta.gis.table.num$sample[meta.gis.table.num$status == 'non']),
            sum(brca.gis.table.num$sample[brca.gis.table.num$status == 'non']),
            sum(icgc.gis.table.num$sample[icgc.gis.table.num$status == 'non'])
            )
        );

    # Calculate effect sizes and variances for meta-analysis
    es.data <- escalc(
        measure = 'OR',
        ai = ai,
        n1i = n1i,
        ci = ci,
        n2i = n2i,
        data = data
        );

    # Perform the meta-analysis using a random-effects model (DerSimonian-Laird method)
    meta.res <- rma.uni(yi, vi, data = es.data, method = 'DL');

    # Extract and store the results: odds ratio, confidence intervals, and p-value
    metafor.cnv.odd <- exp(meta.res$beta);
    metafor.cnv.lower <- exp(meta.res$ci.lb);
    metafor.cnv.upper <- exp(meta.res$ci.ub);
    metafor.cnv.p <- meta.res$pval;

    # Append the results to the corresponding vectors
    metafor.cnv.odd.all <- c(metafor.cnv.odd.all, metafor.cnv.odd);
    metafor.cnv.lower.all <- c(metafor.cnv.lower.all, metafor.cnv.lower);
    metafor.cnv.upper.all <- c(metafor.cnv.upper.all, metafor.cnv.upper);
    metafor.cnv.p.all <- c(metafor.cnv.p.all, metafor.cnv.p);
    }

# Prepare data for visualization (stacked bar plot)

# Adjust p-values for multiple testing using the Benjamini-Hochberg method (FDR)
metafor.cnv.p.all.fdr <- p.adjust(metafor.cnv.p.all, method = 'BH');



### PLOTTING RESULTS ############################################################
# - Use weighted mean
# Initialize an empty vector to store the weighted means
# Initialize vectors to store the weighted sums
brca.sum.vector <- sapply(unique(brca.gis.table.num$new.gis.value), function(x) {
    sum(brca.gis.table.num$sample[brca.gis.table.num$new.gis.value == x])
    })
brca.sum.vector <- rep(brca.sum.vector, each = 2);

meta.sum.vector <- sapply(unique(meta.gis.table.num$new.gis.value), function(x) {
    sum(meta.gis.table.num$sample[meta.gis.table.num$new.gis.value == x])
    })
meta.sum.vector <- rep(meta.sum.vector, each = 2);

icgc.sum.vector <- sapply(unique(icgc.gis.table.num$new.gis.value), function(x) {
    sum(icgc.gis.table.num$sample[icgc.gis.table.num$new.gis.value == x])
    })
icgc.sum.vector <- rep(icgc.sum.vector, each = 2);

# Initialize vector to store the weighted means
all.weighted.mean <- numeric(10);

# Compute weighted means for each status and GIS value combination
for (status in c('non', 'out')) {
    for (i in 1:5) {
        gis.value <- c(-2, -1, 0, 1, 2)[i];

        weights <- c(
            sum(brca.gis.table.num$sample[brca.gis.table.num$new.gis.value == gis.value & brca.gis.table.num$status == status]),
            sum(meta.gis.table.num$sample[meta.gis.table.num$new.gis.value == gis.value & meta.gis.table.num$status == status]),
            sum(icgc.gis.table.num$sample[icgc.gis.table.num$new.gis.value == gis.value & icgc.gis.table.num$status == status])
            );

        values <- c(
            brca.gis.table.num$Freq[brca.gis.table.num$new.gis.value == gis.value & brca.gis.table.num$status == status],
            meta.gis.table.num$Freq[meta.gis.table.num$new.gis.value == gis.value & meta.gis.table.num$status == status],
            icgc.gis.table.num$Freq[icgc.gis.table.num$new.gis.value == gis.value & icgc.gis.table.num$status == status]
            );

        # Compute the weighted average and store it
        index <- ifelse(status == 'non', i, i + 5);
        all.weighted.mean[index] <- weighted.mean(values, weights);
        }
    }


brca.cnv.chr.new.gis.fpkm <- brca.cnv.chr.new.gis[
    brca.cnv.chr.new.gis$Hugo_Symbol %in%
        fpkm.tumor.symbol.filter.brca$Symbol,
    ];
brca.cnv.chr.new.gis.fpkm.order <- brca.cnv.chr.new.gis.fpkm[
    na.omit(match(
        brca.cnv.chr.new.gis.fpkm.order.match.chr$gene_name,
        brca.cnv.chr.new.gis.fpkm$Hugo_Symbol
        )),
    ];
brca.cnv.chr.new.gis.fpkm.order.match <- brca.cnv.chr.new.gis.fpkm.order[
    ,
    c('Hugo_Symbol', colnames(brca.cnv.chr.new.gis.match)[
        colnames(brca.cnv.chr.new.gis.match) %in%
            substr(colnames(outlier.patient.tag.01.brca), 1, 15)
        ])
    ];


# Save these variables for later scripts
cache.multiple.computed.variables(c(
    'all.weighted.mean',
    'brca.gis.table.num',
    'meta.gis.table.num',
    'icgc.gis.table.num',
    'brca.cnv.chr.new.gis.fpkm.order.match',
    'brca.cnv.chr.new.gis.match',
    'icgc.cnv.chr.new.gis.fpkm.order.match',
    'icgc.cnv.chr.new.gis.fpkm.order.match.chr',
    'meta.cnv.chr.new.gis.match',
    'metafor.cnv.lower.all',
    'metafor.cnv.odd.all',
    'metafor.cnv.p.all',
    'metafor.cnv.upper.all',
    'metafor.cnv.p.all.fdr',
    'outlier.patient.tag.01.brca.cnv.match',
    'outlier.patient.tag.01.meta.cnv.match'
    ));

save.session.profile(file.path('output', 'copy.number.analysis.txt'));
