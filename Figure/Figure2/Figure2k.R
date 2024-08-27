### HISTORY #####################################################################
# This script processes gene expression data to compare outlier and non-outlier 
# patients across tumor and normal samples. It generates a boxplot to visualize 
# the mean beta values and performs a Kruskal-Wallis test for statistical analysis.
# Date: 2024-08-14

# Load necessary library
library(BoutrosLab.plotting.general);



# Define functions
calculate.mean <- function(data, symbol, patients) {
    mean(na.omit(as.numeric(data[symbol, patients])));
    }

process.gene <- function(symbol, data) {
    # 1. Outlier patients - tumor
    out.tumor.patients <- colnames(two.outlier.patient.status.merge.filter.500)[which(two.outlier.patient.status.merge.filter.500[symbol,] == 1)];
    out.tumor.mean <- calculate.mean(two.outlier.promoter.symbol.sample.match.merge.filter.500, symbol, out.tumor.patients);
    
    # 2. Outlier patients - normal
    out.normal.patients.meta <- gsub("_", ".", sample.info.norm$samp[sample.info.norm$matched_tumor %in% gsub("\\.", "_", out.tumor.patients)]);
    out.normal.patients.brca <- colnames(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500)[substr(colnames(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500), 1, 12) %in% substr(out.tumor.patients, 1, 12)];
    
    out.normal.meta <- calculate.mean(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500, symbol, out.normal.patients.meta);
    out.normal.brca <- calculate.mean(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500, symbol, out.normal.patients.brca);
    out.normal.mean <- mean(c(out.normal.meta, out.normal.brca));
    
    # 3. Non-outlier patients - tumor
    non.out.tumor.mean <- calculate.mean(two.outlier.promoter.symbol.sample.match.merge.filter.500, symbol, setdiff(colnames(two.outlier.promoter.symbol.sample.match.merge.filter.500), out.tumor.patients));
    
    # 4. Non-outlier patients - normal
    non.out.normal.mean <- calculate.mean(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500, symbol, setdiff(colnames(two.outlier.promoter.symbol.sample.normal.match.merge.filter.500), c(out.normal.patients.meta, out.normal.patients.brca)));
    
    c(out.tumor.mean, out.normal.mean, non.out.tumor.mean, non.out.normal.mean);
    }

# Main analysis
# divide into outlier and non-outlier

outlier.sample.me.two.500 <- list();
non.outlier.sample.me.two.500 <- list();
for (i in 1:nrow(two.outlier.promoter.symbol.sample.match.merge.filter.500)) {
    
    # symbol.name <- rownames(two.outlier.promoter.symbol.sample.match.merge.filter);
    outlier.patient.gene <- two.outlier.patient.status.merge.filter.500[i,];
    # outlier.patient.gene.na <- na.omit(outlier.patient.gene);
    outlier.me <- two.outlier.promoter.symbol.sample.match.merge.filter.500[i,][outlier.patient.gene == 1];
    non.outlier.me <- two.outlier.promoter.symbol.sample.match.merge.filter.500[i,][outlier.patient.gene == 0];
    
    outlier.sample.me.two.500[[i]] <- outlier.me;
    non.outlier.sample.me.two.500[[i]] <- non.outlier.me;
    
    }

outlier.sample.me.two.unlist.500 <- as.numeric(unlist(outlier.sample.me.two.500));
non.outlier.sample.me.two.unlist.500 <- as.numeric(unlist(non.outlier.sample.me.two.500));


outlier.sample.me.two.unlist.mean.500 <- lapply(outlier.sample.me.two.500, function(x) {mean(na.omit(as.numeric(x)))});
non.outlier.sample.me.two.unlist.mean.500 <- lapply(non.outlier.sample.me.two.500, function(x) {mean(na.omit(as.numeric(x)))});

mean.beta.merge.two.500 <- apply(two.outlier.promoter.symbol.sample.match.merge.filter.500, 1, function(x) {mean(na.omit(as.numeric(x)))});
minus.beta.merge.two.500 <- as.numeric(outlier.sample.me.two.unlist.mean.500) - as.numeric(non.outlier.sample.me.two.unlist.mean.500);
mean.minus.ma.merge.two.500 <- data.frame(cbind(as.numeric(mean.beta.merge.two.500), as.numeric(minus.beta.merge.two.500), rownames(two.outlier.promoter.symbol.sample.match.merge.filter.500)));
mean.minus.ma.merge.two.500[,1] <- as.numeric(mean.minus.ma.merge.two.500[,1]);
mean.minus.ma.merge.two.500[,2] <- as.numeric(mean.minus.ma.merge.two.500[,2]);
colnames(mean.minus.ma.merge.two.500) <- c('mean.beta', 'minus.beta', 'Symbol');

two.out.non.tumor.normal.gene.500 <- mean.minus.ma.merge.two.500$Symbol[mean.minus.ma.merge.two.500$Symbol %in% rownames(normal.tumor.beta.comparison.two.minus.order.500)];

results <- lapply(two.out.non.tumor.normal.gene.500, process.gene);
two.out.non.tumor.normal.gene.value.mean.500 <- do.call(rbind, results);
rownames(two.out.non.tumor.normal.gene.value.mean.500) <- two.out.non.tumor.normal.gene.500;
colnames(two.out.non.tumor.normal.gene.value.mean.500) <- c('out.tu', 'out.nor', 'non.tu', 'non.nor');

two.out.non.tumor.normal.gene.value.mean.na.500 <- na.omit(two.out.non.tumor.normal.gene.value.mean.500);

# boxplot data
box.data <- data.frame(
    value = as.vector(two.out.non.tumor.normal.gene.value.mean.na.500),
    sample = rep(c('a', 'b', 'c', 'd'), each = nrow(two.out.non.tumor.normal.gene.value.mean.na.500))
    );

# Kruskal-Wallis test
kruskal.result.tumor.na <- kruskal.test(value ~ sample, data = box.data);
text.pvalue.tumor.na <- display.statistical.result(
    x = kruskal.result.tumor.na$p.value,
    statistic.type = 'p',
    symbol = ' = '
    );

# Create boxplot
tumor.normal.box.plot <- BoutrosLab.plotting.general::create.boxplot(
    formula = value ~ sample,
    data = box.data,
    xlab.label = NULL,
    ylab.label = expression(paste('Mean of ', beta, ' value')),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    outliers = FALSE,
    ylimits = c(-0.05, 1.1),
    yat = seq(0, 1, 0.2),
    key = list(
        text = list(lab = text.pvalue.tumor.na, cex = 1),
        x = 0.55,
        y = 0.95
        ),
    add.stripplot = FALSE,
    points.pch = 1,
    points.cex = 0.8,
    points.col = 'grey50',
    add.rectangle = TRUE,
    xleft.rectangle = c(1.5, 3.5),
    xright.rectangle = c(2.5, 5),
    ybottom.rectangle = -200,
    ytop.rectangle = 200,
    col.rectangle = "grey",
    alpha.rectangle = 0.25,
    lwd = 1.2,
    col = c('red3', 'gold3', 'dodgerblue4', 'darkgreen'),
    alpha = 0.4
    );



# Save the box plot as a PDF
pdf(
    file = generate.filename(
        'merge_tumour_normal_me', 
        'box', 
        'pdf'
        ), 
    width = 4.5, 
    height = 5
    );
tumor.normal.box.plot;
dev.off();

# Save the box plot as a PNG
png(
    file = generate.filename(
        'merge_tumour_normal_me', 
        'box', 
        'png'
        ), 
    width = 4.5, 
    height = 5,
    unit = 'in', 
    res = 1200
    );
tumor.normal.box.plot;
dev.off();


