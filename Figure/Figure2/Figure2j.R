### HISTORY #####################################################################
# This script processes methylation data to analyze outlier and non-outlier 
# genes in  outier and non-outlier patients. 
# Date: 2024-08-14



# Load necessary libraries
library(stats);
library(BoutrosLab.plotting.general);

load('2024-06-02_brca_meta_methylation.RData');


unequal.quan <- rev(seq(0, 0.9, 0.1));

# 1. Analyze outlier probes in outlier patients
percent.beta.out.out.500.both <- NULL;
for (i in 1:nrow(two.outlier.promoter.symbol.sample.match.merge.filter.500.both)) {
    value.vector <- as.numeric(two.outlier.promoter.symbol.sample.match.merge.filter.500.both[i,]);
    percentile <- ecdf(value.vector);
    percent.value <- percentile(outlier.sample.me.two.500.both[[i]]);
    percent.value.mean <- na.omit(percent.value);
    percent.beta.out.out.500.both <- c(percent.beta.out.out.500.both, percent.value.mean);
    }

# 2. Analyze outlier probes in non-outlier patients
percent.beta.non.out.500.both <- NULL;
for (i in 1:nrow(two.outlier.promoter.symbol.sample.match.merge.filter.500.both)) {
    value.vector <- as.numeric(two.outlier.promoter.symbol.sample.match.merge.filter.500.both[i,]);
    percentile <- ecdf(value.vector);
    percent.value <- percentile(non.outlier.sample.me.two.500.both[[i]]);
    percent.value.mean <- na.omit(percent.value);
    percent.beta.non.out.500.both <- c(percent.beta.non.out.500.both, percent.value.mean);
    }

# 3. & 4. Analyze non-outlier probes in outlier and non-outlier patients
two.outlier.patient.status.merge.filter.sum.500.both <- apply(
    two.outlier.patient.status.merge.filter.500.both, 
    2, 
    function(x) { sum(na.omit(as.numeric(x))); }
    );

me.non.out.symbol.two.500.both <- unique(c(
    rownames(brca.outlier.non.promoter.symbol.sample.match.merge.500), 
    rownames(meta.outlier.non.promoter.symbol.sample.match.merge.500)
    ));

me.non.out.symbol.two.500.both <- me.non.out.symbol.two.500.both[
    !(me.non.out.symbol.two.500.both %in% me.out.symbol.two.500.both)
    ];

brca.outlier.non.promoter.symbol.sample.match.merge.500 <- data.frame(
    brca.outlier.non.promoter.symbol.sample.match.merge.500
    );

meta.outlier.non.promoter.symbol.sample.match.merge.500 <- data.frame(
    meta.outlier.non.promoter.symbol.sample.match.merge.500
    );

percent.beta.out.non.500.both <- list();
percent.beta.non.non.500.both <- list();

for (i in 1:length(me.non.out.symbol.two.500.both)) {
    i.symbol <- me.non.out.symbol.two.500.both[i];
    value.vector.brca <- as.numeric(
        brca.outlier.non.promoter.symbol.sample.match.merge.500[i.symbol,]
        );
    value.vector.meta <- as.numeric(
        meta.outlier.non.promoter.symbol.sample.match.merge.500[i.symbol,]
        );
    
    value.vector <- na.omit(c(value.vector.brca, value.vector.meta));
    percentile <- ecdf(value.vector);
    
    percent.value.out.non <- percentile(
        value.vector[two.outlier.patient.status.merge.filter.sum.500.both > 0]
        );
    percent.beta.out.non.500.both[[i]] <- na.omit(percent.value.out.non);
    
    percent.value.non.non <- percentile(
        value.vector[two.outlier.patient.status.merge.filter.sum.500.both == 0]
        );
    percent.beta.non.non.500.both[[i]] <- na.omit(percent.value.non.non);
    }

percent.beta.out.non.500.both <- unlist(percent.beta.out.non.500.both);
percent.beta.non.non.500.both <- unlist(percent.beta.non.non.500.both);

# Calculate density and create data frames
calculate_density <- function(data) {
    dens <- density(na.omit(data), from = 0, to = 1, n = 30);
    df <- as.data.frame(dens[c('x', 'y')]);
    return(list(df = df));
    }

dens_list <- list(
    out_out = calculate_density(percent.beta.out.out.500.both),
    non_out = calculate_density(percent.beta.non.out.500.both),
    out_non = calculate_density(percent.beta.out.non.500.both),
    non_non = calculate_density(percent.beta.non.non.500.both)
    );

# Prepare heatmap data
dens.merge.two.four.group.500.both <- cbind(
    group1 = dens_list$out_out$df$y,
    group2 = dens_list$non_out$df$y,
    group3 = dens_list$out_non$df$y,
    group4 = dens_list$non_non$df$y
    );

heat.df <- data.frame(dens.merge.two.four.group.500.both);
heat.df.rev <- heat.df[rev(seq(nrow(heat.df))),];

# Set up legend
legend.col <- list(
    legend = list(
        colours = c("#2166ac", "white", "#b2182b"),
        title = expression(underline('z-score')),
        labels = c(0, 1.5),
        size = 3,
        label.cex = 1,
        continuous = TRUE,
        height = 3
        )
    );

# Create heatmap
heat.out <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(heat.df.rev),
    clustering.method = 'none',
    colour.scheme = c("#2166ac", "white", "#b2182b"),
    col.colour = 'white',
    grid.row = FALSE,
    grid.col = TRUE,
    yaxis.tck = 0,
    xaxis.tck = 0,
    xaxis.lab = NULL,
    xlab.label = expression('Outlier Genes'),
    yaxis.cex = 1.2,
    xaxis.cex = 1.2,
    yaxis.rot = 0,
    xaxis.rot = 90,
    xlab.cex = 1.2,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    colour.centering.value = 0.75,
    at = seq(0, 1.5, 0.01),
    covariate.legend = legend.col,
    legend.cex = 1,
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );



# Save the scatter plot as a PDF
pdf(
    file = generate.filename(
        'merge_me_quantile', 
        'heatmap', 
        'pdf'
        ), 
    width = 6, 
    height = 4.5
    );
heat.out;
dev.off();

# Save the scatter plot as a PNG
png(
    file = generate.filename(
        'merge_me_quantile', 
        'heatmap', 
        'png'
        ), 
    width = 6, 
    height = 4.5,
    unit = 'in', 
    res = 1200
    );
heat.out;
dev.off();

