### HISTORY #####################################################################
# This script processes methylation data to analyze outlier and non-outlier 
# genes in  outier and non-outlier patients. 
# Date: 2024-08-14

### (DELETE) - Below analysis probably will take a while. I also included the intermediate variable such as percent.beta.out.out.500, percent.beta.non.out.500 and so on. 

# Load necessary libraries
library(stats);
library(BoutrosLab.plotting.general);

load('2024-06-02_brca_meta_methylation.RData');


unequal.quan <- rev(seq(0, 0.9, 0.1));

# Outlier patients on outlier probes
percent.beta.out.out.500 <- NULL;
for (i in 1:nrow(two.outlier.promoter.symbol.sample.match.merge.filter.500)) { 
    value.vector <- as.numeric(two.outlier.promoter.symbol.sample.match.merge.filter.500[i,]);
    percentile <- ecdf(as.numeric(value.vector));
    percent.value <- percentile(outlier.sample.me.two.500[[i]]);
    percent.value.mean <- na.omit(percent.value);
    percent.beta.out.out.500 <- c(percent.beta.out.out.500, percent.value.mean);
}

# Non-outlier patients on outlier probes
percent.beta.non.out.500 <- NULL;
for (i in 1:nrow(two.outlier.promoter.symbol.sample.match.merge.filter.500)) { 
    value.vector <- as.numeric(two.outlier.promoter.symbol.sample.match.merge.filter.500[i,]);
    percentile <- ecdf(as.numeric(value.vector));
    percent.value <- percentile(non.outlier.sample.me.two.500[[i]]);
    percent.value.mean <- na.omit(percent.value);
    percent.beta.non.out.500 <- c(percent.beta.non.out.500, percent.value.mean);
}

# Outlier and non-outlier patients on non-outlier probes
two.outlier.patient.status.merge.filter.sum.500 <- apply(two.outlier.patient.status.merge.filter.500, 2, function(x) { sum(na.omit(as.numeric(x))); });

beta.all.patient.non.probe.com.meta.symbol <- meta.me.com.symbol[rownames(beta.all.patient.non.probe.com.meta),]$name;

me.non.out.symbol.two.500 <- unique(c(rownames(brca.outlier.non.promoter.symbol.sample.match.merge.500), beta.all.patient.non.probe.com.meta.symbol));
me.non.out.symbol.two.500 <- me.non.out.symbol.two.500[!(me.non.out.symbol.two.500 %in% me.out.symbol.two.500)];
brca.outlier.non.promoter.symbol.sample.match.merge.500 <- data.frame(brca.outlier.non.promoter.symbol.sample.match.merge.500);

percent.beta.out.non.500 <- list();
percent.beta.non.non.500 <- list();
for (i in 1:length(me.non.out.symbol.two.500)) { 
    i.symbol <- me.non.out.symbol.two.500[i];
    value.vector.brca <- as.numeric(brca.outlier.non.promoter.symbol.sample.match.merge.500[i.symbol,]);
    
    i.symbol.meta <- rownames(meta.me.com.symbol)[meta.me.com.symbol$name %in% i.symbol];
    value.vector.meta.1 <- beta.all.patient.non.probe.com.meta[i.symbol.meta,];
    value.vector.meta.1.mea <- apply(value.vector.meta.1, 2, function(x) {mean(na.omit(as.numeric(x))); });
    value.vector.meta <- as.numeric(value.vector.meta.1.mea);
    value.vector <- na.omit(c(value.vector.brca, value.vector.meta));
    percentile <- ecdf(as.numeric(value.vector));
    
    percent.value.out.non <- percentile(value.vector[two.outlier.patient.status.merge.filter.sum.500 > 0]);
    percent.value.out.non.mean <- na.omit(percent.value.out.non);
    percent.beta.out.non.500[[i]] <- percent.value.out.non.mean;
    
    percent.value.non.non <- percentile(value.vector[two.outlier.patient.status.merge.filter.sum.500 == 0]);
    percent.value.non.non.mean <- na.omit(percent.value.non.non);
    percent.beta.non.non.500[[i]] <- percent.value.non.non.mean;
    }

percent.beta.non.non.500 <- unlist(percent.beta.non.non.500);
percent.beta.out.non.500 <- unlist(percent.beta.out.non.500);


# Create histograms and calculate percentages for each group
breaks <- seq(0, 1, length.out = 31);

create_histogram_df <- function(data) {
    hist_data <- hist(data, breaks = breaks, plot = FALSE)
    percentages <- hist_data$counts / sum(hist_data$counts) * 100
    data.frame(
        bin_start = hist_data$breaks[-length(hist_data$breaks)],
        bin_end = hist_data$breaks[-1],
        percentage = percentages
        );
    }

percent.beta.out.out.500.percentages.df <- create_histogram_df(percent.beta.out.out.500);
percent.beta.non.out.500.percentages.df <- create_histogram_df(percent.beta.non.out.500);
percent.beta.out.non.500.percentages.df <- create_histogram_df(percent.beta.out.non.500);
percent.beta.non.non.500.percentages.df <- create_histogram_df(percent.beta.non.non.500);

# Combine percentages for all groups
percent.merge.two.four.group.500 <- cbind(
    group1 = percent.beta.out.out.500.percentages.df$percentage,
    group2 = percent.beta.non.out.500.percentages.df$percentage,
    group3 = percent.beta.out.non.500.percentages.df$percentage,
    group4 = percent.beta.non.non.500.percentages.df$percentage
    );

# Prepare data for heatmap
heat.df <- data.frame(percent.merge.two.four.group.500);
heat.df.rev <- heat.df[rev(seq(nrow(heat.df))),];

# Set up legend
legend.col <- list(
    legend = list(
        colours = c("#2166ac", "white", "#b2182b"),
        title = expression(underline('Percentage')),
        labels = c(0, 10),
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
    colour.scheme = c("#2166ac", "white","#b2182b"), 
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
    colour.centering.value = 5,
    at = seq(0, 10, 0.01),
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

