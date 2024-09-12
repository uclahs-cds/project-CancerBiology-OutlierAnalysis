### HISTORY #####################################################################
# This script generates a heatmap to visualize the quantiles of protein 
# abundance for outlier genes using TCGA-BRCA CPTAC data. 
# This code is linked with the analysis in Figure 3a.
# Date: 2024-08-14


# Load necessary library
library(BoutrosLab.plotting.general);

# Calculate quantiles of protein abundance for outlier genes
percent.protein.cptac.quantile <- NULL;
for (i in 1:length(outlier.protein.cptac.list.no.p.na)) { 
    unequal.quan <- rev(seq(0, 0.9, 0.1));
    value.vector <- na.omit(as.numeric(unlist(non.outlier.protein.cptac.list.no.p.na[i])));
    non.value <- quantile(value.vector, p = unequal.quan);
    out.value <- as.numeric(unlist(outlier.protein.cptac.list.no.p.na[i]));
    all.value <- c(out.value, non.value);
    percent.protein.cptac.quantile <- rbind(percent.protein.cptac.quantile, all.value);
    }

# Calculate the mean of quantiles
percent.protein.cptac.quantile.mean <- apply(percent.protein.cptac.quantile, 2, mean);

# Prepare data for heatmap
heat.df <- t(data.frame(percent.protein.cptac.quantile));
rownames(heat.df) <- c('Outliers', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100');
heat.df.rev <- heat.df[rev(seq(nrow(heat.df))), ];

# Define legend for the heatmap
legend.col <- list(
    legend = list(
        colours = c("black", "white", "#b2182b"),
        title = expression(underline('z-score')), 
        labels = c(-3, 0, 3),
        size = 3,
        label.cex = 1, 
        continuous = TRUE,
        height = 3
        )
    );

# Generate the heatmap
heat.out <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(heat.df.rev),
    clustering.method = 'none',
    colour.scheme = c("black", "white", "#b2182b"), 
    col.colour = 'white',
    grid.row = FALSE, 
    grid.col = TRUE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    xaxis.lab = NULL,
    yaxis.lab = rev(c('Outliers', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100')),
    xlab.label = expression('Outlier Genes'),
    yaxis.cex = 1.2,
    xaxis.cex = 1.2,
    yaxis.rot = 0,
    xaxis.rot = 90,
    xlab.cex = 1.2,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    colour.centering.value = 0,
    at = seq(-3, 3, 0.01),
    covariate.legend = legend.col,
    legend.cex = 1,
    colourkey.cex = 1.3,
    print.colour.key = FALSE
    );  


# Save the heatmap as a PDF
pdf(
    file = generate.filename(
        'cptac', 
        'heatmap', 
        'pdf'
        ), 
    width = 6, 
    height = 4.5
    );
print(heat.out);
dev.off();

# Save the heatmap as a PNG
png(
    file = generate.filename(
        'cptac', 
        'heatmap', 
        'png'
        ), 
    width = 6, 
    height = 4.5,
    unit = 'in', 
    res = 1200
    );
print(heat.out);
dev.off();


