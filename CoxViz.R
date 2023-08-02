library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)
library(BoutrosLab.plotting.survival)
load('meta.cox.RData')

# Get the row names from 'meta_significant_genes'
selected_genes <- rownames(meta_significant_genes)

# Subset the 'cox_models_list' to include only models with names in 'selected_genes'
selected_cox_models <- cox.model.results[names(cox.model.results) %in% selected_genes]

# SCATTER PLOTS: x = LOG2(HR) , y = -LOG10(PVAL)

# Get the gene names in the order of meta_significant_genes

# Reorder selected_cox_models to match the order of gene names in meta_significant_genes
selected_cox_models <- selected_cox_models[rownames(meta_significant_genes)]

# Create an empty dataframe to store the results
result_df <- data.frame(GeneName = rownames(meta_significant_genes),
                        PValue = meta_significant_genes$PValue,
                        HR = NA,
                        stringsAsFactors = FALSE)

# Iterate over the rows of meta_significant_genes and fill in HR values from selected_cox_models
for (i in seq_along(result_df$GeneName)) {
    gene_name <- result_df$GeneName[i]
    if (gene_name %in% names(selected_cox_models)) {
        hr <- exp(coef(selected_cox_models[[gene_name]]$model))
        result_df$HR[i] <- hr
    }
}

# Set the row names of the new dataframe to the gene names
create.scatterplot(
    height = 15,
    width = 15,
    formula = -log10(PValue) ~ log2(HR),
    data = result_df,
    main = 'METABRIC: Gene-Wise, Univariate Cox Models',
    main.cex = 1.2,
    ylab.label = expression(log[10]("q-value")),
    ylab.cex = 1
)


