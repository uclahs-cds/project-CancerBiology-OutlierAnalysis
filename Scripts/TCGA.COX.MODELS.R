output_dir <- "tcga_cox_modeling_output";
dir.create(output_dir, showWarnings = FALSE);

library(BoutrosLab.plotting.survival);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);
library(BoutrosLab.statistics.survival);
library(survival);
load('/Users/amaanjsattar/Desktop/2023-07-07_TCGA_BRCA_Outlier.rda')

tcga.fpkm <- fpkm.tumor.symbol.filter
tcga.outlier.tags <- outlier.patient.tag.01.t.p.order
tcga.fpkm <- tcga.fpkm[rownames(tcga.outlier.tags), ]

# Subset the clinical DataFrame to keep only the desired columns
tcga.subset <- brca.clinic[, c('Overall.Survival..Months.', 'Overall.Survival.Status')]

# transpose tcga.fpkm so that patients are rows and genes are columns
tcga.fpkm <- t(tcga.fpkm)

extract_first_12_chars <- function(row_names) {
    substr(row_names, 1, 12)
}

tcga.fpkm.substring <- extract_first_12_chars(rownames(tcga.fpkm))
sum(duplicated(tcga.fpkm.substring))

# Find the indices of the duplicates in tcga.fpkm.substring
duplicate_indices <- which(duplicated(tcga.fpkm.substring) | duplicated(tcga.fpkm.substring, fromLast = TRUE))

# Get all the rows with the duplicates in tcga.fpkm
rows_with_duplicates <- tcga.fpkm[duplicate_indices, ]

# Identify the rows that don't end with "01A" and get their indices
indices_to_remove <- duplicate_indices[!grepl("01A$", rownames(tcga.fpkm[duplicate_indices, ]))]

# Subsetting rows_with_duplicates to show the rows that will be removed
rows_to_remove <- rows_with_duplicates[!grepl("01A$", rownames(rows_with_duplicates)), ]

# Remove rows with indices that don't end with "01A" from the original tcga.fpkm
tcga.fpkm.filtered <- tcga.fpkm[-indices_to_remove, ]

tcga.fpkm.substring <- extract_first_12_chars(rownames(tcga.fpkm.filtered))

# Find the index of the row in tcga.subset that doesn't match any value in tcga.fpkm.substring
index_not_matching <- which(!(rownames(tcga.subset) %in% tcga.fpkm.substring))

# Remove the rows not matching from tcga.subset
tcga.subset.filtered <- tcga.subset[-index_not_matching, ]

# Find the value in tcga.fpkm.substring that is not a row name in tcga.subset.filtered
value_not_matching <- tcga.fpkm.substring[!(tcga.fpkm.substring %in% rownames(tcga.subset.filtered))]
# Find the index of the row in tcga.fpkm.filtered that has the value "Symbol" in tcga.fpkm.substring
index_with_symbol <- which(tcga.fpkm.substring %in% "Symbol")

# Print the row in tcga.fpkm.filtered with the value "Symbol" in tcga.fpkm.substring
symbol <- (tcga.fpkm.filtered[index_with_symbol, ])

# Remove the row from the dataframe
tcga.fpkm.filtered <- tcga.fpkm.filtered[-index_with_symbol, ]

# Store tcga.fpkm.filtered row names separately
tcga.fpkm.names <- rownames(tcga.fpkm.filtered)

# Remove the "Symbol" row from tcga.fpkm.substring
tcga.fpkm.substring <- tcga.fpkm.substring[tcga.fpkm.substring != "Symbol"]

# Set the modified tcga.fpkm.substring as the new row names for tcga.fpkm.filtered
rownames(tcga.fpkm.filtered) <- tcga.fpkm.substring

# Reorder tcga.fpkm.filtered so that it matches the order of tcga.subset.filtered
tcga.fpkm.filtered <- tcga.fpkm.filtered[match(rownames(tcga.subset.filtered), rownames(tcga.fpkm.filtered)), ]

### FUNCTIONS (TAKEN FROM TCGA.FUNCTIONS.R) ########################################################

binary.survival <- function(brca.data, old.col.name, new.col.name) {
    brca.data[[new.col.name]] <- ifelse(
        grepl('^1:', brca.data[[old.col.name]]), 1, 0)
    return(brca.data)
};


create.surv <- function(brca.data, survival.time, survival.status) {
    surv.obj <- Surv(as.numeric(brca.data[[survival.time]]),
                     brca.data[[survival.status]])
    return(surv.obj)
};

####################################################################################################
tcga.subset.filtered <- binary.survival(
    tcga.subset.filtered, 'Overall.Survival.Status', 'OS.Status'
);

tcga.surv.cox <- create.surv(
    tcga.subset.filtered, 'Overall.Survival..Months.', 'OS.Status'
);

tcga.cox.results <- list()
tcga.genes <- colnames(tcga.fpkm.filtered)

for (gene in tcga.genes) {
    # Extract the FPKM values for the current gene (as a vector)
    fpkm.values <- log2(as.numeric(tcga.fpkm.filtered[, gene]) + 1)
    
    # Fit the Cox proportional hazards model using the 'coxph' function from the 'survival' package
    cox.model <- coxph(tcga.surv.cox ~ fpkm.values)
    
    # Store the model and its summary in the cox.model.results list
    tcga.cox.results[[gene]] <- list(model = cox.model, summary = summary(cox.model))
}

gene_names_vector <- c()
p_values_vector <- c()

# Loop through each gene in gene.names and extract the p-value from the summary
for (gene in tcga.genes) {
    # Access the summary object for the specific gene
    summary_for_gene <- tcga.cox.results[[gene]]$summary
    
    # Access the p-value for the gene
    p_value <- summary_for_gene$coefficients["fpkm.values", "Pr(>|z|)"]
    
    # Store the gene name and p-value in the corresponding vectors
    gene_names_vector <- c(gene_names_vector, gene)
    p_values_vector <- c(p_values_vector, p_value)
}

# Create a dataframe with p-value results
tcga_p_values <- data.frame(GeneName = gene_names_vector, PValue = p_values_vector)

# Sort p-values in ascending order
tcga_p_values <- tcga_p_values[order(tcga_p_values$PValue), ]

fdr_corrected_p_values <- p.adjust(tcga_p_values$PValue, method = 'BH')

tcga_p_values$FDR.Corrected <- fdr_corrected_p_values
tcga_significant_genes_unadjusted <- subset(tcga_p_values, PValue <= 0.05)
tcga_significant_genes <- subset(tcga_p_values, FDR.Corrected <= 0.1)

tcga_significant_genes$Ensembl <- sub("\\.\\d+$", "", tcga_significant_genes$GeneName)


# Load the biomaRt library
library(biomaRt)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# Assuming your TCGA gene names are in the column 'GeneName' of the dataframe tcga_significant_genes
tcga_significant_genes$Official_Name <- mapIds(org.Hs.eg.db, keys = tcga_significant_genes$Ensembl, column = "SYMBOL", keytype = "ENSEMBL")



