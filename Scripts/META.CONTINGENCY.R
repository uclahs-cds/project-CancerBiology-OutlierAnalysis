### CONTINGENCY TABLES AND STATISTICAL ANALYSIS ####################################################
# RESULTS LOG: her2.stats, histgrade.spearman

source('~/Desktop/BIGSUMMER.PROJ/META.KM.SUBTYPES.R');
source('~/Desktop/BIGSUMMER.PROJ/META.KM.OUTLIERS.R');
# install.packages('vcd');
library(DescTools)
library(vcd);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);
library(gtsummary)
library(tidyverse)
meta.clinic.sample.metadata <- meta.clinic.sample[1:4,];

meta.clinic.sample <- meta.clinic.sample[-c(1:4), ];

meta.clinic.sample.merged <- merge(
    meta.clinic.sample,
    outlier.counts.df,
    by = 0
    );

rownames(meta.clinic.sample.merged) <- meta.clinic.sample.merged$Row.names;

meta.clinic.sample.merged$Row.names <- NULL

meta.clinic.sample.merged <- merge(
    meta.clinic.sample.merged,
    meta.clinic.patient.rm,
    by = 0
)


### FUNCTION #######################################################################################
# Create contingency table
create.contingency <- function(
        data, rows, cols) {
    cont <- table(data[[rows]], data[[cols]])
    rownames(cont) <- c('0 Outliers', '1 Outlier', '2 Outliers', '3+ Outliers')
    return(cont)
};
### HER2 STATUS AND OUTLIERS #######################################################################
# Remove NA Values
her2.subset <- subset(
    meta.clinic.sample.merged,
    HER2.Status %in% c('Negative', 'Positive')
);

# Acquire only relevant columns
her2.subset <- her2.subset[, c('HER2.Status', 'Outliers')];

# Assuming you have a two-column dataframe 'my_data' with columns 'Outlier_Group' and 'HER2_Status'


# Create the contingency table
her2.contingency <- table(her2.subset$Outliers, her2.subset$HER2.Status);

# Add proper row and column labels
rownames(her2.contingency) <- c("0 Outliers", "1 Outlier", "2 Outliers", "3+ Outliers");
colnames(her2.contingency) <- c("HER2 Negative", "HER2 Positive");

# Add row and column totalsL margin = 1 for row totals, margin = 2 for column totals
her2.contingency.totals <- addmargins(
    her2.contingency, 
    margin = 1
    );

her2.contingency.totals <- addmargins(
    her2.contingency, 
    margin = 2
    );

# Edit row and column names
rownames(her2.contingency.totals)[rownames(her2.contingency.totals) == 'Sum'] <- 'Total';

colnames(her2.contingency.totals)[colnames(her2.contingency.totals) == 'Sum'] <- 'Total';

# Conduct Chi-Squared Test for Independence (Caution: Don't use contingency table with totals)
her2.x.squared <- chisq.test(her2.contingency)$statistic

her2.p.value <- chisq.test(her2.contingency)$p.value

her2.cramer <- assocstats(her2.contingency)$cramer

her2.stats <- c(paste('X-Squared: ', her2.x.squared),
                paste('P-Value: ', her2.p.value),
                paste('Cramér\'s V:', her2.cramer)
    );

### INTERPRETATION: There is a statistically significant association between outlier gene count and 
    # HER2 status.  There is a small-to-moderate effect size indicating the strength of the association.

### #######################################################################
meta.clinic.sample.merged$Neoplasm.Histologic.Grade <- as.numeric(
    meta.clinic.sample.merged$Neoplasm.Histologic.Grade
    );

outlier_ranks <- c('0' = 0, '1' = 1, '2' = 2, '3+' = 3)
meta.clinic.sample.merged$Outliers.Ordinal <- outlier_ranks[as.character(meta.clinic.sample.merged$Outliers)]

histgrade.subset <- meta.clinic.sample.merged[, c('Neoplasm.Histologic.Grade', 'Outliers.Ordinal')];

histgrade.subset <- na.omit(histgrade.subset);

histgrade.spearman <- cor.test(
    histgrade.subset$Outliers.Ordinal,
    histgrade.subset$Neoplasm.Histologic.Grade, 
    method = "spearman");

histgrade.cont <- create.contingency(
    histgrade.subset, 
    'Outliers.Ordinal', 
    'Neoplasm.Histologic.Grade'
    );

chisq.test(histgrade.cont)

# Based on the results, there is a statistically significant positive correlation 
# (p-value < 2.2e-16, rho = 0.2208) between the 'Outliers.Ordinal' and 'Neoplasm.Histologic.Grade' 
# variables in the dataset.This indicates that 
# as the number of outliers increases (0, 1, 2, 3),
# the tumor histologic grade tends to increase (1, 2, 3) as well.


tbl.subset <- meta.clinic.sample.merged[, c('HER2.Status',
                                            'Neoplasm.Histologic.Grade',
                                            'Oncotree.Code', 'PR.Status',
                                            'Tumor.Stage', 'Outliers', 
                                            'ER.Status', 
                                            'Cancer.Type.Detailed',
                                            'Cellularity',
                                            'Pam50...Claudin.low.subtype',
                                            'Overall.Survival.Status',
                                            'X3.Gene.classifier.subtype',
                                            'Patient.s.Vital.Status',
                                            'Primary.Tumor.Laterality',
                                            'Tumor.Other.Histologic.Subtype',
                                            'Relapse.Free.Status'
                                            )]


tbl.subset[tbl.subset == ''] <- NA

# Create the summary table with the 'PBS' level excluded
summary_table <- tbl_summary(tbl.subset,
                             by = Outliers,
                             missing = 'no',
                             label = list(HER2.Status ~ 'HER2 Status',
                                          Neoplasm.Histologic.Grade ~ 'Neoplasm Histologic Grade',
                                          Oncotree.Code ~ 'Oncotree Code',
                                          PR.Status ~ 'Progesterone Receptor Status',
                                          Tumor.Stage ~ 'Tumor Stage',
                                          ER.Status ~ 'Estrogen Receptor Status',
                                          Cancer.Type.Detailed ~ 'Histologic Subtype',
                                          Pam50...Claudin.low.subtype ~ 'PAM50 Molecular Subtype',
                                          Overall.Survival.Status ~ 'Overall Survival Status',
                                          Relapse.Free.Status ~ 'Relapse.Free.Status'
                                          )
);

summary_table <- modify_caption(summary_table,
                                "**Contingency Table: Distribution of Clinical Variables by Outlier Gene Count**")
# Remove the 'PBS' level from the 'Oncotree.Code' variable


summary_table <- italicize_labels(summary_table)
# Display the modified summary table
summary_table


head(her2.subset)

# Custom function to calculate Cramer's V for the entire contingency table
my_cramer_v <- function(data, variable, by, ...) {
    tbl <- table(data[[variable]], data[[by]])
    rstatix::cramer_v(tbl)
}

# Calculate the overall Cramer's V value
cramer_value <- my_cramer_v(her2.subset, variable = "Outliers", by = "HER2.Status")

# Create the tbl_summary object with Chi-square p-values
tbl_cross_ex1 <- her2.subset %>% 
    tbl_summary(by = HER2.Status) %>% 
    add_p()

tbl_cross_ex1

tbl_cross_ex1 <- tbl_cross_ex1 %>% 
    add_stat(
        fns = list(Outliers ~ my_cramer_v)
    ) %>% 
    modify_header(
        add_stat_1 = "**Effect Size**"
    ) %>%
    modify_header(
        label = "**HER2 Status**"
    )

tbl_cross_ex1 <- modify_caption(tbl_cross_ex1,
                                "**Contingency Table and Test for Association: HER2 Status by Outlier Gene Count**")

tbl_cross_ex1







