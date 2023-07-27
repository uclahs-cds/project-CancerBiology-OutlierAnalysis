### CONTINGENCY TABLES AND STATISTICAL ANALYSIS ####################################################
source('~/Desktop/BIGSUMMER.PROJ/META.KM.SUBTYPES.R');
source('~/Desktop/BIGSUMMER.PROJ/META.KM.OUTLIERS.R');
# install.packages('vcd');
library(vcd);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);

meta.clinic.sample.metadata <- meta.clinic.sample[1:4,];

meta.clinic.sample <- meta.clinic.sample[-c(1:4), ];

meta.clinic.sample.merged <- merge(
    meta.clinic.sample,
    outlier.counts.df,
    by = 0
    );

rownames(meta.clinic.sample.merged) <- meta.clinic.sample.merged$Row.names;

meta.clinic.sample.merged$Row.names <- NULL

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


