### HISTORY #####################################################################
# This script generates a scatter and density plot to visualize the distribution
# of outlier genes per patient across Luminal A and Luminal B subtypes,
# comparing mutated and non-mutated PIK3CA cases.
# Date: 2024-08-14

### DESCRIPTION #################################################################
# The script processes data for the number of outlier genes per patient across
# Luminal A and Luminal B breast cancer subtypes. It compares mutated and
# non-mutated PIK3CA cases, combining the data into a single data frame and
# creating scatter and density plots for visualization.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-10-08_Figure1_2_3_4_min_input.rda'));

load.multiple.computed.variables(c(
    'meta.mutation.driver.list.gene.vector.data.convert.na',
    'outlier.patient.tag.01.t.p.order.sum'
    ));

# Helper function to subset and organize data by subtype
process.subtype.data <- function(subtype.id) {
    subtype.data <- subtype.total.outlier.num[subtype.total.outlier.num$subtype == subtype.id, ];
    outlier.patients <- outlier.patient.tag.01.t.p.order.sum[rownames(na.omit(subtype.data))]
    mutation.data <- meta.mutation.driver.list.gene.vector.data.convert.na[, rownames(na.omit(subtype.data))]

    outlier.brca <- outlier.patients[substr(names(outlier.patients), 1, 4) == 'TCGA']
    outlier.meta <- outlier.patients[substr(names(outlier.patients), 1, 2) == 'MB']
    outlier.icgc <- outlier.patients[substr(names(outlier.patients), 1, 2) == 'PR']

    mutation.brca <- mutation.data[, substr(colnames(mutation.data), 1, 4) == 'TCGA']
    mutation.meta <- mutation.data[, substr(colnames(mutation.data), 1, 2) == 'MB']
    mutation.icgc <- mutation.data[, substr(colnames(mutation.data), 1, 2) == 'PR']

    list(
        outlier.patients.all = c(outlier.meta, outlier.brca, outlier.icgc),
        mutation.data_all = data.frame(mutation.meta, mutation.brca, mutation.icgc)
        )
    }

# Process data for Luminal A (luma) and Luminal B (lumb)
luma.data <- process.subtype.data(3)
lumb.data <- process.subtype.data(4)

# Define mutated and non-mutated PIK3CA cases
gene.alias <- 'PIK3CA'
luma.mut <- luma.data$outlier.patients.all[luma.data$mutation.data_all[gene.alias, ] == 'mutation']
luma.non <- luma.data$outlier.patients.all[luma.data$mutation.data_all[gene.alias, ] == 'normal']
lumb.mut <- lumb.data$outlier.patients.all[lumb.data$mutation.data_all[gene.alias, ] == 'mutation']
lumb.non <- lumb.data$outlier.patients.all[lumb.data$mutation.data_all[gene.alias, ] == 'normal']

# Perform statistical tests
wilcox.test(luma.mut, luma.non, alternative = 'two.sided', conf.int = TRUE)
wilcox.test(lumb.mut, lumb.non, alternative = 'two.sided', conf.int = TRUE)

# Create density plots
create.density.data <- function(values, bw.value) {
    density.data <- density(na.omit(values), bw = bw.value, from = 0, to = max(as.numeric(values)))
    return(as.data.frame(density.data[c('x', 'y')]))
    }

bw.value <- 0.6
dens.luma.mut.df <- create.density.data(luma.mut, bw.value)
dens.luma.non.df <- create.density.data(luma.non, bw.value)
dens.lumb.mut.df <- create.density.data(lumb.mut, bw.value)
dens.lumb.non.df <- create.density.data(lumb.non, bw.value)

# Combine density data
combine.density <- rbind(
    dens.luma.mut.df, dens.luma.non.df, dens.lumb.mut.df, dens.lumb.non.df
    )

combine.density.group <- cbind(
    combine.density,
    group = c(
        rep('a', nrow(dens.luma.mut.df)),
        rep('b', nrow(dens.luma.non.df)),
        rep('c', nrow(dens.lumb.mut.df)),
        rep('d', nrow(dens.lumb.non.df))
        )
    )

# Create scatter plot
mutation.density <- BoutrosLab.plotting.general::create.scatterplot(
    y ~ log2(x + 1),
    data = combine.density.group,
    type = 'l',
    xlimits = c(-0.05, 6.4),
    groups = combine.density.group$group,
    ylab.label = expression('Density'),
    xlab.label = expression('Number of outlier genes per patient'),
    ylimits = c(-0.002, 0.45),
    yat = seq(0, 3, 0.2),
    xat = seq(0, 10, 2),
    xaxis.lab = expression(2^0, 2^2, 2^4, 2^6),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.cex = 1,
    yaxis.cex = 1,
    cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.4,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    lwd = 3,
    lty = c('dashed', 'solid', 'dashed', 'solid'),
    col = c('red3', 'red3', 'dodgerblue2', 'dodgerblue2'),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = c('red3', 'red3', 'dodgerblue2', 'dodgerblue2'),
                        pch = 21,
                        cex = 2,
                        fill = c('red3', 'red3', 'dodgerblue2', 'dodgerblue2')
                        ),
                    text = list(
                        lab = c(
                            'Luminal A - PIK3CA mutated patients',
                            'Luminal A - PIK3CA non-mutated patients',
                            'Luminal B - PIK3CA mutated patients',
                            'Luminal B - PIK3CA non-mutated patients'
                            )
                        ),
                    padding.text = c(0, 5, 0, 5),
                    cex = 1
                    )
                ),
            x = 0.2,
            y = 0.97,
            draw = FALSE
            )
        ),
    main = expression('Number of outlier genes per patient')
    )

# Save plot and session profile
save.outlier.figure(
    mutation.density,
    c('Figure2f', 'PIK3CA', 'mutation', 'luma', 'lumb', 'scatter', 'density'),
    width = 6,
    height = 5.5
    )

save.session.profile(file.path('output', 'Figure2f.txt'))
