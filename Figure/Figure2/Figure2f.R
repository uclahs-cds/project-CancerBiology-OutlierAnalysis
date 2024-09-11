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

# PIK3CA density plot
#   1. LumA - mutated
#   2. LumA - non-mutated
#   3. LumB - mutated
#   4. LumB - non-mutated

# Luma subtype patient
subtype.total.outlier.num.luma <- subtype.total.outlier.num[subtype.total.outlier.num$subtype == 3, ];
outlier.patient.tag.01.t.p.order.sum.luma <- outlier.patient.tag.01.t.p.order.sum[rownames(na.omit(subtype.total.outlier.num.luma))];
meta.mutation.driver.list.gene.vector.data.convert.na.luma <- meta.mutation.driver.list.gene.vector.data.convert.na[, rownames(na.omit(subtype.total.outlier.num.luma))];

outlier.patient.tag.01.t.p.order.sum.luma.brca <- 
    outlier.patient.tag.01.t.p.order.sum.luma[substr(names(outlier.patient.tag.01.t.p.order.sum.luma), 1, 4) == 'TCGA'];
outlier.patient.tag.01.t.p.order.sum.luma.meta <- 
    outlier.patient.tag.01.t.p.order.sum.luma[substr(names(outlier.patient.tag.01.t.p.order.sum.luma), 1, 2) == 'MB'];
outlier.patient.tag.01.t.p.order.sum.luma.icgc <- 
    outlier.patient.tag.01.t.p.order.sum.luma[substr(names(outlier.patient.tag.01.t.p.order.sum.luma), 1, 2) == 'PR'];

meta.mutation.driver.list.gene.vector.data.convert.na.luma.brca <- 
    meta.mutation.driver.list.gene.vector.data.convert.na.luma[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma), 1, 4) == 'TCGA'];
meta.mutation.driver.list.gene.vector.data.convert.na.luma.meta <- 
    meta.mutation.driver.list.gene.vector.data.convert.na.luma[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma), 1, 2) == 'MB'];
meta.mutation.driver.list.gene.vector.data.convert.na.luma.icgc <- 
    meta.mutation.driver.list.gene.vector.data.convert.na.luma[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.luma), 1, 2) == 'PR'];

# Lumb subtype patient
subtype.total.outlier.num.lumb <- subtype.total.outlier.num[subtype.total.outlier.num$subtype == 4, ];
outlier.patient.tag.01.t.p.order.sum.lumb <- outlier.patient.tag.01.t.p.order.sum[rownames(na.omit(subtype.total.outlier.num.lumb))];
meta.mutation.driver.list.gene.vector.data.convert.na.lumb <- meta.mutation.driver.list.gene.vector.data.convert.na[, rownames(na.omit(subtype.total.outlier.num.lumb))];

outlier.patient.tag.01.t.p.order.sum.lumb.brca <- 
    outlier.patient.tag.01.t.p.order.sum.lumb[substr(names(outlier.patient.tag.01.t.p.order.sum.lumb), 1, 4) == 'TCGA'];
outlier.patient.tag.01.t.p.order.sum.lumb.meta <- 
    outlier.patient.tag.01.t.p.order.sum.lumb[substr(names(outlier.patient.tag.01.t.p.order.sum.lumb), 1, 2) == 'MB'];
outlier.patient.tag.01.t.p.order.sum.lumb.icgc <- 
    outlier.patient.tag.01.t.p.order.sum.lumb[substr(names(outlier.patient.tag.01.t.p.order.sum.lumb), 1, 2) == 'PR'];

meta.mutation.driver.list.gene.vector.data.convert.na.lumb.brca <- 
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb), 1, 4) == 'TCGA'];
meta.mutation.driver.list.gene.vector.data.convert.na.lumb.meta <- 
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb), 1, 2) == 'MB'];
meta.mutation.driver.list.gene.vector.data.convert.na.lumb.icgc <- 
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb[, substr(colnames(meta.mutation.driver.list.gene.vector.data.convert.na.lumb), 1, 2) == 'PR'];


meta.mutation.driver.list.gene.vector.data.convert.na.luma.all <- data.frame(
    meta.mutation.driver.list.gene.vector.data.convert.na.luma.meta,
    meta.mutation.driver.list.gene.vector.data.convert.na.luma.brca,
    meta.mutation.driver.list.gene.vector.data.convert.na.luma.icgc
    );

outlier.patient.tag.01.t.p.order.sum.luma.all <- c(
    outlier.patient.tag.01.t.p.order.sum.luma.meta,
    outlier.patient.tag.01.t.p.order.sum.luma.brca,
    outlier.patient.tag.01.t.p.order.sum.luma.icgc
    );

meta.mutation.driver.list.gene.vector.data.convert.na.lumb.all <- data.frame(
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb.meta,
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb.brca,
    meta.mutation.driver.list.gene.vector.data.convert.na.lumb.icgc
    );

outlier.patient.tag.01.t.p.order.sum.lumb.all <- c(
    outlier.patient.tag.01.t.p.order.sum.lumb.meta,
    outlier.patient.tag.01.t.p.order.sum.lumb.brca,
    outlier.patient.tag.01.t.p.order.sum.lumb.icgc
    );

i <- 'PIK3CA';

# Create scatter plot
i.luma.mut <- outlier.patient.tag.01.t.p.order.sum.luma.all[meta.mutation.driver.list.gene.vector.data.convert.na.luma.all[i,] == 'mutation'];
i.luma.non <- outlier.patient.tag.01.t.p.order.sum.luma.all[meta.mutation.driver.list.gene.vector.data.convert.na.luma.all[i,] == 'normal'];
i.lumb.mut <- outlier.patient.tag.01.t.p.order.sum.lumb.all[meta.mutation.driver.list.gene.vector.data.convert.na.lumb.all[i,] == 'mutation'];
i.lumb.non <- outlier.patient.tag.01.t.p.order.sum.lumb.all[meta.mutation.driver.list.gene.vector.data.convert.na.lumb.all[i,] == 'normal'];

wilcox.test(i.luma.mut, i.luma.non, alternative = "two.sided", conf.int = TRUE);
wilcox.test(i.lumb.mut, i.lumb.non, alternative = "two.sided", conf.int = TRUE);

bw.value <- 0.6;

dens.luma.mut <- density(
    na.omit(i.luma.mut), 
    bw = bw.value, 
    from = 0, 
    to = max(as.numeric(i.luma.mut))
    );

dens.luma.mut.df <- as.data.frame(
    dens.luma.mut[c('x', 'y')]
    );



dens.luma.non <- density(
    na.omit(i.luma.non), 
    bw = bw.value, 
    from = 0, 
    to = max(as.numeric(i.luma.non))
    );

dens.luma.non.df <- as.data.frame(
    dens.luma.non[c('x', 'y')]
    );



dens.lumb.mut <- density(
    na.omit(i.lumb.mut), 
    bw = bw.value, 
    from = 0, 
    to = max(as.numeric(i.lumb.mut))
    );

dens.lumb.mut.df <- as.data.frame(
    dens.lumb.mut[c('x', 'y')]
    );



dens.lumb.non <- density(
    na.omit(i.lumb.non), 
    bw = bw.value, 
    from = 0, 
    to = max(as.numeric(i.lumb.non))
    );

dens.lumb.non.df <- as.data.frame(
    dens.lumb.non[c('x', 'y')]
    );



combine.density.scatter.old <- rbind(
    dens.luma.mut.df,
    dens.luma.non.df,
    dens.lumb.mut.df,
    dens.lumb.non.df
    );

combine.density.scatter.group.old <- cbind(
    combine.density.scatter.old,
    group = c(
        rep('a', nrow(dens.luma.mut.df)),
        rep('b', nrow(dens.luma.non.df)),
        rep('c', nrow(dens.lumb.mut.df)),
        rep('d', nrow(dens.lumb.non.df))
        )
    );


mutation.density <- BoutrosLab.plotting.general::create.scatterplot(
    y ~ log2(x + 1),
    data = combine.density.scatter.group.old,
    type = 'l',
    xlimits = c(-0.05, 6.4),
    groups = combine.density.scatter.group.old$group,
    ylab.label = expression('Density'),
    xlab.label = expression('Number of outlier genes per patient'),
    ylimits = c(-0.002, 0.45),
    yat = seq(0, 3, 0.2),
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xat = seq(0, 10, 2),
    xaxis.lab = expression(2^0, 2^2, 2^4, 2^6),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.4,
    width = 10,
    height = 10,
    lwd = 3,
    lty = c("dashed", "solid", "dashed", "solid"),
    bandwidth.adjust = 0.1,
    col = c("red3", "red3", "dodgerblue2", "dodgerblue2"),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = c("red3", "red3", "dodgerblue2", "dodgerblue2"),
                        pch = 21,
                        cex = 2,
                        fill = c("red3", "red3", "dodgerblue2", "dodgerblue2")
                        ),
                    text = list(
                        lab = c('Luminal A - PIK3CA mutated patients', 
                                'Luminal A - PIK3CA non-mutated patients', 
                                'Luminal B - PIK3CA mutated patients',
                                'Luminal B - PIK3CA non-mutated patients')
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
    );




# Save the density plot as a PDF
pdf(
    file = generate.filename(
        'PIK3CA_mutation_luma_lumb', 
        'scatter_density', 
        'pdf'
        ), 
    width = 6, 
    height = 5.5
    );
print(mutation.density);
dev.off();

# Save the density plot as a PNG
png(
    file = generate.filename(
        'PIK3CA_mutation_luma_lumb', 
        'scatter_density', 
        'png'
        ), 
    width = 6, 
    height = 5.5,
    unit = 'in', 
    res = 1200
    );
print(mutation.density);
dev.off();


