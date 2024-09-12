### HISTORY #####################################################################
# Gene dependency score of example outlier genes. 
# Date: 2024-08-16

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));


# Filter for genes with dependency score difference > 0.5
gene.dependency.diff.matrix.05.overlap.minus.05 <- gene.dependency.diff.matrix.05.overlap[
    gene.dependency.diff.matrix.05.overlap$diff > 0.5, 
    ];
gene.dependency.diff.matrix.05.overlap.minus.05.order <- gene.dependency.diff.matrix.05.overlap.minus.05[
    order(gene.dependency.diff.matrix.05.overlap.minus.05$diff), 
    ];

# Prepare data for box plot
dependency.score.05.overlap.minus.05 <- gene.dependency.breast.t.num.match.05.na[
    as.numeric(rownames(gene.dependency.diff.matrix.05.overlap.minus.05.order)), 
    ];
outlier.status.05.overlap.minus.05 <- ccle.sample.outlier.status.overlap.na[
    as.numeric(rownames(gene.dependency.diff.matrix.05.overlap.minus.05.order)), 
    ];

# Create gene name box
gene.name.box.05 <- NULL;
for (i in 1:nrow(outlier.status.05.overlap.minus.05)) {
    gene.name <- rep(sub("\\..*", "", rownames(dependency.score.05.overlap.minus.05))[i], ncol(dependency.score.05.overlap.minus.05));
    gene.name.box.05 <- c(gene.name.box.05, gene.name);
}

# Combine dependency score and outlier status into a data frame
dependency.05.box <- data.frame(
    cbind(
        score = as.numeric(unlist(t(dependency.score.05.overlap.minus.05))),
        gene = gene.name.box.05,
        status = as.numeric(unlist(t(outlier.status.05.overlap.minus.05)))
        )
    );
dependency.05.box$score <- as.numeric(dependency.05.box$score);
dependency.05.box$status <- as.numeric(dependency.05.box$status);

# Filter out specific genes from the data
dependency.05.box.part <- dependency.05.box[!(
    dependency.05.box$gene %in% c('CASC3', 'TSEN54', 'MRPL21')
    ), ];

# Add specific genes to the box plot data
gene.dependency.diff.matrix.05.overlap.plus.05 <- gene.dependency.diff.matrix.05.overlap[
    gene.dependency.diff.matrix.05.overlap$symbol %in% c('TACC3', 'CCT2'), 
    ];
dependency.score.05.overlap.plus.05 <- gene.dependency.breast.t.num.match.05.na[
    rownames(gene.dependency.diff.matrix.05.overlap.plus.05), 
    ];
outlier.status.05.overlap.plus.05 <- ccle.sample.outlier.status.overlap.na[
    rownames(gene.dependency.diff.matrix.05.overlap.plus.05), 
    ];
dependency.05.box.plus <- data.frame(
    cbind(
        score = as.numeric(unlist(t(dependency.score.05.overlap.plus.05))),
        gene = c(
            rep(gene.dependency.diff.matrix.05.overlap.plus.05$symbol[1], ncol(dependency.score.05.overlap.minus.05)),
            rep(gene.dependency.diff.matrix.05.overlap.plus.05$symbol[2], ncol(dependency.score.05.overlap.minus.05))
        ),
        status = as.numeric(unlist(t(outlier.status.05.overlap.plus.05)))
        )
    );

# Combine the additional data into the main data frame
dependency.05.box.part <- rbind(dependency.05.box.part, dependency.05.box.plus);

# Convert score and status to numeric
dependency.05.box.part$score <- as.numeric(dependency.05.box.part$score);
dependency.05.box.part$status <- as.numeric(dependency.05.box.part$status);

# Further filter out specific genes from the data
dependency.05.box.part.4 <- dependency.05.box.part[!(
    dependency.05.box.part$gene %in% c('CCT2', 'TACC3', 'MSL1', 'RTN4IP1')
    ), ];

# Set colors for the box plot
dot.colours <- rep('grey70', nrow(dependency.05.box.part.4));
dot.colours[dependency.05.box.part.4$status == 1] <- 'red2';

# Establish an arbitrary but consistent random seed for plotting consistency
set.seed(sum(utf8ToInt('Figure4f')));

# Create the box plot
dependency.05.box.plot <- BoutrosLab.plotting.general::create.boxplot(
    formula = score ~ gene,
    data = dependency.05.box.part.4,
    main = expression('Gene dependency score of outlier genes'),
    outlier = TRUE,
    add.stripplot = TRUE,
    add.rectangle = TRUE,
    xleft.rectangle = seq(0.5, 14.5, 2),
    xright.rectangle = seq(1.5, 15.5, 2),
    ybottom.rectangle = -3,
    ytop.rectangle = 5,
    col.rectangle = "grey",
    alpha.rectangle = 0.25,
    main.cex = 1.5,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Gene dependency score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xaxis.rot = 90,
    sample.order = "increasing",
    points.pch = 16,
    points.cex = 1,
    points.col = dot.colours,
    lwd = 1.2,
    col = c('gold2'),
    alpha = 0.25
    );


save.outlier.figure(
    dependency.05.box.plot,
    c('gene', 'dependency', 'example', 'box'),
    width = 6,
    height = 5
    );
