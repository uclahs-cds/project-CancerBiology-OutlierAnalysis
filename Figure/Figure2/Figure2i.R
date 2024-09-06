### HISTORY #####################################################################
# This script processes NGF and LRP4 gene expression and methylation data, and generates 
# a scatter plot to visualize the relationship between RNA abundance and DNA 
# methylation levels across different patients.
# The codes are connected to Figure 2h.
# Date: 2024-08-14


library(BoutrosLab.plotting.general);



# Two example genes: NGF, LRP4 - run separately
do.plot.2i <- function(i, filename) {

i.me <- two.outlier.promoter.symbol.sample.match.merge.filter.500[
    i, 
    , 
    drop = FALSE
    ];

i.patient <- two.outlier.patient.status.merge.filter.500[
    i, 
    ];

i.me.mean <- mean(
    as.numeric(i.me), 
    na.rm = TRUE
    );

i.me.patient <- i.me[
    , 
    i.patient == 1, 
    drop = FALSE
    ];

i.me.patient.non <- i.me[
    , 
    i.patient == 0, 
    drop = FALSE
    ];

i.me.patient.non.mean <- mean(
    as.numeric(i.me.patient.non), 
    na.rm = TRUE
    );

i.me.patient.mean <- mean(
    as.numeric(i.me.patient), 
    na.rm = TRUE
    );

# FPKM data processing
fpkm.i.brca <- fpkm.tumor.symbol.filter.brca[
    fpkm.tumor.symbol.filter.brca$Symbol %in% i, 
    -ncol(fpkm.tumor.symbol.filter.brca), 
    drop = FALSE
    ];

fpkm.i.meta <- fpkm.tumor.symbol.filter.meta.symbol[
    fpkm.tumor.symbol.filter.meta.symbol$Symbol %in% i, 
    -ncol(fpkm.tumor.symbol.filter.meta.symbol), 
    drop = FALSE
    ];

fpkm.i.two <- c(
    scale(as.numeric(fpkm.i.brca)), 
    scale(as.numeric(fpkm.i.meta))
    );

fpkm.i.two.df <- t(
    data.frame(fpkm.i.two)
    );

rownames(fpkm.i.two.df) <- i;

colnames(fpkm.i.two.df) <- c(
    colnames(fpkm.i.brca), 
    colnames(fpkm.i.meta)
    );

# Prepare scatter plot data
fpkm.i.order <- fpkm.i.two.df[
    , 
    colnames(two.outlier.patient.status.merge.filter), 
    drop = FALSE
    ];

fpkm.i.order <- fpkm.i.order[
    , 
    order(as.numeric(fpkm.i.order[1, ]), decreasing = TRUE), 
    drop = FALSE
    ];

i.me.order <- i.me[
    , 
    colnames(fpkm.i.order)
    ];

i.me.fpkm.scatter <- data.frame(
    fpkm = as.numeric(fpkm.i.order),
    me = as.numeric(i.me.order)
    );

rownames(i.me.fpkm.scatter) <- colnames(i.me.order);

i.patient.order <- i.patient[
    colnames(i.me.order)
    ];

dot.colours <- ifelse(
    i.patient.order == 1, 
    'red2', 
    'black'
    );


i.me.fpkm.scatter.rev <- i.me.fpkm.scatter[rev(seq(nrow(i.me.fpkm.scatter))),]
dot.colours.rev <- rev(dot.colours);

scatter.i <- create.scatterplot(
    formula = fpkm ~ me,
    data = i.me.fpkm.scatter.rev,
    col = dot.colours.rev,
    alpha = .6,
    xlimits = c(-0.06, 1.07),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    add.grid = TRUE,
	grid.colour = 'grey80',
    cex = 0.9,
    add.text = FALSE,
    text.cex = 1.1,
    text.x = i.me.fpkm.scatter[which(i.patient == 1),]$me,
    text.y = log2(i.me.fpkm.scatter[which(i.patient == 1),]$fpkm),
    text.labels = '*outlier patient',
    text.guess.labels = TRUE,
    text.guess.label.position = 45,
    text.guess.radius.factor = 1.5,
    text.fontface = 1,
    text.col = 'red2',
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.6,
    main = i,
    ylab.label = expression(paste('RNA abundance (z-score)')),
    xlab.label = expression(paste('DNA methylation (', beta, ' value)')), 
    type = c('p', 'r', 'g'), 
        legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = i.me.fpkm.scatter$fpkm,
                    y = i.me.fpkm.scatter$me,
                    label.items = c('spearman'),
                    alpha.background = 0,
                    key.cex = 1.1
                ) ),
            x = 0.75,
            y = 0.95,
            corner = c(0,1)
            )
        )
    );

# Save the plot as a PNG
write.plot(
    trellis.object = scatter.i,
    filename = filename,
    width = 6,
    height = 6,
    size.units = 'in',
    resolution = 1200
);

}

output.directory <- get0('output.directory', ifnotfound = 'figures');

do.plot.2i('NGF', file.path(output.directory, 'Figure_2_i_NGF.png'));
do.plot.2i('LRP4', file.path(output.directory, 'Figure_2_i_LRP4.png'));
