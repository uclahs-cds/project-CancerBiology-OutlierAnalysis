### HISTORY ######################################################################
# This script benchmarks the OutSeekR algorithm across different conditions and
# visualizes the resulting AUPRC (Area under Precision-Recall Curve) values using
# heatmaps.

### DESCRIPTION ##################################################################
# The script evaluates the performance of OutSeekR under varying simulation
# settings (sample sizes, gene numbers, XEG prevalence, etc.) and generates
# comparative heatmaps to illustrate the resulting aUPRC scores. It constructs
# multi-panel visualizations, including legends and layout adjustments, to
# provide a clear benchmarking overview across conditions.

### PREAMBLE #####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);


### DATA PREPARATION ############################################################

# Color code
sub.right.2.col <- c('#efedf5', '#bdbddb', '#9f9ac7', '#6b52a2', '#422977');
num.simul.sam.col <- c('#dde1b2', '#84b186', '#47875d');
num.simul.gene.col <- c('#f5e5b2', '#f1a935');
sub.right.1.col <- colorRampPalette(c('white', 'deepskyblue4'))(10)[c(2, 4, 6, 8, 10)];
auprc.col <- c("#107090", "white","#b2402b");

# Figure legend
legend.sample.grob <- BoutrosLab.plotting.general:::legend.grob(
    list(
        legend = list(
            title = expression(underline('Number of simulated samples')),
            colours = num.simul.sam.col,
            labels = c('100', '500', '1000'),
            size = 2,
            label.cex = 0.8,
            continuous = FALSE,
            height = 3,
            padding.text = 2
            ),
        legend = list(
            title = expression(underline('Sample number with XEGs')),
            colours = sub.right.2.col,
            labels = c('1', '5', '10', '20', '50'),
            size = 2,
            label.cex = 0.8,
            continuous = FALSE,
            height = 3,
            padding.text = 2
            ),
        legend = list(
            title = expression(underline('XEG prevalence')),
            colours = sub.right.1.col,
            labels = c('0.1%', '1%', '5%', '10%', '20%'),
            size = 2,
            label.cex = 0.8,
            continuous = FALSE,
            height = 3,
            padding.text = 2
            ),
        legend = list(
            title = expression(underline('Number of simulated genes')),
            colours = num.simul.gene.col,
            labels = c('10,000', '20,000'),
            size = 2,
            label.cex = 0.8,
            continuous = FALSE,
            height = 3,
            padding.text = 2
            ),
        legend = list(
            title = expression(underline('z-score')),
            continuous = TRUE,
            colours = auprc.col,
            total.colours = 100,
            labels = c('0', '0.5', '1.0'),
            cex = 0.8,
            height = 3
            )),
    label.cex = 0.8,
    title.cex = 0.8,
    title.just = 'left',
    title.fontface = 'plain'
    );

# First main heatmap
pr.100.all.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = pr.100.all,
    clustering.method = 'none',
	colour.scheme = auprc.col,
	colour.alpha = 0.8,
    grid.row = TRUE, 
    grid.col = TRUE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.lab = NULL,
    xaxis.lab = NULL,
	axes.lwd = 0,
	row.colour = c('black', rep('white',9), 'black'),
    col.colour = c('black', 'white', 'black'),
    row.lwd = c(2, rep(1,9), 2),
    col.lwd = c(2, 1, 2),
    yaxis.cex = 0,
    xaxis.cex = 0,
    yaxis.rot = 0,
    xaxis.rot = 0,
    ylab.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    colour.centering.value = 0.5,
    at = seq(0, 1, 0.001),
	cell.text =  format(round(unlist(pr.100.all), digits = 2), nsmall = 2),
	text.cex = 0.8,
	text.fontface = 1,
	col.pos = rep(1:ncol(pr.100.all), each = nrow(pr.100.all)),
	row.pos = rep(nrow(pr.100.all):1, times = ncol(pr.100.all)),
	same.as.matrix = TRUE,
    print.colour.key = FALSE
    );
pr.100.all.hetmap;

# Second main heatmap
pr.500.all.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = pr.500.all,
    clustering.method = 'none',
	colour.scheme = auprc.col,
	colour.alpha = 0.8,
	row.colour = c('black', rep('white',19), 'black'),
    col.colour = c('black', 'white', 'black'),
    row.lwd = c(2, rep(1,19), 2),
    col.lwd = c(2, 1, 2),
    grid.row = TRUE, 
    grid.col = TRUE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.lab = NULL,
    xaxis.lab = NULL,
	axes.lwd = 0,
    yaxis.cex = 0,
    xaxis.cex = 0,
    yaxis.rot = 0,
    xaxis.rot = 0,
    ylab.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    colour.centering.value = 0.5,
    at = seq(0, 1, 0.001),
	cell.text =  format(round(unlist(pr.500.all), digits = 2), nsmall = 2),
	text.cex = 0.8,
	text.fontface = 1,
	col.pos = rep(1:ncol(pr.500.all), each = nrow(pr.500.all)),
	row.pos = rep(nrow(pr.500.all):1, times = ncol(pr.500.all)),
	same.as.matrix = TRUE,
    print.colour.key = FALSE
    );
pr.500.all.hetmap;

# Third main heatmap
pr.1000.all.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = pr.1000.all,
    clustering.method = 'none',
	colour.scheme = auprc.col,
	colour.alpha = 0.8,
    # total.colours = 27,
    # row.colour = 'black',
	row.colour = c('black', rep('white',24), 'black'),
    col.colour = c('black', 'white', 'black'),
    row.lwd = c(2, rep(1,24), 2),
    col.lwd = c(2, 1, 2),
    grid.row = TRUE, 
    grid.col = TRUE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.lab = NULL,
    xaxis.lab = NULL,
	axes.lwd = 2,
    yaxis.cex = 0,
    xaxis.cex = 0,
    yaxis.rot = 0,
    xaxis.rot = 0,
    ylab.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    colour.centering.value = 0.5,
    at = seq(0, 1, 0.001),
	cell.text =  format(round(unlist(pr.1000.all), digits = 2), nsmall = 2),
	text.cex = 0.8,
	text.fontface = 1,
	col.pos = rep(1:ncol(pr.1000.all), each = nrow(pr.1000.all)),
	row.pos = rep(nrow(pr.1000.all):1, times = ncol(pr.1000.all)),
	same.as.matrix = TRUE,
    print.colour.key = FALSE
    );
pr.1000.all.hetmap;


# creat empty plot
empty.data <- matrix(1, nrow = 1, ncol = 1);
rownames(empty.data) <- 1;
colnames(empty.data) <- 1;

empty.plot <- BoutrosLab.plotting.general:::create.heatmap(
    x = data.frame(empty.data),
    clustering.method = 'none',
    colour.scheme = c("white", "white"),  
    total.colours = 1,
    print.colour.key = FALSE,
    axes.lwd = 0,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xaxis.cex = 0,
    yaxis.cex = 0,
    grid.row = FALSE,
    grid.col = FALSE
    )
empty.plot;

# First column heatmap
heat.all.1 <-  BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(empty.plot, pr.100.all.hetmap),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(1, 2),
    panel.heights = c(1, 1.5),
    ylab.padding = -1,
    xlab.to.xaxis.padding = 0,
	axes.lwd = 2, 
    y.spacing = 0.1,
    main.cex = 0,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    yaxis.lab = NULL,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.tck = 0,
	remove.all.border.lines = TRUE
    );
heat.all.1;

# Second column heatmap
heat.all.2 <-  BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(empty.plot, pr.500.all.hetmap),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(1, 2),
    panel.heights = c(1, 0.25),
    ylab.padding = 0.5,
    xlab.to.xaxis.padding = -1.5,
	axes.lwd = 2, 
    y.spacing = 0.1,
    main.cex = 0,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.tck = 0,
	remove.all.border.lines = TRUE
    );
heat.all.2;

# Third column heatmap
heat.all.3 <-  BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(empty.plot, pr.1000.all.hetmap),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(1, 2),
    panel.heights = c(1, 0),
    ylab.padding = 0.5,
    xlab.to.xaxis.padding = -1.5,
	axes.lwd = 2, 
    y.spacing = 0.1,
    main.cex = 0,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.tck = 0,
	remove.all.border.lines = TRUE
    );
heat.all.3;



sub.right.1.col.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = data.frame(t(sub.right.1)),
    clustering.method = 'none',
    colour.scheme = sub.right.1.col,
    total.colours = 6,
    row.colour = 'black',
    col.colour = 'black',
    row.lwd = c(2, rep(1,24), 2),
    col.lwd = c(2, 2),
    grid.row = TRUE,
    grid.col = TRUE,
    yaxis.tck = 0,
    xaxis.tck = 0,
	axes.lwd = 2,
    xaxis.fontface = 1,
    xaxis.rot = 90,
    xaxis.cex = 1,
    print.colour.key = FALSE
    );
sub.right.1.col.hetmap;


sub.right.2.col.hetmap <- BoutrosLab.plotting.general:::create.heatmap(
    x = data.frame(t(sub.right.2)),
    clustering.method = 'none',
    colour.scheme = sub.right.2.col,
    total.colours = 6,
    row.colour = 'black',
    col.colour = 'black',
    row.lwd = c(2, rep(1,4), 2),
    col.lwd = c(2, 2),
    grid.row = TRUE,
    grid.col = TRUE,
    yaxis.tck = 0,
    xaxis.tck = 0,
	axes.lwd = 2,
    xaxis.fontface = 1,
    xaxis.rot = 90,
    xaxis.cex = 1,
    print.colour.key = FALSE
    );
sub.right.2.col.hetmap;




right.heatmap.1 <-  BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(empty.plot, empty.plot, sub.right.1.col.hetmap, sub.right.2.col.hetmap),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    layout.skip = c(FALSE,FALSE, FALSE, FALSE),
    plot.layout = c(2, 2),
    panel.heights = c(1, 0),
    panel.widths = c(1, 1),
    ylab.padding = 0.5,
    xlab.to.xaxis.padding = -1.5,
	axes.lwd = 1, 
    y.spacing = 0.1,
    x.spacing = 0,
    main.cex = 0,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.tck = 0,
	remove.all.border.lines = TRUE
    );
right.heatmap.1;

top.heatmap.1 <-  BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(sub.top.second.col.hetmap, sub.top.first.1.col.hetmap),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(1, 2),
    panel.heights = c(1, 1),
    ylab.padding = 0.5,
    xlab.to.xaxis.padding = -1.5,
	axes.lwd = 1, 
    y.spacing = -1,
    main.cex = 0,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.tck = 0,
	remove.all.border.lines = FALSE
    );
top.heatmap.1;

top.heatmap.2 <-  BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(sub.top.second.col.hetmap, sub.top.first.2.col.hetmap),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(1, 2),
    panel.heights = c(1, 1),
    ylab.padding = 0.5,
    xlab.to.xaxis.padding = -1.5,
	axes.lwd = 1,
    y.spacing = -1,
    main.cex = 0,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.tck = 0,
	remove.all.border.lines = FALSE
    );
top.heatmap.2;

top.heatmap.3 <-  BoutrosLab.plotting.general:::create.multiplot(
    plot.objects = list(sub.top.second.col.hetmap, sub.top.first.3.col.hetmap),
    x.relation = 'sliced',
    y.relation = 'sliced',
    main = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    layout.skip = c(FALSE, FALSE),
    plot.layout = c(1, 2),
    panel.heights = c(1, 1),
    ylab.padding = 0.5,
    xlab.to.xaxis.padding = -1.5,
	axes.lwd = 1, 
    y.spacing = -1,
    main.cex = 0,
    xaxis.cex = 0,
    xaxis.lab = NULL,
    yaxis.cex = 0,
    yaxis.tck = 0,
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.tck = 0,
	remove.all.border.lines = FALSE
    );
top.heatmap.3;

heat.combine <-  BoutrosLab.plotting.general:::create.multipanelplot(
    # main = NULL,
    plot.objects = list(
        top.heatmap.1,
        top.heatmap.2,
        top.heatmap.3,
        heat.all.1,
        heat.all.2,
        heat.all.3,
        right.heatmap.1
        ),
    layout.height = 2,
    layout.width = 4,
    plot.objects.heights = c(0.115, 1),
    plot.objects.widths = c(1, 1, 1, 0.53),
    right.legend.padding = 0,
	x.spacing = -1.5, 
	y.spacing = -0.7, 
    layout.skip = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    legend = list(right = list(fun = legend.sample.grob))
    );
heat.combine;


pdf(file = generate.filename('Benchmakring', 'heatmap', 'pdf'), width = 6.5, height = 11);
heat.combine;
dev.off();
    
png(file = generate.filename('Benchmakring', 'heatmap', 'png'), width = 6.5, height = 11, unit = 'in', res = 1200);
heat.combine;
dev.off();
