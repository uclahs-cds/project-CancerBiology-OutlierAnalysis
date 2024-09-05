### HISTORY #####################################################################
# This script generates a boxplot for the gene SYK, comparing z-scores of drug 
# response in outlier vs non-outlier samples.
# Date: 2024-08-16

i <- 'SYK';

drug.target <- sanger.zscore.drug.breast.match.info.df[sanger.zscore.drug.breast.match.info.df$repurposing_target %in% i,]$Drug.Name;

# Get patient status and drug target information
patient.status <- ccle.sample.outlier.status.overlap.na.samger.match.dup.filter[i,];
drug.target <- depmap.drug.info.match.sanger.dup[depmap.drug.info.match.sanger.dup$repurposing_target %in% i,];

# Get IC50 z-score values for the selected gene
ic50.value.i <- sanger.drug.match.dup.zscore[match(drug.target$Drug.Name, rownames(sanger.drug.match.dup.zscore)), colnames(patient.status), drop =FALSE];

# Create a data frame for plotting
i.drug.box <- data.frame(
    score = unlist(data.frame(t(ic50.value.i))),
    status = rep(unlist(patient.status), 3),
    group = c(
        rep('a',length(patient.status)),
        rep('b',length(patient.status)),
        rep('c',length(patient.status))
        )
    );

# Define colors for the points based on the patient status
dot.colours <- rep('grey70', nrow(i.drug.box));
dot.colours[i.drug.box$status ==1] <- 'green4';

# Generate the boxplot
i.drug.box.plot <- BoutrosLab.plotting.general::create.boxplot(
    formula = score ~ group,
    data = i.drug.box,
    main = expression('Gene effect score of outlier genes'),
    outlier = TRUE,
    add.stripplot = TRUE,
    jitter.factor = 10,
    main.cex = 1.5,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('z-score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.lab = i.drug,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    xaxis.rot = 90,
    ylimits = c(-1.7,2.5),
    yat = seq(-2,3,1),
    sample.order = "none",
    points.pch = 16,
    points.cex = 1,
    points.col = dot.colours,
    lwd = 1.2,
    col = 'gold2',
    alpha = 0.25
    );


# Save the box plot as a PDF
pdf(
    file = generate.filename(
        i, 
        'box_drug', 
        'pdf'
        ), 
    width = 5, 
    height = 6
    );
i.drug.box.plot;
dev.off();

# Save the box plot as a PNG
png(
    file = generate.filename(
        i, 
        'box_drug', 
        'png'
        ), 
    width = 5, 
    height = 6,
    unit = 'in', 
    res = 1200
    );
i.drug.box.plot;
dev.off();


