### HISTORY #####################################################################
# This script generates a boxplot for the gene SYK, comparing z-scores of drug 
# response in outlier vs non-outlier samples.
# Date: 2024-08-16

source(file.path(dirname(dirname(parent.frame(2)$ofile)), 'common_functions.R'));

i <- 'SYK';

drug.target <- sanger.zscore.drug.breast.match.info.df[sanger.zscore.drug.breast.match.info.df$repurposing_target %in% i,]$Drug.Name;

# Get patient status and drug target information
patient.status <- ccle.sample.outlier.status.overlap.na.samger.match.dup.filter[i,];
drug.target <- depmap.drug.info.match.sanger.dup[depmap.drug.info.match.sanger.dup$repurposing_target %in% i,];

# Get IC50 z-score values for the selected gene
ic50.value.i <- sanger.drug.match.dup.zscore[match(drug.target$Drug.Name, rownames(sanger.drug.match.dup.zscore)), colnames(patient.status), drop =FALSE];
#remove duplicated row
ic50.value.i <- ic50.value.i[c(1, 3, 4),];

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
    # jitter.factor = 10,
    main.cex = 1.5,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('z-score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1.1,
    xaxis.lab = rownames(ic50.value.i),
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2,0),
    xaxis.tck = c(0.2,0),
    xaxis.rot = 90,
    # ylimits = c(-1.7,2.5),
    ylimits = c(-2.05, 3.3),
    yat = seq(-2, 3, 1),
    sample.order = "none",
    points.pch = 16,
    points.cex = 1,
    points.col = dot.colours,
    lwd = 1.2,
    col = 'gold2',
    alpha = 0.25
    );

save.outlier.figure(
    i.drug.box.plot,
    c(i, 'box', 'drug'),
    width = 5,
    height = 6
    );
