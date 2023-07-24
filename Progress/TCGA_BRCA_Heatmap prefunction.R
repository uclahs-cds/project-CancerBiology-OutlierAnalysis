### TCGA-BRCA Heatmapping #########################################################################

### PREAMBLE ######################################################################################
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(BoutrosLab.statistics.survival);
#loading TCGA_BRCA
load(
    file = '2023-07-07_TCGA_BRCA_Outlier.rda',
    );
#loading Metabric
load(
    file = '2023-07-07_Metabric_Outlier.rda',
    );


create.multipanelplot(
    filename = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Multiplot5',
        extension = 'tiff'
    ),
    layout.width = 2,
    layout.height = 6,
    height = 20,
    width = 25,
    resolution = 300,
    plot.objects = list(patient.outlier,zscore.map,gene.outlier,age.map,Tumor.stage.map,PFS.map,subtype.map),
    plot.objects.heights = c(2,10,1,1,1,1),
    plot.objects.widths = c(10,1),
    layout.skip = c(FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE),
    x.spacing = -1,
    y.spacing = c(-3,-2,-5,-5,-5,-5)
    );


zscore.map <- create.heatmap(
    x = Z.score,
    #clustering.method = 'diana',
    clustering.method = 'none',
    #plot.dendrograms = 'none',
    force.clustering = TRUE,
    #filename = generate.filename(
    #   project.stem = 'Cancer',
    #  file.core = 'Heatmean4',
    # extension = 'tiff'
    #),
    resolution = 300,
    height = 48,
    width = 48,
    print.colour.key = FALSE
    );



age.colours <- as.character(force.colour.scheme(c(t(outlier.brca.clinic['Diagnosis.Age'])), scheme = "age.categorical.prostate"))
brca.clinic.color <- cbind(outlier.brca.clinic,age.colours)
age.map <- create.heatmap(
    x = subset(brca.clinic.color,select = c('age.colours')),
    input.colours = TRUE,
    #filename = generate.filename(
    #   project.stem = 'cancer',
    #  file.core = 'AgeHeat',
    # extension = 'tiff'
    #),
    clustering.method = 'none',
    same.as.matrix = TRUE,
    force.grid.col = FALSE,
    col.colour = 'black',
    xaxis.tck = 0,
    xlab.label = ' ',
    xaxis.lab = rep('',nrow(brca.clinic.color)),
    xaxis.cex = 0,
    yaxis.tck = 0,
    ylab.label = ' ',
    yaxis.lab = rep('',nrow(brca.clinic.color)),
    legend.cex = 0,
    print.colour.key = FALSE,
    );

subtype.color <- replace(
    as.numeric(c(t(brca.clinic['Subtype']))),
    ('BRCA_Basal' == brca.clinic['Subtype']),
    'powderblue'
    );
subtype.color <- replace(
    subtype.color,
    ('BRCA_Her2' == brca.clinic['Subtype']),
    'lightskyblue'
    );
subtype.color <- replace(
    subtype.color,
    ('BRCA_LumA' == brca.clinic['Subtype']),
    'mediumpurple1'
    );
subtype.color <- replace(
    subtype.color,
    ('BRCA_LumB' == brca.clinic['Subtype']),
    'palevioletred2'
    );
subtype.color <- replace(
    subtype.color,
    ('BRCA_Normal' == brca.clinic['Subtype']),
    'yellow3'
    );
brca.clinic.color <- cbind(brca.clinic,subtype.color)
subtype.map <- create.heatmap(
    x = subset(brca.clinic.color,select = c('subtype.color')),
    input.colours = TRUE,
    clustering.method = 'none',
    same.as.matrix = TRUE,
    force.grid.col = FALSE,
    col.colour = 'black',
    xaxis.tck = 0,
    xlab.label = ' ',
    xaxis.lab = rep('',nrow(brca.clinic.color)),
    xaxis.cex = 0,
    yaxis.tck = 0,
    ylab.label = ' ',
    yaxis.lab = rep('',nrow(brca.clinic.color)),
    legend.cex = 0,
    print.colour.key = FALSE,
    );


PFS.color <- replace(
    as.numeric(c(t(brca.clinic['Progress.Free.Survival..Months.']))),
    (60 <= brca.clinic['Progress.Free.Survival..Months.']),
    'white'
    );
PFS.color <- replace(
    PFS.color,
    (60 > brca.clinic['Progress.Free.Survival..Months.'] & 24 < brca.clinic['Progress.Free.Survival..Months.']),
    'darkorchid1'
    );
PFS.color <- replace(
    PFS.color,
    (24 >= brca.clinic['Progress.Free.Survival..Months.'] & 12 < brca.clinic['Progress.Free.Survival..Months.']),
    'mediumorchid4'
    );
PFS.color <- replace(
    PFS.color,
    (12 >= brca.clinic['Progress.Free.Survival..Months.']),
    'purple4'
    );

brca.clinic.color <- cbind(brca.clinic,PFS.color)
PFS.map <- create.heatmap(
    x = subset(brca.clinic.color,select = c('PFS.color')),
    input.colours = TRUE,
    clustering.method = 'none',
    same.as.matrix = TRUE,
    force.grid.col = FALSE,
    col.colour = 'black',
    xaxis.tck = 0,
    xlab.label = ' ',
    xaxis.lab = rep('',nrow(brca.clinic.color)),
    xaxis.cex = 0,
    yaxis.tck = 0,
    ylab.label = ' ',
    yaxis.lab = rep('',nrow(brca.clinic.color)),
    legend.cex = 0,
    print.colour.key = FALSE,
    );


stage.colours <- as.character(force.colour.scheme(outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code, scheme = "stage"))
brca.clinic.color <- cbind(outlier.brca.clinic,stage.colours)
Tumor.stage.map <- create.heatmap(
    x = subset(brca.clinic.color,select = c('stage.colours')),
    #x = stage.color,
    input.colours = TRUE,
    clustering.method = 'none',
    same.as.matrix = TRUE,
    force.grid.col = FALSE,
    col.colour = 'black',
    colour.scheme = force.colour.scheme(scheme = 'stage'),
    xaxis.tck = 0,
    xlab.label = ' ',
    xaxis.lab = rep('',nrow(brca.clinic.color)),
    xaxis.cex = 0,
    yaxis.tck = 0,
    ylab.label = ' ',
    yaxis.lab = rep('',nrow(brca.clinic.color)),
    legend.cex = 0,
    print.colour.key = FALSE,
    height = 4,
    width = 24
    );

patient.outlier <- create.barplot(
    formula = Outlier.totals ~ Patient.ID,
    data = outlier.brca.clinic,
    xlab.label = '',
    xat = '',
    xaxis.tck = 0,
    ylab.label = '# Outliers',
    yaxis.tck = 0,
    yat = c(0,25,50),
    ylab.cex = 1.5,
    yaxis.cex = 1.5
    );
gene.outlier <- create.barplot(
    formula = Gene.ID ~ Outlier.totals,
    data = outlier.totals.by.gene,
    height = 30,
    ylab.label = '',
    yaxis.lab = rep('',nrow(outlier.totals.by.gene)),
    yat = 0,
    xat = c(0,2,4),
    xlab.top.y = 1,
    xaxis.tck = 0,
    xlab.label = '',
    xlab.top.label = '# Outlier Patients',
    plot.horizontal = T,
    xlab.top.cex = 1.5,
    xaxis.cex = 1.5
    );
