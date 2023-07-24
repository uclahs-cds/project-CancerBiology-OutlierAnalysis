### TCGA-BRCA Heatmapping #########################################################################

### Preamble ######################################################################################
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

### fpkm.normilization ############################################################################
#Description
    #
#Input Variables
#data
    # the name of the data frame that contains non-normalized fpkm values
#Output Variables
    #
fpkm.normalization.zscore <- function(data) {
    number.rows <- nrow(data)
    number.cols <- ncol(data)
    fpkm.stats <- data.frame(
        Genes = row.names(data),
        Means.by.row = rep(NA,number.rows),
        SD.by.row = rep(NA,number.rows)
        );
    for (i in 1:number.rows) {
        fpkm.stats[i,2]<- mean(as.numeric(c(t(data[i,]))),na.rm = TRUE)
        fpkm.stats[i,3]<- sd(as.numeric(c(t(data[i,]))),na.rm = TRUE)
        }
    norm.fpkm <- data.frame(matrix(ncol = number.cols,nrow = number.rows,dimnames = list(row.names(data),colnames(data))))
    for (r in 1:number.rows) {
        for (c in 1:number.cols) {
            norm.fpkm[r,c] <- (as.numeric(data[r,c]) - fpkm.stats[r,2]) / fpkm.stats[r,3]
            print(r,c)
            }
        }
    result <- list(fpkm.stats,norm.fpkm)
    return(norm.fpkm)
    }

### fpkm.normilization ############################################################################
#Description
    #
#Input Variables
    #data:
    #the name of the data frame that contains non-normalized fpkm values
#Output Variables
    #median.norm.fpkm
    #data frame with normalized fpkm values utilizing the median and IQR

fpkm.robust.zscore <- function(data) {
    number.rows <- nrow(data)
    number.cols <- ncol(data)
    fpkm.stats <- data.frame(
        Genes = row.names(data),
        Median.by.row = rep(NA,number.rows),
        IQR.by.row = rep(NA,number.rows)
    );
    for (i in 1:number.rows) {
        fpkm.stats[i,2]<- median(as.numeric(c(t(data[i,]))),na.rm = TRUE)
        fpkm.stats[i,3]<- IQR(as.numeric(c(t(data[i,]))),na.rm = TRUE) * 0.7413
    }
    #empty df for results
    median.norm.fpkm <- data.frame(
        matrix(ncol = number.cols,
        nrow = number.rows,
        dimnames = list(row.names(data),colnames(data))))
    for (r in 1:number.rows) {
        for (c in 1:number.cols) {
            median.norm.fpkm[r,c] <- (as.numeric(data[r,c]) - fpkm.stats[r,2]) / fpkm.stats[r,3]
        }
    }
    #result <- list(fpkm.stats,median.norm.fpkm)
    return(median.norm.fpkm)
}

### Data Analysis #################################################################################
#saveRDS(median.normalized.fpkm,file = "Robust-Median-Zscore.RDS")
#saveRDS(fpkm.norm,file = 'Z-score.RDS')
fpkm.nonnorm <- fpkm.tumor.symbol.filter
### Fixing Column names in fpkm data to match clinical id
colnames(fpkm.nonnorm) <- substr(
    x = colnames(fpkm.nonnorm),
    start = 1,
    stop = nchar(colnames(fpkm.nonnorm))-4
    );
#replacing '.' with '-'
colnames(fpkm.nonnorm) <- gsub(
    pattern = '\\.',
    replacement = '-',
    x = colnames(fpkm.nonnorm)
    );


### TCGA-BRCA Heatmapping #########################################################################

###
###
###
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
### fpkm.normilization ############################################################################
#Description
#
#Input Variables
#data
#Output Variables
#
fpkm.normalization <- function(data) {
    number.rows <- nrow(data)
    number.cols <- ncol(data)
    fpkm.stats <- data.frame(
        Genes = row.names(data),
        Means.by.row = rep(NA,number.rows),
        SD.by.row = rep(NA,number.rows)
    );
    for (i in 1:number.rows) {
        fpkm.stats[i,2]<- mean(as.numeric(c(t(data[i,]))),na.rm = TRUE)
        fpkm.stats[i,3]<- sd(as.numeric(c(t(data[i,]))),na.rm = TRUE)
    }
    norm.fpkm <- data.frame(matrix(ncol = number.cols,nrow = number.rows,dimnames = list(row.names(data),colnames(data))))
    for (r in 1:number.rows) {
        for (c in 1:number.cols) {
            norm.fpkm[r,c] <- (as.numeric(data[r,c]) - fpkm.stats[r,2]) / fpkm.stats[r,3]
            print(r,c)
        }
    }
    result <- list(fpkm.stats,norm.fpkm)
    return(norm.fpkm)
}


fpkm.robust <- function(data) {
    number.rows <- nrow(data)
    number.cols <- ncol(data)
    fpkm.stats <- data.frame(
        Genes = row.names(data),
        Median.by.row = rep(NA,number.rows),
        IQR.by.row = rep(NA,number.rows)
    );
    for (i in 1:number.rows) {
        fpkm.stats[i,2]<- median(as.numeric(c(t(data[i,]))),na.rm = TRUE)
        fpkm.stats[i,3]<- IQR(as.numeric(c(t(data[i,]))),na.rm = TRUE) * 0.7413
    }
    median.norm.fpkm <- data.frame(matrix(ncol = number.cols,nrow = number.rows,dimnames = list(row.names(data),colnames(data))))
    for (r in 1:number.rows) {
        for (c in 1:number.cols) {
            median.norm.fpkm[r,c] <- (as.numeric(data[r,c]) - fpkm.stats[r,2]) / fpkm.stats[r,3]
        }
    }
    result <- list(fpkm.stats,median.norm.fpkm)
    return(median.norm.fpkm)
}

fpkm.nonnorm <- fpkm.tumor.symbol.filter
### Fixing Column names in fpkm data to match clinical id
colnames(fpkm.nonnorm) <- substr(
    x = colnames(fpkm.nonnorm),
    start = 1,
    stop = nchar(colnames(fpkm.nonnorm))-4
);
#replacing '.' with '-'
colnames(fpkm.nonnorm) <- gsub(
    pattern = '\\.',
    replacement = '-',
    x = colnames(fpkm.nonnorm)
);

sub.fpkm <- fpkm.nonnorm[1:5,]


age.color <- replace(
    as.numeric(c(t(brca.clinic['Diagnosis.Age']))),
    (50 >= brca.clinic['Diagnosis.Age']),
    'white'
);
age.color <- replace(
    age.color,
    (59 >= brca.clinic['Diagnosis.Age'] & 50 <= brca.clinic['Diagnosis.Age']),
    'lightpink1'
);
age.color <- replace(
    age.color,
    (69 >= brca.clinic['Diagnosis.Age'] & 60 <= brca.clinic['Diagnosis.Age']),
    'palevioletred2'
);
age.color <- replace(
    age.color,
    (79 >= brca.clinic['Diagnosis.Age'] & 70 <= brca.clinic['Diagnosis.Age']),
    'hotpink2'
);
age.color <- replace(
    age.color,
    (80 <= brca.clinic['Diagnosis.Age']),
    'violetred'
);
age.colours <- as.character(force.colour.scheme(outlier.brca.clinic$Diagnosis.Age, scheme = "age.categorical.prostate"))
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

stage.color <- substring(text = c(t(brca.clinic['American.Joint.Committee.on.Cancer.Tumor.Stage.Code'])),first = 1,last = 2)
stage.color <- replace(
    stage.color,
    ('T1' == brca.clinic['American.Joint.Committee.on.Cancer.Tumor.Stage.Code']),
    'white'
);
stage.color <- replace(
    stage.color,
    ('T2' == brca.clinic['American.Joint.Committee.on.Cancer.Tumor.Stage.Code']),
    'darkseagreen'
);
stage.color <- replace(
    stage.color,
    ('T3' == brca.clinic['American.Joint.Committee.on.Cancer.Tumor.Stage.Code']),
    'darkolivegreen'
);
stage.color <- replace(
    stage.color,
    ('T4' == brca.clinic['American.Joint.Committee.on.Cancer.Tumor.Stage.Code']),
    'darkgreen'
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
create.heatmap(
    x = subset(brca.clinic,select = c('American.Joint.Committee.on.Cancer.Tumor.Stage.Code')),
    input.colours = TRUE,
    clustering.method = 'none',
    colour.scheme = force.colour.scheme(scheme = 'stage'),
    same.as.matrix = TRUE,
    force.grid.col = FALSE,
    col.colour = 'black',
    xaxis.tck = 0,
    xlab.label = ' ',
    xaxis.lab = rep('',nrow(brca.clinic)),
    xaxis.cex = 0,
    yaxis.tck = 0,
    ylab.label = ' ',
    yaxis.lab = rep('',nrow(brca.clinic)),
    legend.cex = 0,
    print.colour.key = FALSE,
);
#Heat map for fpkm normalized by robust z-scores
create.heatmap(
    x = median.normalized.fpkm,
    clustering.method = 'diana',
    filename = generate.filename(
        project.stem = 'Cancer',
        file.core = 'Heatmedian2',
        extension = 'tiff'
    ),
    resolution = 300,
    height = 48,
    width = 48
    );
## Heatmap for normalized fpkm by z-score
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
    width = 48
    );

### Bar plot with total number of outlier genes per patient
patient.outlier <- create.barplot(
    formula = Outlier.totals ~ Patient.ID,
    data = outlier.brca.clinic,
    #filename = generate.filename(
     #   project.stem = 'Cancer',
      #  file.core = 'TOPHIST2',
       # extension = 'tiff'
        #),
    xlab.label = '',
    xat = '',
    xaxis.tck = 0,
    ylab.label = '# Outliers',
    yaxis.tck = 0,
    yat = c(0,25,50),
    )
gene.outlier <- create.barplot(
    formula = Gene.ID ~ Outlier.totals,
    data = outlier.totals.by.gene,
    #filename = generate.filename(
     #   project.stem = 'Cancer',
      #  file.core = 'SIDEHIST2',
       # extension = 'tiff'
    #)
    height = 30,
    ylab.label = '',
    yaxis.lab = rep('',nrow(outlier.totals.by.gene)),
    yat = 0,
    xat = c(0,2,4),
    xlab.top.y = 1,
    xaxis.tck = 0,
    xlab.label = '',
    xlab.top.label = '# Outlier Patients',
    plot.horizontal = T
)

outlier.brca.clinic.ID <- outlier.brca.clinic[order(outlier.brca.clinic$Patient.ID,decreasing = T),]

create.multipanelplot(
    layout.height = 4,
    layout.width = 1,
    plot.objects.heights = C(1,1,1,1),
    plot.objects = list(age.map,subtype.map,PFS.map,Tumor.stage.map)
    )
create.multipanelplot(
    plot.objects = list(age.map,subtype.map,PFS.map,Tumor.stage.map)
)
create.multipanelplot(
    filename = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Multiplot',
        extension = 'tiff'
        ),
    layout.height = 6,
    layout.width = 2,
    plot.objects.heights = C(2,10,1,1,1,1),
    plot.objects.widths = c(10,2),
    plot.objects = list(zscore.map)
    );
create.multipanelplot(
    filename = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Multiplot3',
        extension = 'tiff'
    ),
    layout.width = 2,
    layout.height = 6,
    plot.objects = list(patient.outlier,zscore.map,gene.outlier,age.map,Tumor.stage.map,PFS.map,subtype.map),
    plot.objects.heights = c(2,10,1,1,1,1),
    layout.skip = c(FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE)
);
create.multipanelplot(
    plot.objects = list(patient.outlier,gene.outlier)
    );
##### Working mutlipanelplot for covariates
covariates <- create.multipanelplot(
    plot.objects = list(age.map,Tumor.stage.map,PFS.map,subtype.map),y.spacing = -2.3
)
