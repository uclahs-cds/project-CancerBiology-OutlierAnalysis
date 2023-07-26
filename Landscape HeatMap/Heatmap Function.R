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

### create.landscape.heatmap ######################################################################
#Description
    #
#Input Variables
    # (1) age
    #column name in data that contains ages
    # (2) stage
    #column name in data that contains stages
    # (3) progress.free
    #column name in data that contains progress free survival in months
    # (4) subtype
    #column name in data that contains subtypes
    # (5) patient.ID
    #column name in data that contains patients ID number
    # (6) data
    #data frame of clinical data
    # (7) normalized.fpkm
    #data frame of normalized fpkm values produced by fpkm.normalization or fpkm.robust.zscore
    # (8) patient.total
    #column name in outlier.totals.patient that contains number of outliers
    # (9) outlier.totals.patient
    #data frame with Patient.ID and outlier totals
    # (10) gene.ID
    #column in outlier.totals.per.gene with gene id
    # (11) gene.totals
    #column in outlier.totals.per.gene with number of outliers per gene
    # (12) outlier.totals.per.gene
    #data frame with gene id and outlier totals per gene
    # (13) file.name
    #name for plot that will be produced
#Output
    # (1) file output in wd
    # mutipanel heatmap produced in working directory
create.landscape.heatmap <- function(age,stage,progress.free,subtype,patient.ID,data,normalized.fpkm,patient.total,outlier.totals.patients,gene.ID,gene.totals,outlier.totals.per.gene,file.name) { 
    zscore.map <- create.heatmap(
        x = normalized.fpkm,
        clustering.method = 'none',
        force.clustering = TRUE,
        resolution = 300,
        height = 48,
        width = 48,
        print.colour.key = FALSE,
        
        );
    patient.outlier <- create.barplot(
        formula = as.formula(
            paste0(
                paste(
                    patient.total, 
                    ' ~ ',
                    sep = ''
                ),
                paste(
                    patient.ID,
                    sep = ''
                    )
                )
            ),
        data = outlier.totals.patients,
        xlab.label = '',
        xat = '',
        xaxis.tck = 0,
        ylab.label = '# Outliers',
        yaxis.tck = 0,
        yat = c(0,25,50),
        ylab.cex = 1,
        yaxis.cex = 1
        );
    
    gene.outlier <- create.barplot(
        formula = as.formula(
            paste0(
                paste(
                    gene.ID, 
                    ' ~ ',
                    sep = ''
                    ),
                paste(
                    gene.totals,
                    sep = ''
                    )
                )
            ),
        #formula = Gene.ID ~ Outlier.totals,
        data = outlier.totals.by.gene,
        #filename = generate.filename(
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
        xlab.top.cex = 1,
        xaxis.cex = 1
        );
    
    age.colours <- as.character(force.colour.scheme(c(t(data[age])), scheme = "age.categorical.prostate"))
    data.color <- cbind(data,age.colours)
    age.map <- create.heatmap(
        x = subset(data.color,select = c('age.colours')),
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
    
    subtype.color <- replace(
        as.numeric(c(t(data[subtype]))),
        ('BRCA_Basal' == data[subtype] | 'Basal' == data[subtype]),
        'powderblue'
    );
    subtype.color <- replace(
        subtype.color,
        ('BRCA_Her2' == data[subtype] | 'Her2' == data[subtype]),
        'lightskyblue'
    );
    subtype.color <- replace(
        subtype.color,
        ('BRCA_LumA' == data[subtype] | 'LumA' == data[subtype]),
        'mediumpurple1'
    );
    subtype.color <- replace(
        subtype.color,
        ('BRCA_LumB' == data[subtype] | 'LumB' == data[subtype]),
        'palevioletred2'
    );
    subtype.color <- replace(
        subtype.color,
        ('BRCA_Normal' == data[subtype] | 'Normal' == data[subtype]),
        'yellow3'
    );
    data.color <- cbind(data,subtype.color)
    subtype.map <- create.heatmap(
        x = subset(data.color,select = c('subtype.color')),
        input.colours = TRUE,
        clustering.method = 'none',
        same.as.matrix = TRUE,
        force.grid.col = FALSE,
        col.colour = 'black',
        xaxis.tck = 0,
        xlab.label = ' ',
        xaxis.lab = rep('',nrow(data)),
        xaxis.cex = 0,
        yaxis.tck = 0,
        ylab.label = ' ',
        yaxis.lab = rep('',nrow(data)),
        legend.cex = 0,
        print.colour.key = FALSE,
    );
    
    
    PFS.color <- replace(
        as.numeric(c(t(data[progress.free]))),
        (60 <= data[progress.free]),
        'white'
    );
    PFS.color <- replace(
        PFS.color,
        (60 > data[progress.free] & 24 < data[progress.free]),
        'darkorchid1'
    );
    PFS.color <- replace(
        PFS.color,
        (24 >= data[progress.free] & 12 < data[progress.free]),
        'mediumorchid4'
    );
    PFS.color <- replace(
        PFS.color,
        (12 >= data[progress.free]),
        'purple4'
    );
    
    data.color <- cbind(data,PFS.color)
    PFS.map <- create.heatmap(
        x = subset(data.color,select = c('PFS.color')),
        input.colours = TRUE,
        clustering.method = 'none',
        same.as.matrix = TRUE,
        force.grid.col = FALSE,
        col.colour = 'black',
        xaxis.tck = 0,
        xlab.label = ' ',
        xaxis.lab = rep('',nrow(data.color)),
        xaxis.cex = 0,
        yaxis.tck = 0,
        ylab.label = ' ',
        yaxis.lab = rep('',nrow(data.color)),
        legend.cex = 0,
        print.colour.key = FALSE,
    );
    #stage heatmap
    stage.colours <- as.character(force.colour.scheme(c(t(data[stage])), scheme = "stage"))
    data.color <- cbind(data,stage.colours)
    Tumor.stage.map <- create.heatmap(
        x = subset(data.color,select = c('stage.colours')),
        #x = stage.color,
        input.colours = TRUE,
        clustering.method = 'none',
        same.as.matrix = TRUE,
        force.grid.col = FALSE,
        col.colour = 'black',
        colour.scheme = force.colour.scheme(scheme = 'stage'),
        xaxis.tck = 0,
        xlab.label = ' ',
        xaxis.lab = rep('',nrow(data.color)),
        xaxis.cex = 0,
        yaxis.tck = 0,
        ylab.label = ' ',
        yaxis.lab = rep('',nrow(data.color)),
        legend.cex = 0,
        print.colour.key = FALSE,
        height = 4,
        width = 24
    );
    #stage.color.key <- data.frame(
     #   stage = unique(c(t(data[stage]))),
      #  colour = unique(stage.colours)
    #   )
    main.legend <- legend.grob(
        list(
            legend = list(
                title = expression(underline('Z-score')),
                colours = c('red3','white','blue'),
                labels = c(-5,5),
                size = 1.5,
                border = 'black',
                continuous = TRUE
                ),
            legend = list(
                title = expression(underline('Subtype')),
                colours = c('powderblue','lightskyblue','mediumpurple1','palevioletred2','yellow3'),
                labels = c('Basal', 'Her2','LumA','LumB','Normal'),
                size = 1.5,
                border = 'black'
                ),
            legend = list(
                title = expression(underline('Progression Free Survival')),
                colours = c('white','darkorchid1','mediumorchid4','purple4'),
                labels = c('> 60', '24-60','12-24','< 12'),
                size = 1.5,
                border = 'black'
                ),
            legend = list(
                title = expression(underline('Tumor Stage')),
                colours = c('green'),
                labels = c('BAd'),
                size = 1.5,
                border = 'black'
                ),
            legend = list(
                title = expression(underline('Age')),
                colours = c('gray100','gray75','gray50','gray25','gray0'),
                labels = c('<40','40-50', '50-65','65-70','>70'),
                size = 1.5,
                border = 'black'
                )
            )
        )
    
    create.multipanelplot(
        filename = file.name,
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
        y.spacing = c(-3,-2,-5,-5,-5,-5),
        legend = list(
            right = list(
                fun = main.legend #for heat maps
                )
            )
        );
    }

### DATA Analysis #################################################################################
#prep for stage 
outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code[grep('X',outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)]<- 'X'
outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code[grep('IV',outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)]<- 'IV'
outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code[grep('III',outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)]<- 'III'
outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code[grep('II',outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)]<- 'II'
outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code[grep('STAGE I',outlier.brca.clinic$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)]<- 'I'

#prep for bar plot of number of outliers per gene
outlier.totals.by.gene <- data.frame(
    Gene.ID = row.names(outlier.patient.tag.01.t.p.order),
    Outlier.totals = apply(
        X = outlier.patient.tag.01.t.p.order,
        MARGIN = 1,
        FUN = sum
        )
    );
outlier.genes <- data.frame(
    Gene.ID = row.names(outlier.patient.tag.01.t.p.order),
    outlier.patient.tag.01.t.p.order
)
Outlier.genes <- merge(
    y = outlier.genes,
    x = outlier.totals.by.gene,
    by = 'Gene.ID'
)

create.landscape.heatmap(
    age = 'Diagnosis.Age',
    stage = 'Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code',
    progress.free = 'Progress.Free.Survival..Months.',
    subtype = 'Subtype',
    patient.ID = 'Patient.ID',
    data = outlier.brca.clinic,
    normalized.fpkm = Z.score,
    patient.total = 'Outlier.totals',
    outlier.totals.patients = outlier.totals ,
    gene.ID = 'Gene.ID',
    gene.totals = 'Outlier.totals',
    outlier.totals.per.gene = outlier.totals.by.gene,
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Multiplot_function_3',
        extension = 'tiff'
        )
    );


save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Heatmap-Function',
    extension = 'txt')
    );