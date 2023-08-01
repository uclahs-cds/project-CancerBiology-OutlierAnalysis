### TCGA_LandscapeHeatmap #########################################################################

### Preamble ######################################################################################
library(stats);
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(BoutrosLab.statistics.survival);

#mean.sub.TCGA.zscore.fpkm <- readRDS(
#"C:/Users/jre83/OneDrive/UCLA B.I.G. Summer/Paul_Lab/Survival_Analysis_TCGA-BRCA/TCGAZscore.RDS")

### Analysis ######################################################################################
# - order the patients
distance.z.score.abundance.patient <- dist(t(median.sub.TCGA.robust.zscore.fpkm), method = "euclidean")
cluster.distance.z.score.abundance.patient <- hclust(distance.z.score.abundance.patient, method = "ward.D2")
col_order <- cluster.distance.z.score.abundance.patient$order
z.score.abundance.patient.order <- median.sub.TCGA.robust.zscore.fpkm[, col_order]
outlier.by.patient <- outlier.patient.tag.01.t.p.order[,col_order]
#totaling number of outlier genes per patient
outlier.totals.by.patient <- data.frame(
    Patient.ID = seq(1,ncol(outlier.by.patient)),
    Outlier.totals = apply(
        X = outlier.by.patient,
        MARGIN = 2,
        FUN = sum
        )
    );

# - order the genes
distance.z.score.abundance.gene <- dist(median.sub.TCGA.robust.zscore.fpkm, method = "euclidean")
cluster.distance.z.score.abundance.gene <- hclust(distance.z.score.abundance.gene, method = "ward.D2")
row_order <- cluster.distance.z.score.abundance.gene$order
z.score.abundance.patient.gene.order <- median.sub.TCGA.robust.zscore.fpkm[row_order, ]
outlier.by.gene <- outlier.patient.tag.01.t.p.order[row_order,]
#prep for bar plot of number of outliers per gene
outlier.totals.by.gene <- data.frame(
    Gene.ID = seq(1,nrow(outlier.by.gene)),
    Outlier.totals = apply(
        X = outlier.by.gene,
        MARGIN = 1,
        FUN = sum
        )
    );

# Order clinical data
brca.sample.clinic <- gsub("-", ".", outlier.brca.clinic$Patient.ID, fixed = TRUE);
brca.clinic.sort <- brca.clinic[order(brca.clinic$Subtype),];
brca.clinic.order <- brca.clinic.sort[substr(colnames(outlier.by.patient), 1, 12),]


# other clinical status
#   1. subtype
#   2. Progression free survival: Progress.Free.Survival..Months.
#   3. tumor stage: American.Joint.Committee.on.Cancer.Tumor.Stage.Code
#   4. sgage

# 1
brca.clinic.order.data <- data.frame(brca.clinic.order$Subtype);
brca.clinic.order.data[is.na(brca.clinic.order.data$brca.clinic.order.Subtype),] <- 6;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_Basal',] <- 1;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_Her2',] <- 2;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_LumA',] <- 3;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_LumB',] <- 4;
brca.clinic.order.data[brca.clinic.order.data$brca.clinic.order.Subtype == 'BRCA_Normal',] <- 5;
brca.clinic.order.data.num <- data.frame(as.numeric(brca.clinic.order.data$brca.clinic.order.Subtype));
rownames(brca.clinic.order.data.num) <- colnames(outlier.patient.tag.01.t.p.order);

# 2
# brca.clinic.order.data.pfs <- data.frame(brca.clinic.order$Progress.Free.Survival..Months.);
# brca.clinic.order.data.pfs[is.na(brca.clinic.order.data.pfs$brca.clinic.order.Progress.Free.Survival..Months.),] <- 6;
# brca.clinic.order.data.pfs[brca.clinic.order.data.pfs$brca.clinic.order.Progress.Free.Survival..Months. <= 12,] <- 10;
# brca.clinic.order.data.pfs[brca.clinic.order.data.pfs$brca.clinic.order.Progress.Free.Survival..Months.> 12 &
#                                brca.clinic.order.data.pfs$brca.clinic.order.Progress.Free.Survival..Months. <= 24,] <- 9;
# brca.clinic.order.data.pfs[brca.clinic.order.data.pfs$brca.clinic.order.Progress.Free.Survival..Months.> 24 &
#                                brca.clinic.order.data.pfs$brca.clinic.order.Progress.Free.Survival..Months. <= 60,] <- 8;
# brca.clinic.order.data.pfs[brca.clinic.order.data.pfs$brca.clinic.order.Progress.Free.Survival..Months. > 60,] <- 7;
# brca.clinic.order.data.pfs.num <- data.frame(as.numeric(brca.clinic.order.data.pfs$brca.clinic.order.Progress.Free.Survival..Months.));
# rownames(brca.clinic.order.data.pfs.num) <- colnames(outlier.patient.tag.01.t.p.order);

# 3 
# T2 / T3 / T4
brca.clinic.order.data.stage <- brca.clinic.order$American.Joint.Committee.on.Cancer.Tumor.Stage.Code;
brca.clinic.order.data.stage.sub <- substr(brca.clinic.order.data.stage, 1, 2)
grade.sub.stage <- data.frame();
for (i in 1:length(brca.clinic.order.data.stage.sub)) {
    if (is.na(brca.clinic.order.data.stage.sub[i])) {
        grade.sub.stage[i,1] <- 6;
    }
    else if (brca.clinic.order.data.stage.sub[i] == 'T1') {
        grade.sub.stage[i,1] <- 11;
    }
    else if (brca.clinic.order.data.stage.sub[i] == 'T2') {
        grade.sub.stage[i,1] <- 12;
    }
    else if (brca.clinic.order.data.stage.sub[i] == 'T3') {
        grade.sub.stage[i,1] <- 13;
    }
    else if (brca.clinic.order.data.stage.sub[i] == 'T4') {
        grade.sub.stage[i,1] <- 14;
    }
}
grade.sub.stage <- data.frame(grade.sub.stage);
rownames(grade.sub.stage) <- colnames(outlier.patient.tag.01.t.p.order);

# 4
# 40~49 / 50~59 / 60~69 / 70<
brca.clinic.order.data.age <- brca.clinic.order$Diagnosis.Age;
age.stage <- data.frame();
for (i in 1:length(brca.clinic.order.data.age)) {
    if (is.na(brca.clinic.order.data.age[i])) {
        age.stage[i,1] <- 6;
    }
    else if (brca.clinic.order.data.age[i] < 50) {
        age.stage[i,1] <- 15;
    }
    else if (brca.clinic.order.data.age[i] >= 50 & brca.clinic.order.data.age[i] < 60) {
        age.stage[i,1] <- 16;
    }
    else if (brca.clinic.order.data.age[i] >= 60 & brca.clinic.order.data.age[i] < 70) {
        age.stage[i,1] <- 17;
    }
    else if (brca.clinic.order.data.age[i] >= 70 & brca.clinic.order.data.age[i] < 80) {
        age.stage[i,1] <- 18;
    }
    else if (brca.clinic.order.data.age[i] >= 80) {
        age.stage[i,1] <- 19;
    }
}
rownames(age.stage) <- colnames(outlier.patient.tag.01.t.p.order);

# Color code
chr.name <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
              '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
              '21', '22', 'MT', 'X', 'Y')
bio <- c("DNA", "RNA", "Protein", "Carbohydrate", "Lipid");
chr.colours <- BoutrosLab.plotting.general:::force.colour.scheme(bio, scheme = 'biomolecule');
sub.col <- c(chr.colours, 'grey90');
cnv.col <- c('dodgerblue4', 'white', 'darkred');

age.ColourFunction <- colorRamp(c('white', 'deeppink3'), space = 'Lab');
age.col <- rgb(age.ColourFunction(seq(0, 1, length.out = 5)), maxColorValue = 255);
stage.ColourFunction <- colorRamp(c('white', 'darkgreen'), space = 'Lab');
stage.col <- rgb(stage.ColourFunction(seq(0, 1, length.out = 4)), maxColorValue = 255);
#rt.color <- c('lightsalmon3', 'darkolivegreen3');
all.col <- c(sub.col, stage.col, age.col);

# Covariate heatmap
all.clinical.order <- cbind(
    brca.clinic.order.data.num,
    grade.sub.stage,
    age.stage
    );
subtype.heat <-  BoutrosLab.plotting.general:::create.heatmap(
    x = all.clinical.order,
    clustering.method = 'none',
    colour.scheme = all.col, 
    total.colours = 16,
    row.colour = 'black',
    col.colour = 'black',
    grid.row = TRUE, 
    grid.col = FALSE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.lab = 'Subtype',
    ylab.label = '',
    yaxis.cex = 0,
    print.colour.key = FALSE,
    # force.grid.row = TRUE,
    # force.grid.col = TRUE
    );  

#top bar plot with total number of outliers per patient
patient.outlier <- create.barplot(
    formula = Outlier.totals ~ Patient.ID,
    data = outlier.totals.by.patient,
    #filename = generate.filename(
    #   project.stem = 'Cancer',
    #  file.core = 'TOPHIST2',
    # extension = 'tiff'
    #),
    xlab.label = '',
    xat = '',
    xaxis.tck = 0,
    ylab.label = '# Outlier',
    ylab.cex = 0.75,
    yaxis.tck = 0,
    yat = c(0,25,50),
    yaxis.cex = 1,
    ylimits = c(0,55)
    );

# right barplot with number of outlier patients per gene
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
    xaxis.cex = 1,
    xlab.top.y = 1,
    xaxis.tck = 0,
    xlab.label = '',
    xlab.top.label = '# Outlier Patients',
    xlab.top.cex = 0.75,
    plot.horizontal = T
    );


#Make legend
legend.clinic <- BoutrosLab.plotting.general:::legend.grob(
    list(
        legend = list(
            colours = cnv.col,
            title = expression(underline('Robust Z-Score')), 
            labels = c('-5','5'),
            size = 2,
            label.cex = 1, 
            continuous = TRUE,
            height = 3
        ),
        # create legend for subtype
        legend = list(
            colours = sub.col[1:5],
            title = expression(underline('Subtype')), 
            labels = c('Basal', 'Her2', 'LumA', 'LumB', 'Normal'),
            size = 2,
            title.cex = 0.5,
            label.cex = 0.5, 
            border = 'black'
        ),
        # create legend for stage
        legend = list(
            colours = stage.col,
            title = expression(underline('Tumour stage')), 
            labels = c('T1', 'T2', 'T3', 'T4'),
            size = 2,
            title.cex = 0.5,
            label.cex = 0.5, 
            border = 'black'
        ),
        # create legend for age
        legend = list(
            colours = age.col,
            title = expression(underline('Age')), 
            labels = c(expression('\u2264 '*'50'), '50 - 59', '60 - 69','70 - 79', expression('\u2265 '*'80')),
            size = 2,
            title.cex = 0.5,
            label.cex = 0.5, 
            border = 'black'
        )
    ),
    title.just = 'left',
    title.fontface = 'plain',
    size = 3);
## Warning in FUN(X[[i]], ...): 'x' is NULL so the result will be NULL

# Main heatmap
#   - show the outlier patient status
patient.heat <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(median.sub.TCGA.robust.zscore.fpkm),
    clustering.method = 'ward.D2',
    #clustering.method = 'none',
    cluster.dimensions = 'both',
    plot.dendrograms = 'none',
    force.clustering = TRUE, 
    yaxis.lab = NULL,
    yaxis.cex = 0,
    xaxis.lab = NULL,
    xaxis.cex = 0,
    #main = paste('Outlier status of patients'),
    main.cex = 1.7,
    grid.col = FALSE,
    print.colour.key = FALSE,
    ylab.label = 'Outlier',
    ylab.cex = 1.3,
    #xlab.label = 'Patients',
    xlab.cex = 1.3,
    yaxis.tck = 0,
    xaxis.tck = 0,
    grid.lwd = 0,
    axes.lwd = 0,
    colour.scheme = cnv.col,
    # colour.scheme = c('white', 'red3'),
    colour.centering.value = 0,
    at = seq(-5,5),
    colourkey.cex = 1,
    # covariate.legend = legend.col.p,
    legend.side = 'right',
    legend.title.just = 'left',
    rows.distance.method = 'euclidean',
    cols.distance.method = 'euclidean'
    );  
## CAUTION: grid.lwd is DEPRECATED!  Use row.lwd/col.lwd. Using: 0

#old heatmap
#BoutrosLab.plotting.general:::
# create.multiplot(
#     filename = 'WOW10.tiff',
#     plot.objects = list(subtype.heat,patient.heat,gene.outlier,patient.outlier),
#     x.relation = "free",
#     y.relation = "free",
#     main = 'Outlier status',
#     xlab.label = 'Patients',
#     ylab.label = 'Outliers',
#     plot.layout = c(2, 3),
#     layout.skip = c(F,T,F,F,F,T),
#     panel.heights = c(1, 0.15),
#     panel.widths = c(1,.1),
#     ylab.padding = 2,
#     xlab.to.xaxis.padding = 1.1,
#     y.spacing = -0.7,
#     x.spacing = 0.3,
#     main.cex = 1.7,
#     xaxis.cex = 0,
#     xaxis.lab = NULL,
#     yaxis.lab = list(
#         NULL,NULL,NULL,NULL,NULL,
#         seq(0,50,25)
#         ),
#     yaxis.cex = 0,
#     yaxis.tck = 0,
#     ylab.cex = 1.5,
#     xlab.cex = 1.5,
#     xaxis.rot = 90,
#     xaxis.tck = 0,
#     legend = list(right = list(fun = legend.clinic)),
#     print.new.legend = TRUE,
#     resolution = 500,
#     height = 10,
#     width = 15
#     );

# Final Heatmap
create.multipanelplot(
    filename = generate.filename(
        project.stem = 'CancerBiology.OutlierAnalysis',
        file.core = 'TCGA LandscapeMap',
        extension = 'tiff'
        ),
    xlab.label = 'Patient',
    main = 'Outlier Status',
    main.cex = 2,
    layout.width = 2,
    layout.height = 3,
    height = 10,
    width = 15,
    resolution = 300,
    plot.objects = list(patient.outlier,patient.heat,gene.outlier,subtype.heat),
    plot.objects.heights = c(2.5,11,2),
    plot.objects.widths = c(10,1),
    layout.skip = c(FALSE,TRUE,FALSE,FALSE,FALSE,TRUE),
    x.spacing = -1.75,
    y.spacing = c(-10,-6,-1),
    legend = list(
        right = list(
            fun =legend.clinic #for heat maps
            )
        ),
    right.padding = 1,
    right.legend.padding = 1
    );

save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA Landscape',
    extension = 'txt')
    );

