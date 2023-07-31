### TCGA_LandscapeHeatmap #########################################################################

### Preamble ######################################################################################
library(stats);
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(BoutrosLab.statistics.survival);

#mean.sub.META.zscore.fpkm <- readRDS(
#"C:/Users/jre83/OneDrive/UCLA B.I.G. Summer/Paul_Lab/Survival_Analysis_TCGA-BRCA/METAzscore.RDS")

### Analysis ######################################################################################

# - order the patients
meta.distance.z.score.abundance.patient <- dist(t(mean.sub.META.zscore.fpkm), method = "euclidean")
meta.cluster.distance.z.score.abundance.patient <- hclust(meta.distance.z.score.abundance.patient, method = "ward.D2")
meta_col_order <- meta.cluster.distance.z.score.abundance.patient$order
z.score.abundance.patient.order <- mean.sub.META.zscore.fpkm[, meta_col_order]
outlier.by.patient <- outlier.patient.tag.01[,meta_col_order]
#totaling number of outlier genes per patient
meta.outlier.totals.by.patient <- data.frame(
    Patient.ID = seq(1,ncol(outlier.by.patient)),
    Outlier.totals = apply(
        X = outlier.by.patient,
        MARGIN = 2,
        FUN = sum
    )
    );

# - order the genes
distance.z.score.abundance.gene <- dist(mean.sub.META.zscore.fpkm, method = "euclidean")
cluster.distance.z.score.abundance.gene <- hclust(distance.z.score.abundance.gene, method = "ward.D2")
row_order <- cluster.distance.z.score.abundance.gene$order
z.score.abundance.patient.gene.order <- mean.sub.TCGA.zscore.fpkm[row_order, ]
meta.outlier.by.gene <- outlier.patient.tag.01[row_order,]
#prep for bar plot of number of outliers per gene
meta.outlier.totals.by.gene <- data.frame(
    Gene.ID = seq(1,nrow(meta.outlier.by.gene)),
    Outlier.totals = apply(
        X = meta.outlier.by.gene,
        MARGIN = 1,
        FUN = sum
    )
    );



# Order clinical data
meta.sample.clinic <- gsub("-", ".", outlier.meta.patient$Patient.ID, fixed = TRUE);
rownames(outlier.meta.clinic) <- meta.sample.clinic;
meta.clinic.sort <- meta.clinic[order(meta.clinic$Subtype),];
meta.clinic.order <- meta.clinic.sort[substr(colnames(outlier.patient.tag.01), 1, 12),]


# other clinical status
#   1. subtype
#   2. Progression free survival: Progress.Free.Survival..Months.
#   3. tumor stage: American.Joint.Committee.on.Cancer.Tumor.Stage.Code
#   4. sgage

# 1
meta.clinic.order.data <- data.frame(meta.clinic.patient$Pam50...Claudin.low.subtype);
meta.clinic.order.data[is.na(meta.clinic.order.data$meta.clinic.patient.Pam50...Claudin.low.subtype),] <- 6;
meta.clinic.order.data[meta.clinic.order.data$meta.clinic.patient.Pam50...Claudin.low.subtype == 'Basal',] <- 1;
meta.clinic.order.data[meta.clinic.order.data$meta.clinic.patient.Pam50...Claudin.low.subtype == 'Her2',] <- 2;
meta.clinic.order.data[meta.clinic.order.data$meta.clinic.patient.Pam50...Claudin.low.subtype == 'LumA',] <- 3;
meta.clinic.order.data[meta.clinic.order.data$meta.clinic.patient.Pam50...Claudin.low.subtype == 'LumB',] <- 4;
meta.clinic.order.data[meta.clinic.order.data$meta.clinic.patient.Pam50...Claudin.low.subtype == 'Normal',] <- 5;
meta.clinic.order.data.num <- data.frame(as.numeric(meta.clinic.order.data$meta.clinic.patient.Pam50...Claudin.low.subtype));
rownames(meta.clinic.order.data.num) <- colnames(outlier.patient.tag.01);

# 2
# meta.clinic.order.data.pfs <- data.frame(meta.clinic.order$Progress.Free.Survival..Months.);
# meta.clinic.order.data.pfs[is.na(meta.clinic.order.data.pfs$meta.clinic.order.Progress.Free.Survival..Months.),] <- 6;
# meta.clinic.order.data.pfs[meta.clinic.order.data.pfs$meta.clinic.order.Progress.Free.Survival..Months. <= 12,] <- 10;
# meta.clinic.order.data.pfs[meta.clinic.order.data.pfs$meta.clinic.order.Progress.Free.Survival..Months.> 12 &
#                                meta.clinic.order.data.pfs$meta.clinic.order.Progress.Free.Survival..Months. <= 24,] <- 9;
# meta.clinic.order.data.pfs[meta.clinic.order.data.pfs$meta.clinic.order.Progress.Free.Survival..Months.> 24 &
#                                meta.clinic.order.data.pfs$meta.clinic.order.Progress.Free.Survival..Months. <= 60,] <- 8;
# meta.clinic.order.data.pfs[meta.clinic.order.data.pfs$meta.clinic.order.Progress.Free.Survival..Months. > 60,] <- 7;
# meta.clinic.order.data.pfs.num <- data.frame(as.numeric(meta.clinic.order.data.pfs$meta.clinic.order.Progress.Free.Survival..Months.));
# rownames(meta.clinic.order.data.pfs.num) <- colnames(outlier.patient.tag.01.t.p.order);
# 
# # 3 
# # T2 / T3 / T4
# meta.clinic.order.data.stage <- meta.clinic.order$American.Joint.Committee.on.Cancer.Tumor.Stage.Code;
# meta.clinic.order.data.stage.sub <- substr(meta.clinic.order.data.stage, 1, 2)
# grade.sub.stage <- data.frame();
# for (i in 1:length(meta.clinic.order.data.stage.sub)) {
#     if (is.na(meta.clinic.order.data.stage.sub[i])) {
#         grade.sub.stage[i,1] <- 6;
#     }
#     else if (meta.clinic.order.data.stage.sub[i] == 'T1') {
#         grade.sub.stage[i,1] <- 11;
#     }
#     else if (meta.clinic.order.data.stage.sub[i] == 'T2') {
#         grade.sub.stage[i,1] <- 12;
#     }
#     else if (meta.clinic.order.data.stage.sub[i] == 'T3') {
#         grade.sub.stage[i,1] <- 13;
#     }
#     else if (meta.clinic.order.data.stage.sub[i] == 'T4') {
#         grade.sub.stage[i,1] <- 14;
#     }
# }
# grade.sub.stage <- data.frame(grade.sub.stage);
# rownames(grade.sub.stage) <- colnames(outlier.patient.tag.01.t.p.order);

# 4
# 40~49 / 50~59 / 60~69 / 70<
meta.clinic.order.data.age <- meta.clinic.patient$Age.at.Diagnosis;
age.stage <- data.frame();
for (i in 1:length(meta.clinic.order.data.age)) {
    if (is.na(meta.clinic.order.data.age[i])) {
        age.stage[i,1] <- 6;
        }
    else if (meta.clinic.order.data.age[i] < 50) {
        age.stage[i,1] <- 15;
        }
    else if (meta.clinic.order.data.age[i] >= 50 & meta.clinic.order.data.age[i] < 60) {
        age.stage[i,1] <- 16;
        }
    else if (meta.clinic.order.data.age[i] >= 60 & meta.clinic.order.data.age[i] < 70) {
        age.stage[i,1] <- 17;
        }
    else if (meta.clinic.order.data.age[i] >= 70 & meta.clinic.order.data.age[i] < 80) {
        age.stage[i,1] <- 18;
        }
    else if (meta.clinic.order.data.age[i] >= 80) {
        age.stage[i,1] <- 19;
        }
    }
rownames(age.stage) <- colnames(outlier.patient.tag.01);



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
pnf.ColourFunction <- colorRamp(c('white', 'purple4'), space = 'Lab');
pnf.col <- rgb(pnf.ColourFunction(seq(0, 1, length.out = 4)), maxColorValue = 255);
stage.ColourFunction <- colorRamp(c('white', 'darkgreen'), space = 'Lab');
stage.col <- rgb(stage.ColourFunction(seq(0, 1, length.out = 4)), maxColorValue = 255);
# rt.color <- c('lightsalmon3', 'darkolivegreen3');
#all.col <- c(sub.col, pnf.col, stage.col, age.col);

all.col <- c(sub.col, age.col);



# Covariate heatmap
all.clinical.order <- cbind(
    meta.clinic.order.data.num,
    age.stage
    );

meta.subtype.heat <-  BoutrosLab.plotting.general:::create.heatmap(
    x = all.clinical.order,
    clustering.method = 'none',
    colour.scheme = all.col, 
    total.colours = 10,
    row.colour = 'black',
    col.colour = 'black',
    grid.row = TRUE, 
    grid.col = FALSE, 
    yaxis.tck = 0, 
    xaxis.tck = 0,
    yaxis.lab = '',
    ylab.label = '',
    yaxis.cex = 0,
    print.colour.key = FALSE,
    # force.grid.row = TRUE,
    # force.grid.col = TRUE
    );  

# top bar plot for number of outlier genes per gene
meta.patient.outlier <- create.barplot(
    formula = Outlier.totals ~ Patient.ID,
    data = meta.outlier.totals.by.patient,
    #filename = generate.filename(
    #   project.stem = 'Cancer',
    #  file.core = 'TOPHIST2',
    # extension = 'tiff'
    #),
    xlab.label = '',
    xat = '',
    xaxis.tck = 0,
    ylab.label = '# Outlier',
    ylab.cex = 1,
    yaxis.tck = 0,
    yaxis.cex = 1,
    yat = c(0,50,100),
    ylimits = c(0,100)
    );

# right bar plot for number of outliers by gene
meta.gene.outlier <- create.barplot(
    formula = Gene.ID ~ Outlier.totals,
    data = meta.outlier.totals.by.gene,
    #filename = generate.filename(
    #   project.stem = 'Cancer',
    #  file.core = 'SIDEHIST2',
    # extension = 'tiff'
    #)
    height = 30,
    ylab.label = '',
    yaxis.lab = rep('',nrow(outlier.totals.by.gene)),
    yat = 0,
    xat = c(0,5,10,15),
    xaxis.cex = 1,
    xlimits = c(0,18),
    xlab.top.y = 1,
    xaxis.tck = 0,
    xlab.label = '',
    xlab.top.label = '# Outlier Patients',
    xlab.top.cex = 0.75,
    plot.horizontal = T
    );


#Make legend
meta.legend.clinic <- BoutrosLab.plotting.general:::legend.grob(
    list(
        legend = list(
            colours = cnv.col,
            title = expression(underline('Outlier status')), 
            labels = c('-5','5'),
            size = 2,
            label.cex = 1, 
            continuous = TRUE,
            height = 3
            ),
        legend = list(
            colours = 'darkred', 
            labels = c(expression('\u2265 '*'5')),
            size = 2,
            label.cex = 1, 
            continuous = TRUE,
            height = 1
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
        # create legend for pfs
        #legend = list(
         #   colours = pnf.col,
          #  title = expression(underline('Progression free survival')), 
           # labels = c(expression('\u2265 '*'60'), '24 - 60', '12 - 24', expression('\u2264 '*'12')),
            #size = 2,
            #title.cex = 0.5,
            #label.cex = 0.5, 
            #border = 'black'
            #),
        # create legend for stage
        #legend = list(
         #   colours = stage.col,
          #  title = expression(underline('Tumour stage')), 
           # labels = c('T1', 'T2', 'T3', 'T4'),
            #size = 2,
            #title.cex = 0.5,
            #label.cex = 0.5, 
            #border = 'black'
            #),
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
    size = 3
    );

## Warning in FUN(X[[i]], ...): 'x' is NULL so the result will be NULL
# Main heatmap
#   - show the outlier patient status
meta.patient.heat <- BoutrosLab.plotting.general:::create.heatmap(
    x = t(mean.sub.META.zscore.fpkm),
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
#     filename = 'WOW11.tiff',
#     plot.objects = list(subtype.heat,patient.heat,meta.gene.outlier,meta.patient.outlier),
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
#     xaxis.lab = ,
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

#Final Heatmap
create.multipanelplot(
    filename = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'META Landscape',
        extension = 'tiff'
        ),
    xlab.label = 'Patient',
    #ylab.label = 'Outliers',
    main = 'Outlier Status',
    main.cex = 2,
    layout.width = 2,
    layout.height = 3,
    height = 10,
    width = 15,
    resolution = 300,
    plot.objects = list(meta.patient.outlier,meta.patient.heat,meta.gene.outlier,meta.subtype.heat),
    plot.objects.heights = c(2.5,11,1.5),
    plot.objects.widths = c(10,1),
    layout.skip = c(FALSE,TRUE,FALSE,FALSE,FALSE,TRUE),
    x.spacing = -1.75,
    y.spacing = c(-10,-6,-1),
    legend = list(
        right = list(
            fun =meta.legend.clinic #for heat maps
            )
        ),
    right.padding = 1,
    right.legend.padding = 1
    );

save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'META Landscape',
    extension = 'txt')
    );
