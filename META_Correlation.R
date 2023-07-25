### META_Correlation ###########################################################################

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

### Data Analysis #################################################################################
meta.data.patient <- meta.clinic.patient[5:2513,1:23]
meta.data.sample <- meta.clinic.sample[5:2513,1:12]

#adding or adjusting ID for merging and analysing 
meta.data.patient <- data.frame(Patient.ID = row.names(meta.data.patient),meta.data.patient)
colnames(meta.data.sample)[1] <- 'Patient.ID'

#totaling number of outlier genes per patient
outlier.totals <- data.frame(
    Patient.ID = colnames(outlier.patient.tag.01),
    Outlier.totals = apply(
        X = outlier.patient.tag.01,
        MARGIN = 2,
        FUN = sum
        )
    );

#creating group identifiers '0','1','2',and '3+' ##REQUIRED###
outlier.totals$Outlier.totals <- replace(
    x = outlier.totals$Outlier.totals,
    list = (3 <= outlier.totals$Outlier.totals),
    values = '3+'
    );
#fixing names for merge
outlier.totals$Patient.ID <- gsub(
    pattern = '\\.',
    replacement = '-',
    x = outlier.totals$Patient.ID
    );

meta.patient.2 <- meta.data.patient
outlier.meta.patient <- merge(
    x = meta.patient.2,
    y = outlier.totals,
    by = 'Patient.ID'
    );

meta.sample.2 <- meta.data.sample
outlier.meta.sample <- merge(
    x = meta.sample.2,
    y = outlier.totals,
    by = 'Patient.ID'
    );

### Lymph Nodes ###
lymph.nodes <- kruskal.continous(
    continous.var = 'Lymph.nodes.examined.positive',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Lymph-Nodes-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Lymph.nodes.examined.positive',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient,
    kruskal = lymph.nodes,
    ylimits = c(-1,40),
    text.x = 3.8,
    text.y = 38
    );

### Nottingham ###
nottingham <- kruskal.continous(
    continous.var = 'Nottingham.prognostic.index',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Nottingham-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Nottingham.prognostic.index',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient,
    kruskal = nottingham,
    ylimits = c(-1,10),
    text.x = 3.8,
    text.y = 9
    );

### Age ###
age <- kruskal.continous(
    continous.var = 'Age.at.Diagnosis',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'AgeDiagnosis-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Age.at.Diagnosis',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient,
    kruskal = age,
    ylimits = c(-1,110),
    text.x = 3.8,
    text.y = 105
    );

### Survival ###
survival <- kruskal.continous(
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Survival-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient,
    kruskal = survival,
    ylimits = c(-20,400),
    text.x = 3.8,
    text.y = 390
    );

### Relapse.free ###
relapse <- kruskal.continous(
    continous.var = 'Relapse.Free.Status..Months.',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'RelapseFREE-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Relapse.Free.Status..Months.',
    subgroups = 'Outlier.totals',
    data = outlier.meta.patient,
    kruskal = relapse,
    ylimits = c(-20,400),
    text.x = 3.8,
    text.y = 390
    );

### Tumor.Size ###
size <- kruskal.continous(
    continous.var = 'Tumor.Size',
    subgroups = 'Outlier.totals',
    data = outlier.meta.sample
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Tumor-size-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Tumor.Size',
    subgroups = 'Outlier.totals',
    data = outlier.meta.sample,
    kruskal = size,
    ylimits = c(-20,200),
    text.x = 3.8,
    text.y = 190
    );

### TMB..nonsynonymous. ###
TMB <- kruskal.continous(
    continous.var = 'TMB..nonsynonymous.',
    subgroups = 'Outlier.totals',
    data = outlier.meta.sample
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'TMB-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'TMB..nonsynonymous.',
    subgroups = 'Outlier.totals',
    data = outlier.meta.sample,
    kruskal = TMB,
    ylimits = c(-5,50),
    text.x = 3.8,
    text.y = 48
    );
