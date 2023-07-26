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
outlier.meta.patient <- modified.clinic(
    data = meta.clinic.patient,
    outlier.data = outlier.patient.tag.01,
    TCGA = FALSE,
    META = TRUE,
    sample = FALSE)

outlier.meta.sample <- modified.clinic(
    data = meta.clinic.sample,
    outlier.data = outlier.patient.tag.01,
    TCGA = FALSE,
    META = TRUE,
    sample = TRUE)

### Lymph Nodes ###
lymph.nodes <- kruskal.continous(
    continous.var = 'Lymph.nodes.examined.positive',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Lymph-Nodes-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Lymph.nodes.examined.positive',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient,
    kruskal = lymph.nodes,
    ylimits = c(-1,40),
    text.x = 3.8,
    text.y = 38
    );

### Nottingham ###
nottingham <- kruskal.continous(
    continous.var = 'Nottingham.prognostic.index',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Nottingham-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Nottingham.prognostic.index',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient,
    kruskal = nottingham,
    ylimits = c(-1,10),
    text.x = 3.8,
    text.y = 9
    );

### Age ###
age <- kruskal.continous(
    continous.var = 'Age.at.Diagnosis',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'AgeDiagnosis-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Age.at.Diagnosis',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient,
    kruskal = age,
    ylimits = c(-1,110),
    text.x = 3.8,
    text.y = 105
    );

### Survival ###
survival <- kruskal.continous(
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Survival-Correlation2',
        extension = 'tiff'
        ),
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient,
    kruskal = survival,
    ylimits = c(-20,400),
    text.x = 3.8,
    text.y = 390
    );

### Relapse.free ###
relapse <- kruskal.continous(
    continous.var = 'Relapse.Free.Status..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'RelapseFREE-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Relapse.Free.Status..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.patient,
    kruskal = relapse,
    ylimits = c(-20,400),
    text.x = 3.8,
    text.y = 390
    );

### Tumor.Size ###
size <- kruskal.continous(
    continous.var = 'Tumor.Size',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.sample
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Tumor-size-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Tumor.Size',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.sample,
    kruskal = size,
    ylimits = c(-20,200),
    text.x = 3.8,
    text.y = 190
    );

### TMB..nonsynonymous. ###
TMB <- kruskal.continous(
    continous.var = 'TMB..nonsynonymous.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.sample
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'TMB-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'TMB..nonsynonymous.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.meta.sample,
    kruskal = TMB,
    ylimits = c(-5,50),
    text.x = 3.8,
    text.y = 48
    );
