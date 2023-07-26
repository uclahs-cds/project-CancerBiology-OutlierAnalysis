### TCGA_BRCA_Correlation #########################################################################

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
#####GET THE FUNCTIONS FORM Continous_Variable_Clinical_Association_Function.R

### Data Analysis #################################################################################
outlier.brca.clinic <- modified.clinic(
    data = brca.clinic,
    outlier.data = outlier.patient.tag.01.t.p.order,
    TCGA = TRUE,
    META = FALSE,
    sample = FALSE);

### Adjusting Purity Data
Purity.scores <- read.csv('41467_2015_BFncomms9971_MOESM1236_ESM.csv')
colnames(Purity.scores) <- Purity.scores[3,1:ncol(Purity.scores)]
Purity.scores <- Purity.scores[4:nrow(Purity.scores),]
Purity.scores$`Sample ID` <- substr(
    x = Purity.scores$`Sample ID`,
    start = 1,
    stop = nchar(Purity.scores$`Sample ID`)-1
    )
colnames(Purity.scores)[1] <- 'Sample.ID'
### merging with clinical data
outlier.brca.clinic <- merge(
    x = outlier.brca.clinic,
    y = Purity.scores,
    by = 'Sample.ID'
    )
#a concatenation of all possible continuous variables
#continous <- c('Aneuploidy.Score','Buffa.Hypoxia.Score','Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date','Birth.from.Initial.Pathologic.Diagnosis.Date','Disease.Free..Months.','Months.of.disease.specific.survival','Fraction.Genome.Altered','MSIsensor.Score','Mutation.Count','Overall.Survival..Months.','Progress.Free.Survival..Months.','Ragnum.Hypoxia.Score');

###### PLOTTING ###################################################################################

### Aneuploidy Score
aneuploidy <- kruskal.continous(
    continous.var = 'Aneuploidy.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Aneoploidy-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Aneuploidy.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    kruskal = aneuploidy,
    ylimits = c(-1,40),
    text.x = 3.8,
    text.y = 38
    );

### Buffa Hypoxia Score
buffa <- kruskal.continous(
    continous.var = 'Buffa.Hypoxia.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Buffa.Hypoxia-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Buffa.Hypoxia.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-50,60),
    kruskal = buffa,
    text.x = 3.8,
    text.y = 57
    );

### Last Communication 
last.comm <- kruskal.continous(
    continous.var = 'Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Last.Communication-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-100,5000),
    kruskal = last.comm,
    text.x = 3.8,
    text.y = 4800,
    ylab.cex = 1
    );

### 'Birth.from.Intitial.Pathological.Diagnosis.Date'
birth <- kruskal.continous(
    continous.var = 'Birth.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Birth.to.Diagnosis-Correlation',
        extension = 'tiff'
        ),
    continous.var ='Birth.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-1,4000),
    kruskal = birth,
    text.x = 3.8,
    text.y = 20
    );

### Disease.Free..Months.
Disease.free <- kruskal.continous(
    continous.var = 'Disease.Free..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Disease.Free.Months-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Disease.Free..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200),
    kruskal = Disease.free,
    text.x = 3.8,
    text.y = 190
    );

### Months of disease 
Disease <- kruskal.continous(
    continous.var = 'Months.of.disease.specific.survival',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Months.with.Disease-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Months.of.disease.specific.survival',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200),
    kruskal = Disease,
    text.x = 3.8,
    text.y = 190
    );

### Fractrion of genome altered
altered.genome <- kruskal.continous(
    continous.var = 'Fraction.Genome.Altered',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Genome.Altered-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Fraction.Genome.Altered',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-0.1,1.3),
    kruskal = altered.genome,
    text.x = 3.8,
    text.y = 1.2
    );

### MSIsensor.Score
sensor <- kruskal.continous(
    continous.var = 'MSIsensor.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'MSIsensor-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'MSIsensor.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-0.9,6),
    kruskal = sensor,
    text.x = 3.8,
    text.y = 5.6
    );

### MSI.MANTIS.Score
mantis <- kruskal.continous(
    continous.var = 'MSI.MANTIS.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'MSI.MANTIS-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'MSI.MANTIS.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(0,1),
    kruskal = mantis,
    text.x = 3.8,
    text.y = 0.95
    );

### 'Mutation.Count'
mutation <- kruskal.continous(
    continous.var = 'Mutation.Count',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Mutation.Count-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Mutation.Count',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,260),
    kruskal = mutation,
    text.x = 3.8,
    text.y = 250
    );

###'Overall.Survival..Months.'
survival <- kruskal.continous(
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Overall.Survival-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200),
    kruskal = survival,
    text.x = 3.8,
    text.y = 190
    );

### 'Progress.Free.Survival..Months.'
progress.free <- kruskal.continous(
    continous.var = 'Progress.Free.Survival..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Progress.Free.Survival-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Progress.Free.Survival..Months.',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200),
    kruskal = progress.free,
    text.x = 3.8,
    text.y = 190
    );

### 'Ragnum.Hypoxia.Score'
ragnum <- kruskal.continous(
    continous.var = 'Ragnum.Hypoxia.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Ragnum.Hypoxia-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Ragnum.Hypoxia.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-30,40),
    kruskal = ragnum,
    text.x = 3.8,
    text.y = 37
    );

### 'Winter.Hypoxia.Score'
winter <- kruskal.continous(
    continous.var = 'Winter.Hypoxia.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Winter.Hypoxia-Correlation',
        extension = 'tiff'
    ),
    continous.var = 'Winter.Hypoxia.Score',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-90,90),
    kruskal = winter,
    text.x = 3.8,
    text.y = 86
    );
####              ####
#### TUMOR PURITY ####
####              ####

### 'Estimate'
estimate <- kruskal.continous(
    continous.var = 'ESTIMATE',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Estimate-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'ESTIMATE',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(0,1.5),
    kruskal = estimate,
    text.x = 3.8,
    text.y = 1.45
    );

### 'Absolute'
absolute <- kruskal.continous(
    continous.var = 'ABSOLUTE',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Absolute-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'ABSOLUTE',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(0,1.5),
    kruskal = absolute,
    text.x = 3.8,
    text.y = 1.45
    );

### 'Lump'
lump <- kruskal.continous(
    continous.var = 'LUMP',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Lump-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'LUMP',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(0,1.5),
    kruskal = lump,
    text.x = 3.8,
    text.y = 1.45
    );

### 'IHC'
ihc <- kruskal.continous(
    continous.var = 'IHC',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'IHC-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'IHC',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(0,1.5),
    kruskal = ihc,
    text.x = 3.8,
    text.y = 1.45
    );

### 'CPE'
cpe <- kruskal.continous(
    continous.var = 'CPE',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'CPE-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'CPE',
    subgroups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(0,1.5),
    kruskal = cpe,
    text.x = 3.8,
    text.y = 1.45,
    y.lab.label = 'CPE Purity'
    );

save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Correlation',
    extension = 'txt')
    );