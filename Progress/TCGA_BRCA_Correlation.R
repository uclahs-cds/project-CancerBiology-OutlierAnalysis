### Exploring TCGA-BRCA ###########################################################################

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

### kruskal.continous funtion #####################################################################
#Description:
    #calculating the correlation between number of outliers and continuous variable
#Input Variables
    #continous.var:
        # NAME OF column that contains the continuous variable of interest
    #subgroups:
        #name of column that contains the subgroup classifications
    #data:
        #the name of the data frame that contains your clinical data and groups '0','1','2','3+' 
        #for outliers merge clinical data with a data frame with ID and outlier groups 
    #posthoc
        #for use with only one continuous variable to display the differences between groups
#Output Variables
    # results:
        #contains the Kruskal wallis test result and pair wise Wilcox test

kruskal.continous <- function(continous.var,subgroups,data,posthoc = FALSE) {
    number.variables = length(continous.var)
    number.subgroups = length(subgroups)
    number.datasets = length(data)
    x.variable.names <- gsub( 
        pattern = '\\.',
        replacement = ' ',
        x = continous.var
        );
    result <- array(
        dim = c(number.variables,1,2),
        dimnames = c(list(continous.var,subgroups,c('p.value','effect.size'))))
    for (i in 1:number.variables) {
        correlation <- kruskal.test(
            x = as.vector(data[,continous.var[i]]),
            g = as.vector(c(t(data[subgroups]))),
            );
        result[i,1,1] <- correlation$p.value;
        result[i,1,2] <- correlation$statistic;
        if (TRUE == posthoc) {
            single.comparison <- pairwise.wilcox.test(
                x = as.numeric(c(t(data[,continous.var[i]]))),
                g = as.vector(c(t(data[subgroups]))),
                p.adjust.method = 'bonf'
                );
            results <- list(
                Kruskal.Wall.Test = result,
                Paired.wise.Wilcox = single.comparison
                );
            }
    }
    if (FALSE == posthoc) {
        results <- list(Kruskal.Wall.Test = result);
        }
    return(results)
    }

### outlier.create.boxplot function ###############################################################
#Description
    #produces graphical representations for continuous variables and number of outlier groups
#Input Variables
    #continuous.var
        #name of column with continuous variable
    #subgroups
        #name of column with subgroup classification
    #data
        #dataframe with clinical data and groups '0','1','2','3+' for outliers
        #merge clinical data with a data frame with ID and outlier groups 
    #kruskal
        #variable name that contains the results from kruskal.continous function without posthoc
    #file.name
        #desired file name
#Output Variables
    #plot
        #box plot 
outlier.create.boxplot <- function(continous.var,subgroups,data,kruskal,ylimits=c(0,100),file.name,text.x,text.y) {
    #easy variable name fix
    x.variable.name <- gsub( 
        pattern = '\\.',
        replacement = ' ',
        x = continous.var
    );
    #creating a subset of data to remove NA's
    subdata <- data.frame(target.var = as.numeric(c(t(data[continous.var]))),data[subgroups])
    subdata <- na.omit(subdata)
    text.pval <- paste('p-value:',format(kruskal$Kruskal.Wall.Test[1],digits = 4))
    #plotting
    plot <- create.boxplot(
        filename = file.name, #fix when plots are good
        formula = as.formula(
            paste0(
                paste(
                    'target.var', 
                    ' ~ ',
                    sep = ''
                    ),
                paste(
                    subgroups,
                    sep = ''
                    )
                )
            ),
        data = subdata,
        ylimits = ylimits,
        ylab.label = x.variable.name,
        add.stripplot = TRUE,
        #strip.col = 'black',
        points.col = 'slategrey',
        xlab.label = 'Number of Outliers per Patient',
        resolution = 300,
        add.text = TRUE,
        text.labels = text.pval,
        text.x = text.x,
        text.y = text.y,
        );
    return(plot);
    }

### Data Analysis #################################################################################
#totaling number of outlier genes per patient
outlier.totals <- data.frame(
    Patient.ID = colnames(outlier.patient.tag.01.t.p.order),
    Outlier.totals = apply(
        X = outlier.patient.tag.01.t.p.order,
        MARGIN = 2,
        FUN = sum
        )
    );
##fixing the names of patients to match across the outlier data to the clinical data
#removing trailing characters and '.'
outlier.totals$Patient.ID <- substr(
    x = outlier.totals$Patient.ID,
    start = 1,
    stop = nchar(outlier.totals$Patient.ID)-4
    );
#replacing '.' with '-'
outlier.totals$Patient.ID <- gsub(
    pattern = '\\.',
    replacement = '-',
    x = outlier.totals$Patient.ID
    );
#creating the group identifiers ###NOT NECESSARY###
outlier.classification <- outlier.totals$Outlier.totals
outlier.classification <- replace(
    x = outlier.classification,
    list = (0 != outlier.classification),
    values = 'Outlier'
    );
outlier.classification <- replace(
    x = outlier.classification,
    list = (0 == outlier.classification),
    values = 'No.Outlier'
    );

#creating group identifiers '0','1','2',and '3+' ###REQUIRED###
outlier.subgroups <- outlier.totals$Outlier.totals;
outlier.subgroups <- replace(
    x = outlier.subgroups,
    list = (3 <= outlier.subgroups),
    values = '3+'
    );

outlier.totals.subgroups <- data.frame(
    Class = outlier.classification,
    outlier.totals,
    Subgroups = outlier.subgroups
    );
brca.clinic.2 <- brca.clinic
outlier.brca.clinic <- merge(
    x = brca.clinic.2,
    y = outlier.totals.subgroups,
    by = 'Patient.ID'
    );
#a concatenation of all possible continuous variables
continous <- c('Aneuploidy.Score','Buffa.Hypoxia.Score','Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date','Birth.from.Initial.Pathologic.Diagnosis.Date','Disease.Free..Months.','Months.of.disease.specific.survival','Fraction.Genome.Altered','MSIsensor.Score','Mutation.Count','Overall.Survival..Months.','Progress.Free.Survival..Months.','Ragnum.Hypoxia.Score');

###### PLOTTING ###################################################################################

### Aneuploidy Score
aneuploidy <- kruskal.continous(
    continous.var = 'Aneuploidy.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Aneoploidy-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Aneuploidy.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    kruskal = aneuploidy,
    ylimits = c(-1,40),
    text.x = 3.8,
    text.y = 38
    );

### Buffa Hypoxia Score
buffa <- kruskal.continous(
    continous.var = 'Buffa.Hypoxia.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Buffa.Hypoxia-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Buffa.Hypoxia.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-50,60),
    kruskal = buffa,
    text.x = 3.8,
    text.y = 57
    );

### Last Communication 
last.comm <- kruskal.continous(
    continous.var = 'Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Last.Communication-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-100,5000),
    kruskal = last.comm,
    text.x = 3.8,
    text.y = 4800
    );

### 'Birth.from.Intitial.Pathological.Diagnosis.Date'
birth <- kruskal.continous(
    continous.var = 'Birth.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Birth.to.Diagnosis-Correlation',
        extension = 'tiff'
        ),
    continous.var ='Birth.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-1,4000),
    kruskal = birth,
    text.x = 3.8,
    text.y = 20
    );

### Disease.Free..Months.
Disease.free <- kruskal.continous(
    continous.var = 'Disease.Free..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Disease.Free.Months-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Disease.Free..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200),
    kruskal = Disease.free,
    text.x = 3.8,
    text.y = 190
    );

### Months of disease 
Disease <- kruskal.continous(
    continous.var = 'Months.of.disease.specific.survival',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Months.with.Disease-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Months.of.disease.specific.survival',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200),
    kruskal = Disease,
    text.x = 3.8,
    text.y = 190
    );

### Fractrion of genome altered
altered.genome <- kruskal.continous(
    continous.var = 'Fraction.Genome.Altered',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Genome.Altered-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Fraction.Genome.Altered',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-0.1,1.3),
    kruskal = altered.genome,
    text.x = 3.8,
    text.y = 1.2
    );

### MSIsensor.Score
sensor <- kruskal.continous(
    continous.var = 'MSIsensor.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'MSIsensor-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'MSIsensor.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-0.9,6),
    kruskal = sensor,
    text.x = 3.8,
    text.y = 5.6
    );

### MSI.MANTIS.Score
mantis <- kruskal.continous(
    continous.var = 'MSI.MANTIS.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
);

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'MSI.MANTIS-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'MSI.MANTIS.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(0,1),
    kruskal = mantis,
    text.x = 3.8,
    text.y = 0.95
    );

### 'Mutation.Count'
mutation <- kruskal.continous(
    continous.var = 'Mutation.Count',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Mutation.Count-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Mutation.Count',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,260),
    kruskal = mutation,
    text.x = 3.8,
    text.y = 250
    );

###'Overall.Survival..Months.'
survival <- kruskal.continous(
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Overall.Survival-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200),
    kruskal = survival,
    text.x = 3.8,
    text.y = 190
    );

### 'Progress.Free.Survival..Months.'
progress.free <- kruskal.continous(
    continous.var = 'Progress.Free.Survival..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Progress.Free.Survival-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Progress.Free.Survival..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200),
    kruskal = progress.free,
    text.x = 3.8,
    text.y = 190
    );

### 'Ragnum.Hypoxia.Score'
ragnum <- kruskal.continous(
    continous.var = 'Ragnum.Hypoxia.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    file.name = generate.filename(
        project.stem = 'CancerBiology-OutlierAnalysis',
        file.core = 'Ragnum.Hypoxia-Correlation',
        extension = 'tiff'
        ),
    continous.var = 'Ragnum.Hypoxia.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-30,40),
    kruskal = ragnum,
    text.x = 3.8,
    text.y = 37
    );

save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Correlation',
    extension = 'txt')
    );