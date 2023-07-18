### Exploring TCGA-BRCA ###########################################################################

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
### kruskal.continous funtion #####################################################################
#Description:
    #calculating the correlation between number of outliers and continuous variable
#Input Variables
    #continous.var:
        # NAME OF column that contains the continuous variable of interest
    #subgroups:
        #name of column that contains the subgroup classifications
    #data:
        #the name of the data frame that contains your clinical data
#Output Variables
    # results:
        #contains the Kruskal wallis test result and pair wise Wilcox test

#WORKING FOR ONE CONTINOUS VARABLE AT A TIME
kruskal.continous <- function(continous.var,subgroups,data) {
    number.variables = length(continous.var)
    number.subgroups = length(subgroups)
    number.datasets = length(data)
    x.variable.names <- gsub( 
        pattern = '\\.',
        replacement = ' ',
        x = continous.var
        );
    y.variable.names <- rev(unique(subgroups))
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
        single.comparison <- pairwise.wilcox.test(
            x = as.numeric(c(t(data[,continous.var[i]]))),
            g = as.vector(c(t(data[subgroups]))),
            p.adjust.method = 'bonf'
        )
        }
    results <- list(Kruskal.Wall.Test = result, Paired.wise.Wilcox = single.comparison)
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
        #dataframe with clinical data
    #file.name
        #desired file name
#Output Variables
    #plot
        #bar plot 
outlier.create.boxplot <- function(continous.var,subgroups,data,ylimits = c(0,100),file.name) {
    x.variable.name <- gsub( 
        pattern = '\\.',
        replacement = ' ',
        x = continous.var
    );
    continous.variable <- as.numeric(c(t(data[continous.var])))
    subdata <- data.frame('continous.var' = continous.variable,data[subgroups])
    subdata <- na.omit(subdata)
    plot <- create.boxplot(
        #filename = file.name,
        formula = as.formula(
            paste0(
                paste(
                    continous.var, 
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
        resolution = 300
        );
    return(plot);
    }
replace(x = outlier.brca.clinic,list = "NA",values = NA)

#Experimeting with boxplotting
create.boxplot(
    formula = Buffa.Hypoxia.Score ~ Subgroups,
    data = outlier.brca.clinic,
    ylimits = c(-60,60),
    #filename = generate.filename(
     #   project.stem = 'CancerBiology-OutlierAnalysis',
      #  file.core = 'Correlation-Boxplot',
       # extension = '.tiff'),
    resolution = 300,
    add.stripplot = TRUE,
    strip.col = 'black'
)




#colors <- c('rosybrown2','paleturquoise3','plum3','slateblue3')


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

outlier.subgroups <- outlier.totals$Outlier.totals;
outlier.subgroups <- replace(
    x = outlier.subgroups,
    list = (3 <= outlier.subgroups),
    values = '3+'
    );

outlier.totals.subgroups <- data.frame(Class = outlier.classification, outlier.totals, Subgroups = outlier.subgroups)
brca.clinic.2 <- brca.clinic
outlier.brca.clinic <- merge(
    x = brca.clinic.2,
    y = outlier.totals.subgroups,
    by = 'Patient.ID'
    );

continous <- c('Aneuploidy.Score','Buffa.Hypoxia.Score','Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date','Birth.from.Initial.Pathologic.Diagnosis.Date','Disease.Free..Months.','Months.of.disease.specific.survival','Fraction.Genome.Altered','MSIsensor.Score','Mutation.Count','Overall.Survival..Months.','Progress.Free.Survival..Months.','Ragnum.Hypoxia.Score');


### Aneuploidy Score BAD
kruskal.continous(
    continous.var = 'Aneuploidy.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Aneuploidy.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(0,50)
    );

### Buffa Hypoxia Score BAD
kruskal.continous(
    continous.var = 'Buffa.Hypoxia.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Buffa.Hypoxia.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-1,200)
    );

### Last Communication BAD
kruskal.continous(
    continous.var = 'Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-1,200)
    );

### 'Birth.from.Intitial.Pathological.Diagnosis.Date' BAD
kruskal.continous(
    continous.var = 'Birth.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var ='Birth.from.Initial.Pathologic.Diagnosis.Date',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-1,200)
    );

### Disease.Free..Months. BAD
kruskal.continous(
    continous.var = 'Disease.Free..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Disease.Free..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-1,200)
    );

### Months of disease GOOD NOT sig
kruskal.continous(
    continous.var = 'Months.of.disease.specific.survival',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Months.of.disease.specific.survival',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200)
    );

### Fractrion of genome altered BAD
kruskal.continous(
    continous.var = 'Fraction.Genome.Altered',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Fraction.Genome.Altered',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-1,200)
    );

### MSIsensor.Score GOOD SIG
kruskal.continous(
    continous.var = 'MSIsensor.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'MSIsensor.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-0.9,6)
    );

### 'Mutation.Count' GOOD SIG
kruskal.continous(
    continous.var = 'Mutation.Count',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Mutation.Count',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,250)
    );

###'Overall.Survival..Months.' GOOD not sig
kruskal.continous(
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Overall.Survival..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200)
    );

### 'Progress.Free.Survival..Months.' GOOD not SIG
kruskal.continous(
    continous.var = 'Progress.Free.Survival..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Progress.Free.Survival..Months.',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-5,200)
    );

### 'Ragnum.Hypoxia.Score' GOOD
kruskal.continous(
    continous.var = 'Ragnum.Hypoxia.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic
    );

outlier.create.boxplot(
    continous.var = 'Ragnum.Hypoxia.Score',
    subgroups = 'Subgroups',
    data = outlier.brca.clinic,
    ylimits = c(-30,40)
    );

#### TESTING
kruskal.continous(continous.var = 'Buffa.Hypoxia.Score',subgroups = 'Subgroups',data = outlier.brca.clinic);

kruskal.continous(continous.var = 'MSIsensor.Score',subgroups = 'Subgroups',data = outlier.brca.clinic);

outlier.create.boxplot(continous.var = 'MSIsensor.Score',subgroups = 'Subgroups',data = outlier.brca.clinic,ylimits = c(-0.9,6));



save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Correlation',
    extension = 'txt')
    );