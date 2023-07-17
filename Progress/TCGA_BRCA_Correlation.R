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
### continuous.cor funtion ########################################################################
#
#
####
#
#
#
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
        }

    results <- list(result)
    return(results)
    }

boxing <- function(continous.var,subgroups,data) {
    
}

create.boxplot(
    formula = Subgroups ~ Aneuploidy.Score,
    data = outlier.brca.clinic,
    xlimits = c(-5,40)
)

create.boxplot(
    formula = Subgroups ~ Aneuploidy.Score,
    data = outlier.brca.clinic,
    xlimits = c(-5,40)
)

create.boxplot(
    formula = Subgroups ~ Buffa.Hypoxia.Score,
    data = outlier.brca.clinic,
    xlimits = c(-60,60)
)

create.boxplot(
    formula = Subgroups ~ Aneuploidy.Score,
    data = outlier.brca.clinic
)
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
)
################# NOT NEEDED #########################################################################################

#reordering data for outlier totals from most to least
#outlier.totals.descend <- outlier.totals[order(outlier.totals$Outlier.totals,decreasing = TRUE),];
#Outlier group patient id
outlier.group <- outlier.totals[outlier.totals$Outlier.totals >= 1,];
outlier.ID <- outlier.group$Patient.ID;
#Not outlier group patinet id
not.outlier.group <- outlier.totals[outlier.totals$Outlier.totals == 0,];
not.outlier.ID <- not.outlier.group$Patient.ID;
#combining outlier and non outlier info per patient
groups <- data.frame( 
    group = c(rep('No outlier',length(not.outlier.ID)),rep('Outlier',length(outlier.ID))),
    Patient.ID = c(not.outlier.ID,outlier.ID)
    );
#merging patient info with extra column for if they had at least one outlier
outlier.classification.brca.clinic <- merge(
    x = brca.clinic,
    y = groups,
    by = 'Patient.ID',
);
outlier.brca.clinic <- merge(
    x = outlier.classification.brca.clinic,
    y = outlier.totals.subgroups,
    by = 'Patient.ID'
)
#seperate dataframes
outlier <- data.frame(outlier.brca.clinic[outlier.brca.clinic["group"] == 'Outlier',])
no.outlier <- data.frame(outlier.brca.clinic[outlier.brca.clinic["group"] == 'No outlier',])

create.boxplot(
    formula = as.numeric(Diagnosis.Age) ~ as.numeric(Buffa.Hypoxia.Score),
    data = outlier.brca.clinic
)

create.dotmap()



continous <- c('Aneuploidy.Score','Buffa.Hypoxia.Score','Last.Communication.Contact.from.intial.Pathologic.Diagnostic.Date','Birth.from.Intitial.Pathological.Diagnosis.Date','Disease.Free.Months','Months.of.disease.specefic.survival','Fraction.Genome.Altered','MSIsensor.Score','Mutation.Count','Overall.Survival..Months','Progress.Free.Survival..Months','Ragnum.Hypoxia.Score','')


####Trouble shooting groups of outliers
zeros <- replace(x = outlier.totals$Outlier.totals,list = (0 != outlier.totals$Outlier.totals),values = NA)
zeros <- outlier.totals[0 == outlier.totals$Outlier.totals,]

ones <- replace(x = outlier.totals$Outlier.totals,list = (1 != outlier.totals$Outlier.totals),values = NA)
ones <- outlier.totals[1 == outlier.totals$Outlier.totals,]

twos <- replace(x = outlier.totals$Outlier.totals,list = (2 != outlier.totals$Outlier.totals),values = NA)
twos <- outlier.totals[2 == outlier.totals$Outlier.totals,]

three.plus <- replace(x = outlier.totals$Outlier.totals,list = (3 >= outlier.totals$Outlier.totals),values = NA)
three.plus <- outlier.totals[3 <= outlier.totals$Outlier.totals,]

number.outliers <- data.frame( 
    group = c(rep('No outlier',length(not.outlier.ID)),rep('Outlier',length(outlier.ID))),
    Patient.ID = c(not.outlier.ID,outlier.ID)
);


outlier.subgroups <- outlier.totals$Outlier.totals
outlier.subgroups <- replace(
    x = outlier.subgroups,
    list = (3 <= outlier.subgroups),
    values = 3)
outlier.totals.subgroups <- data.frame(outlier.totals, subgroups = outlier.subgroups)
