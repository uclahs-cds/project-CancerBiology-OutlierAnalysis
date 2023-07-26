### Continuous_Variables_Clinical_Correlation_Functions ###########################################

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
a<-read.csv(file = '41467_2015_BFncomms9971_MOESM1236_ESM.csv')

### modified.clinic function ######################################################################
##Description
    # adding additional columns to clinical data for analysis
##Input Variables
    # (1) data
    # data frame with clinincal variables
    # (2) outlier.data
    # data frame with genes on the row and patients on column, 1s indicating an outlier gene
    # (3) TCGA
    # binary - is data from TCGA
    # (4) META
    # binary - is data from META
    # (5) sample
    # binary - is data from META sample
##Output Variables
    # (1) outlier.clinic
    # a data frame with the original data containing three to four more columns that contain
    # outlier classification, total outliers per patient, and outlier subgroups 0,1,2,3+

modified.clinic <- function(data,outlier.data,TCGA=TRUE,META=FALSE,sample=FALSE){
    #totaling number of outlier genes per patient
    outlier.totals <- data.frame(
        Patient.ID = colnames(outlier.data),
        Outlier.totals = apply(
            X = outlier.data,
            MARGIN = 2,
            FUN = sum
        )
    );
    #creating the group identifiers ###NOT NECESSARY###
    #outlier.classification <- outlier.totals$Outlier.totals
    outlier.classification <- replace(
        x = outlier.totals$Outlier.totals,
        list = (0 != outlier.totals$Outlier.totals),
        values = 'Outlier'
    );
    outlier.classification <- replace(
        x = outlier.classification,
        list = (0 == outlier.classification),
        values = 'No.Outlier'
    );
    #creating group identifiers '0','1','2',and '3+' ###REQUIRED###
    #outlier.subgroups <- outlier.totals$Outlier.totals;
    outlier.subgroups <- replace(
        x = outlier.totals$Outlier.totals,
        list = (3 <= outlier.totals$Outlier.totals),
        values = '3+'
    );
    
    outlier.subgroups.collection <- data.frame(
        Outlier.Class = outlier.classification,
        outlier.totals,
        Outlier.Subgroups = outlier.subgroups
    );
    #merge
    if (TRUE == TCGA & FALSE == META) {
        outlier.subgroups.collection$Patient.ID<- substr(
            x = outlier.subgroups.collection$Patient.ID,
            start = 1,
            stop = nchar(outlier.subgroups.collection$Patient.ID)-4
        );
        #replacing '.' with '-'
        outlier.subgroups.collection$Patient.ID <- gsub(
            pattern = '\\.',
            replacement = '-',
            x = outlier.subgroups.collection$Patient.ID
        );
        
        brca.clinic.2 <- brca.clinic
        outlier.clinic <- merge(
            x = brca.clinic.2,
            y = outlier.subgroups.collection,
            by = 'Patient.ID'
        );
    } else if (FALSE == TCGA & TRUE == META){
        #fixing names for merge
        outlier.subgroups.collection$Patient.ID <- gsub(
            pattern = '\\.',
            replacement = '-',
            x = outlier.subgroups.collection$Patient.ID
        );
        if (FALSE == sample) {
            #meta.patient
            data <- data[5:2513,1:23]
            data <- data.frame(Patient.ID = row.names(data),data)
            meta.patient.2 <- data
            outlier.clinic <- merge(
                x = meta.patient.2,
                y = outlier.subgroups.collection,
                by = 'Patient.ID'
            );
        } else {
            #meta.sample
            data <- data[5:2513,1:12]
            colnames(data)[1] <- 'Patient.ID'
            meta.sample.2 <- data
            outlier.clinic <- merge(
                x = meta.sample.2,
                y = outlier.subgroups.collection,
                by = 'Patient.ID'
            );
        }
    }
    #modified.clinic <- binary.survival(
    #   brca.data = outlier.brca.clinic,
    #  str.survival.col.name = 'Overall.Survival.Status',
    # binary.survival.col.name = 'Overall.Survival.Binary'
    #);
    return(outlier.clinic)
}

### Formatting Oulier Data ########################################################################
#Description:
    #formatting outlier data for mergeing with clinical data for analysis
#Input Variables
    # (1) data:
    #data frame with outlier data per patient (row) and by gene (column)
    # (2) by.gene:
    #binary - sums outlier totals per gene for heatmap 
    # (3) by.patient:
    #binary - sums outlier totals per patient for clinical significance 
#Output Variables
    # (1) results:
    # a data frame that contains the ID and outlier totals by patient or gene can be used to merge 
    #with clincal data 

create.outlier.totals <- function(data,by.gene = F,by.patient = T){
    if (FALSE == by.gene & TRUE == by.patient) {
        outlier.totals <- data.frame(
            Patient.ID = colnames(data),
            Outlier.totals = apply(
                X = data,
                MARGIN = 2,
                FUN = sum
                )
            );
        }
    if (TRUE == by.gene & FALSE == by.patient) {
        outlier.totals <- data.frame(
            Gene.ID = rownames(data),
            Outlier.totals = apply(
                X = data,
                MARGIN = 1,
                FUN = sum
                )
            );
        }
    result <- outlier.totals
    return(result)
    }

### kruskal.continous funtion #####################################################################
#Description:
    #calculating the correlation between number of outliers and continuous variable
#Input Variables
    # (1) continous.var:
    # NAME OF column that contains the continuous variable of interest
    # (2) subgroups:
    #name of column that contains the subgroup classifications
    # (3) data:
    #the name of the data frame that contains your clinical data and groups '0','1','2','3+' 
    #for outliers merge clinical data with a data frame with ID and outlier groups 
    # (4) posthoc
    #for use with only one continuous variable to display the differences between groups
#Output Variables
    # (1) results:
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
    # (1) continuous.var
    #name of column with continuous variable
    # (2) subgroups
    #name of column with subgroup classification
    # (3) data
    #dataframe with clinical data and groups '0','1','2','3+' for outliers
    #merge clinical data with a data frame with ID and outlier groups 
    # (4) kruskal
    #variable name that contains the results from kruskal.continous function without posthoc
    # (5) file.name
    #desired file name
#Output Variables
    # (1) plot
    #box plot 
outlier.create.boxplot <- function(continous.var,subgroups,data,kruskal,ylimits=c(0,100),file.name,text.x,text.y,y.lab.label = NULL) {
    #easy variable name fix
    if (is.null(y.lab.label)) {
        y.lab.label <- gsub( 
            pattern = '\\.',
            replacement = ' ',
            x = continous.var
            );
    } else {
        y.lab.label = y.lab.label
    }
    #creating a subset of data to remove NA's
    subdata <- data.frame(target.var = as.numeric(c(t(data[continous.var]))),data[subgroups]);
    subdata <- na.omit(subdata);
    text.pval <- paste('p-value:',format(kruskal$Kruskal.Wall.Test[1],digits = 4));
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
        ylab.label = y.lab.label,
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


save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Correlation-Function',
    extension = 'txt')
    );