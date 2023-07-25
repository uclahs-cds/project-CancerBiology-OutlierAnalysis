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

### Formatting Oulier Data ########################################################################
#Description:
    #formatting outlier data for mergeing with clinical data for analysis
#Input Variables
    #data:
    #data frame with outlier data per patient (row) and by gene (column)
    #by.gene:
    #binary - sums outlier totals per gene for heatmap 
    #by.patient:
    #binary - sums outlier totals per patient for clinical significance 
#Output Variables
    # results:
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


save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Correlation-Function',
    extension = 'txt')
    );