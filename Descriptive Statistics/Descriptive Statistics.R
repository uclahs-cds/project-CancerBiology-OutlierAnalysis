### Descriptive Statistics ########################################################################

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
### median.CI function ############################################################################
#Description
# calculatign 95% confidence intervals for medians of survival time
#Input Variables
    # (1) data
    # data 
#Output Variables
    # (1) result
    # matrix with upper and lower limit
median.CI <- function(data) {
    data <- na.omit(as.numeric(data))
    n = length(data);
    data2 <- sort(as.numeric(data));
    median <- 0.5;
    zscore <- 1.96;
    lower <- n * median - zscore * sqrt((n * median * (1 - median)));
    upper <- n * median + zscore * sqrt((n * median * (1 - median)));
    lower.limit <- data2[ceiling(lower)];
    upper.limit <- data2[ceiling(upper)];
    result <- data.frame(lower = lower.limit,upper = upper.limit);
    return(result);
    }

### subgroup.event function ########################################################################
#Description
# calculatign 95% confidence intervals for medians of survival time
#Input Variables
    # (1) data
    # data frame with survival analytics
    # (2) subgroup
    # name of the column with the subtypes
    # (3) subgroup.interest
    # which subgroup is of interest
    # (4) META
    # binary - is the data METABRIC
#Output Variables
    # (1) result
    # list with all subtype vs all versions of survival
subgroup.event <- function(data,subgroup,subgroup.interest,META=FALSE) {
    if (FALSE == META){
        overall <- data.frame(table(outlier.brca.clinic$Overall.Survival.Status[outlier.brca.clinic[subgroup] == subgroup.interest]));
        over.med <- median(
            as.numeric(
                outlier.brca.clinic$Overall.Survival..Months.[outlier.brca.clinic[subgroup] == subgroup.interest]),
            na.rm = T
            );
        overall.CI <- median.CI(outlier.brca.clinic$Overall.Survival..Months.[outlier.brca.clinic[subgroup] == subgroup.interest]);
        
        PFS <- data.frame(table(outlier.brca.clinic$Progression.Free.Status[outlier.brca.clinic[subgroup] == subgroup.interest]));
        PFS.med <- median(
            as.numeric(
                outlier.brca.clinic$Progress.Free.Survival..Months.[outlier.brca.clinic[subgroup] == subgroup.interest]),
            na.rm = T
            );
        PFS.CI <- median.CI(outlier.brca.clinic$Progress.Free.Survival..Months.[outlier.brca.clinic[subgroup] == subgroup.interest]);
        
        Disease.Free <- data.frame(table(outlier.brca.clinic$Disease.Free.Status[outlier.brca.clinic[subgroup] == subgroup.interest]));
        dmed <- median(
            as.numeric(
                outlier.brca.clinic$Disease.Free..Months.[outlier.brca.clinic[subgroup] == subgroup.interest]),
            na.rm = T
            );
        Disease.Free.CI <- median.CI(outlier.brca.clinic$Disease.Free..Months.[outlier.brca.clinic[subgroup] == subgroup.interest]);
        
        #Disease <- data.frame(table(outlier.brca.clinic$Disease.specific.Survival.status[outlier.brca.clinic[subgroup] == subgroup.interest]));
        #Disease.CI <- median.CI(outlier.brca.clinic$Months.of.disease.specific.survival[outlier.brca.clinic[subgroup] == subgroup.interest]);
        result <- list(
            Overall = overall,
            over.med,
            Overall.CI = overall.CI,
            PFS = PFS,
            PFS.med,
            PFS.CI = PFS.CI,
            Disease.Free = Disease.Free,
            dmed,
            DF.CI = Disease.Free.CI
            #Disease.Specific = Disease,
            #DS.CI = Disease.CI
            );
        return(result);
        }
    else if(TRUE == META){
        overall <- data.frame(table(outlier.meta.patient$Overall.Survival.Status[outlier.meta.patient[subgroup] == subgroup.interest]));
        over.med <- median(
            as.numeric(
                outlier.meta.patient$Overall.Survival..Months.[outlier.meta.patient[subgroup] == subgroup.interest]),
            na.rm = T
            );
        overall.CI <- median.CI(outlier.meta.patient$Overall.Survival..Months.[outlier.meta.patient[subgroup] == subgroup.interest]);
        
        recurrence <- data.frame(table(outlier.meta.patient$Relapse.Free.Status[outlier.meta.patient[subgroup] == subgroup.interest]));
        re.med <- median(
            as.numeric(
                outlier.meta.patient$Relapse.Free.Status..Months.[outlier.meta.patient[subgroup] == subgroup.interest]),
            na.rm = T
            );
        recurrence.CI <- median.CI(outlier.meta.patient$Relapse.Free.Status..Months.[outlier.meta.patient[subgroup] == subgroup.interest]);
        
        result <- list(
            Overall = overall,
            over.med,
            Overall.CI = overall.CI,
            Relapse = recurrence,
            re.med,
            Relapse.CI = recurrence.CI
            );
        return(result);
        }
    }
### Data Analysis #################################################################################
table(outlier.brca.clinic$Overall.Survival.Status);
median(brca.clinic$Overall.Survival..Months.,na.rm = T);
median.CI(brca.clinic$Overall.Survival..Months.);

table(outlier.brca.clinic$Progression.Free.Status);
median(brca.clinic$Progress.Free.Survival..Months.,na.rm = T);
median.CI(brca.clinic$Progress.Free.Survival..Months.);

table(outlier.brca.clinic$Disease.Free.Status);
median(brca.clinic$Disease.Free..Months.,na.rm = T);
median.CI(brca.clinic$Disease.Free..Months.);

table(outlier.brca.clinic$Disease.specific.Survival.status);
median.CI(outlier.brca.clinic$Months.of.disease.specific.survival);

subgroup.event(outlier.brca.clinic,'Subtype','BRCA_Basal');
subgroup.event(outlier.brca.clinic,'Subtype','BRCA_Her2');
subgroup.event(outlier.brca.clinic,'Subtype','BRCA_LumA');
subgroup.event(outlier.brca.clinic,'Subtype','BRCA_LumB');
subgroup.event(outlier.brca.clinic,'Subtype','BRCA_Normal');

subgroup.event(outlier.brca.clinic,'Outlier.Subgroups','0');
subgroup.event(outlier.brca.clinic,'Outlier.Subgroups','1');
subgroup.event(outlier.brca.clinic,'Outlier.Subgroups','2');
subgroup.event(outlier.brca.clinic,'Outlier.Subgroups','3+');


table(outlier.meta.patient$Overall.Survival.Status);
median(as.numeric(outlier.meta.patient$Overall.Survival..Months.),na.rm = T)
median.CI(outlier.meta.patient$Overall.Survival..Months.);

table(outlier.meta.patient$Relapse.Free.Status);
median(as.numeric(outlier.meta.patient$Relapse.Free.Status..Months.,na.rm = T));
median.CI(outlier.meta.patient$Relapse.Free.Status..Months.);

subgroup.event(outlier.meta.patient,"Pam50...Claudin.low.subtype",'Basal',META = T);
subgroup.event(outlier.meta.patient,"Pam50...Claudin.low.subtype",'Her2',META = T);
subgroup.event(outlier.meta.patient,"Pam50...Claudin.low.subtype",'LumA',META = T);
subgroup.event(outlier.meta.patient,"Pam50...Claudin.low.subtype",'LumB',META = T);
subgroup.event(outlier.meta.patient,"Pam50...Claudin.low.subtype",'Normal',META = T);

subgroup.event(outlier.meta.patient,'Outlier.Subgroups','0',META = T);
subgroup.event(outlier.meta.patient,'Outlier.Subgroups','1',META = T);
subgroup.event(outlier.meta.patient,'Outlier.Subgroups','2',META = T);
subgroup.event(outlier.meta.patient,'Outlier.Subgroups','3+',META = T);

save.session.profile(filename = generate.filename(
    project.stem = 'CancerBiology-OutlierAnalysis',
    file.core = 'TCGA-BRCA-Descriptive-Stats',
    extension = 'txt')
    );