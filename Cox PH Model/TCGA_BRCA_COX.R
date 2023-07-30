### TCGA-BRCA COX #########################################################################

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



cox.uni <- function(survival.time,survival.event,groups,data,reference.group) {
    
    survobj <- Surv(
        as.numeric(c(t(data[survival.time]))),
        as.numeric(survival.event)
        );
    group <- as.factor(c(t(data[groups])));
    group <- relevel(group,ref = reference.group)
    
    cox <- fit.coxmodel(
        groups = as.factor(group),
        survobj = survobj,
        return.cox.model = TRUE,
        );
    return(cox);
}

### Data Analysis #################################################################################
overall.outliers.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = outlier.progress.binary,
    groups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    reference.group = '0')

#LumA reference group 
overall.subtypes.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = outlier.progress.binary,
    groups = 'Subtype',
    data = outlier.brca.clinic,
    reference.group = 'BRCA_LumA'
    )


### Example
# will have to do outside of function
outlier.progress.binary <- replace(
    x = outlier.brca.clinic$Overall.Survival.Status,
    list = ('0:LIVING' == outlier.brca.clinic$Overall.Survival.Status),
    values = 0
    );
outlier.progress.binary <- replace(
    x = outlier.progress.binary,
    list = ('1:DECEASED' == outlier.progress.binary),
    values = 1
    );
#in functiom
survobj <- Surv(
    as.numeric(c(t(outlier.brca.clinic['Overall.Survival..Months.']))),
    as.numeric(outlier.progress.binary)
    );

fit.coxmodel(
    groups = as.factor(Outlier),
    survobj = survobj,
    return.cox.model = TRUE,
    );


#survobj <- Surv(
 #   as.numeric(c(t(outlier.brca.clinic['Overall.Survival..Months.']))),
  #  as.factor(c(t(outlier.brca.clinic['Progression.Free.Status']))))



outlier.brca.clinic$Subtype <- as.factor(c(t(outlier.brca.clinic['Subtype'])))
outlier.brca.clinic$Subtype <- relevel(outlier.brca.clinic$Subtype,ref = 'BRCA_LumA')
Subgroups <- as.factor(c(t(outlier.brca.clinic['Subtype'])))
Subgroups <- relevel(Subgroups,ref = 'BRCA_LumA')

fit.coxmodel(
    groups = as.factor(Outlier),
    survobj = survobj,
    return.cox.model = TRUE,
    );
fit.coxmodel(
    groups = as.factor(Subgroups),
    survobj = survobj,
    return.cox.model = TRUE,
    );

#Basal
basal <- data.frame(outlier.brca.clinic[outlier.brca.clinic$Subtype == 'BRCA_Basal',])
basal.censor.binary <- replace(
    x = basal$Overall.Survival.Status,
    list = ('0:LIVING' == basal$Overall.Survival.Status),
    values = 0
    );
basal.censor.binary <- replace(
    x = basal.censor.binary,
    list = ('1:DECEASED' == basal.censor.binary),
    values = 1
    );
basal.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = basal.censor.binary,
    groups = 'Outlier.Subgroups',
    data = basal,
    reference.group = '0'
    );
cox.zph(basal.cox)
#HER2
Her2 <- data.frame(outlier.brca.clinic[outlier.brca.clinic$Subtype == 'BRCA_Her2',])
her2.censor.binary <- replace(
    x = Her2$Overall.Survival.Status,
    list = ('0:LIVING' == Her2$Overall.Survival.Status),
    values = 0
    );
her2.censor.binary <- replace(
    x = her2.censor.binary,
    list = ('1:DECEASED' == her2.censor.binary),
    values = 1
    );
her2.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = her2.censor.binary,
    groups = 'Outlier.Subgroups',
    data = Her2,
    reference.group = '0'
    );
cox.zph(her2.cox)
#LUMA
LumA <- data.frame(outlier.brca.clinic[outlier.brca.clinic$Subtype == 'BRCA_LumA',])
LumA.censor.binary <- replace(
    x = LumA$Overall.Survival.Status,
    list = ('0:LIVING' == LumA$Overall.Survival.Status),
    values = 0
    );
LumA.censor.binary <- replace(
    x = LumA.censor.binary,
    list = ('1:DECEASED' == LumA.censor.binary),
    values = 1
    );
LumA.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = LumA.censor.binary,
    groups = 'Outlier.Subgroups',
    data = LumA,
    reference.group = '0'
    );
cox.zph(LumA.cox)
#LUMB
LumB <- data.frame(outlier.brca.clinic[outlier.brca.clinic$Subtype == 'BRCA_LumB',])
LumB.censor.binary <- replace(
    x = LumB$Overall.Survival.Status,
    list = ('0:LIVING' == LumB$Overall.Survival.Status),
    values = 0
    );
LumB.censor.binary <- replace(
    x = LumB.censor.binary,
    list = ('1:DECEASED' == LumB.censor.binary),
    values = 1
    );
LumB.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = LumB.censor.binary,
    groups = 'Outlier.Subgroups',
    data = LumB,
    reference.group = '0'
    );
cox.zph(LumB.cox)
#NORMAL
Normal <- data.frame(outlier.brca.clinic[outlier.brca.clinic$Subtype == 'BRCA_Normal',])
norm.censor.binary <- replace(
    x = Normal$Overall.Survival.Status,
    list = ('0:LIVING' == Normal$Overall.Survival.Status),
    values = 0
    );
norm.censor.binary <- replace(
    x = norm.censor.binary,
    list = ('1:DECEASED' == norm.censor.binary),
    values = 1
    );
normal.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = norm.censor.binary,
    groups = 'Outlier.Subgroups',
    data = Normal,
    reference.group = '0'
    );
cox.zph(normal.cox)#Problem

