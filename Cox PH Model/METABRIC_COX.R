### METABRIC COX #########################################################################

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
### Example
# will have to do outside of function
outlier.meta.progress.binary <- replace(
    x = outlier.meta.patient$Overall.Survival.Status,
    list = ('0:LIVING' == outlier.meta.patient$Overall.Survival.Status),
    values = 0
);
outlier.meta.progress.binary <- replace(
    x = outlier.progress.binary,
    list = ('1:DECEASED' == outlier.progress.binary),
    values = 1
);

overall.outliers.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = outlier.progress.binary,
    groups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    reference.group = '0'
    );

#LumA reference group 
overall.subtypes.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = outlier.progress.binary,
    groups = 'Subtype',
    data = outlier.brca.clinic,
    reference.group = 'BRCA_LumA'
    );


### Example
# will have to do outside of function
outlier.meta.progress.binary <- replace(
    x = outlier.meta.patient$Overall.Survival.Status,
    list = ('0:LIVING' == outlier.meta.patient$Overall.Survival.Status),
    values = 0
);
outlier.meta.progress.binary <- replace(
    x = outlier.progress.binary,
    list = ('1:DECEASED' == outlier.progress.binary),
    values = 1
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
meta.basal <- data.frame(outlier.meta.patient[outlier.meta.patient$Pam50...Claudin.low.subtype == 'Basal',])
meta.basal.censor.binary <- replace(
    x = meta.basal$Overall.Survival.Status,
    list = ('0:LIVING' == meta.basal$Overall.Survival.Status),
    values = 0
    );
meta.basal.censor.binary <- replace(
    x = meta.basal.censor.binary,
    list = ('1:DECEASED' == meta.basal.censor.binary),
    values = 1
    );
meta.basal.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = meta.basal.censor.binary,
    groups = 'Outlier.Subgroups',
    data = meta.basal,
    reference.group = '0'
    );
cox.zph(meta.basal.cox)
#HER2
meta.her2 <- data.frame(outlier.meta.patient[outlier.meta.patient$Pam50...Claudin.low.subtype == 'Her2',])
meta.her2.censor.binary <- replace(
    x = meta.her2$Overall.Survival.Status,
    list = ('0:LIVING' == meta.her2$Overall.Survival.Status),
    values = 0
    );
meta.her2.censor.binary <- replace(
    x = meta.her2.censor.binary,
    list = ('1:DECEASED' == meta.her2.censor.binary),
    values = 1
    );
meta.her2.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = meta.her2.censor.binary,
    groups = 'Outlier.Subgroups',
    data = meta.her2,
    reference.group = '0'
    );
cox.zph(meta.her2.cox)
#LUMA
meta.LumA <- data.frame(outlier.meta.patient[outlier.meta.patient$Pam50...Claudin.low.subtype == 'LumA',])
meta.LumA.censor.binary <- replace(
    x = meta.LumA$Overall.Survival.Status,
    list = ('0:LIVING' == meta.LumA$Overall.Survival.Status),
    values = 0
    );
meta.LumA.censor.binary <- replace(
    x = meta.LumA.censor.binary,
    list = ('1:DECEASED' == meta.LumA.censor.binary),
    values = 1
    );
meta.LumA.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = meta.LumA.censor.binary,
    groups = 'Outlier.Subgroups',
    data = meta.LumA,
    reference.group = '0'
    );
cox.zph(meta.LumA.cox)
#LUMB
meta.LumB <- data.frame(outlier.meta.patient[outlier.meta.patient$Pam50...Claudin.low.subtype == 'LumB',])
meta.LumB.censor.binary <- replace(
    x = meta.LumB$Overall.Survival.Status,
    list = ('0:LIVING' == meta.LumB$Overall.Survival.Status),
    values = 0
    );
meta.LumB.censor.binary <- replace(
    x = meta.LumB.censor.binary,
    list = ('1:DECEASED' == meta.LumB.censor.binary),
    values = 1
    );
meta.LumB.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = meta.LumB.censor.binary,
    groups = 'Outlier.Subgroups',
    data = meta.LumB,
    reference.group = '0'
    );
cox.zph(meta.LumB.cox)
#NORMAL
meta.norm <- data.frame(outlier.meta.patient[outlier.meta.patient$Pam50...Claudin.low.subtype == 'Normal',])
meta.norm.censor.binary <- replace(
    x = meta.norm$Overall.Survival.Status,
    list = ('0:LIVING' == meta.norm$Overall.Survival.Status),
    values = 0
    );
meta.norm.censor.binary <- replace(
    x = meta.norm.censor.binary,
    list = ('1:DECEASED' == meta.norm.censor.binary),
    values = 1
    );
meta.norm.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = meta.norm.censor.binary,
    groups = 'Outlier.Subgroups',
    data = meta.norm,
    reference.group = '0'
    );
cox.zph(meta.norm.cox)
