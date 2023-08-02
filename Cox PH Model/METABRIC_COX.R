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

### cox.uni function ##############################################################################
#Description
# calculating cox model for univariate stats
#Input Variables
    # (1) survival.time
    #   column name in data that contains survival time
    # (2) survival.event
    #   binary form of event 
    # (3) groups
    #   column name of variable in data for comparison
    # (4) data
    #   name of data frame with clinical variables
    # (5) reference.group
    #   value of variable within grroups that is used as reference
#Output Variable
    # (1) cox
    #   result from fit.coxmodel from BoutrosLab.statistics.survival
cox.uni <- function(survival.time,survival.event,groups,data,reference.group,cox.model = TRUE) {
    
    survobj <- Surv(
        as.numeric(c(t(data[survival.time]))),
        as.numeric(survival.event)
        );
    group <- as.factor(c(t(data[groups])));
    group <- relevel(group,ref = reference.group);
    
    cox <- fit.coxmodel(
        groups = as.factor(group),
        survobj = survobj,
        return.cox.model = cox.model,
        );
    return(cox);
    }

### mod.cox.forest function #######################################################################
#Description
# reconstructing coxph model inorder to utilize forest.plot package
#Input Variables
# (1) data
#   data object that contains a coxph 
# (2) subtype.name
#   name of the subtype
#Output Variable
# (1) new.mod
#   result of reconstruction for plotting in forest.plot
mod.cox.forest <- function(data,subtype.name) {
    mod <- broom::tidy(data);
    term <- HR <- lb <- ub <- p <- vector(length = nrow(data));
    for (i  in 1:nrow(data)) {
        term[i] = mod$term[i]
        HR[i] = round(mod$estimate[i],4);
        lb[i] = mod$estimate[i] - 1.96 * mod$std.error[i];
        ub[i] = mod$estimate[i] + 1.96 * mod$std.error[i];
        p[i] = round(mod$p.value[i],4);
        } 
    new.mod <- data.frame(
        Subtypes = subtype.name,
        feature = term,
        coef = as.numeric(HR),
        lower95.coef = as.numeric(lb),
        higher95.coef = as.numeric(ub),
        p = as.numeric(p)
        );
    return(new.mod);
    }
### Data Analysis #################################################################################
outlier.meta.clinic.patient <- modified.clinic(
    data = meta.clinic.patient,
    outlier.data = outlier.patient.tag.01,
    TCGA = FALSE,
    META = TRUE,
    sample = FALSE)
outlier.meta.progress.binary <- replace(
    x = outlier.meta.clinic.patient$Overall.Survival.Status,
    list = ('0:LIVING' == outlier.meta.clinic.patient$Overall.Survival.Status),
    values = 0
    );
outlier.meta.progress.binary <- replace(
    x = outlier.meta.progress.binary,
    list = ('1:DECEASED' == outlier.meta.progress.binary),
    values = 1
    );

meta.overall.outliers.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = outlier.meta.progress.binary,
    groups = 'Outlier.Subgroups',
    data = outlier.meta.clinic.patient,
    reference.group = '0',
    cox.model = TRUE 
    );
meta.overall.outliers.ph <- cox.zph(meta.overall.outliers.cox)
#LumA reference group 
meta.overall.subtypes.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = outlier.meta.progress.binary,
    groups = 'Pam50...Claudin.low.subtype',
    data = outlier.meta.clinic.patient,
    reference.group = 'LumA'
    );
meta.overall.subtypes.ph <- cox.zph(meta.overall.subtypes.cox)

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
meta.basal.ph <- cox.zph(meta.basal.cox);

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
meta.her2.ph <- cox.zph(meta.her2.cox);

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
meta.LumA.ph <- cox.zph(meta.LumA.cox);

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
meta.LumB.ph <- cox.zph(meta.LumB.cox);

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
meta.normal.ph <- cox.zph(meta.norm.cox);

meta.assumption <- list(
    Overall.Outliers.model = meta.overall.outliers.cox,
    Overall.Outliers.assumption = meta.overall.outliers.ph,
    Overall.Subtypes.model = meta.overall.subtypes.cox,
    Overall.Subtypes.assumption = meta.overall.subtypes.ph,
    Basal.model = meta.basal.cox,
    Basal.assumption = meta.basal.ph,
    Her2.model = meta.her2.cox,
    Her2.assumption = meta.her2.ph,
    LumA.model = meta.LumA.cox,
    LumA.assumption = meta.LumA.ph,
    LumB.model = meta.LumB.cox,
    LumB.assumption = meta.LumB.ph,
    Normal.model = meta.norm.cox,
    Normal.assumption = meta.normal.ph
    );

meta.requested.model <- coxph(formula = Surv(as.numeric(c(t(Overall.Survival..Months.))),as.numeric(outlier.meta.progress.binary)) ~ Outlier.Subgroups + Pam50...Claudin.low.subtype + Outlier.Subgroups * Pam50...Claudin.low.subtype, data=outlier.meta.patient)
cox.zph(meta.requested.model)

meta.basal.mod.forest <- mod.cox.forest(meta.basal.cox,'Basal')
meta.her2.mod.forest <- mod.cox.forest(meta.her2.cox,'Her2')
meta.lumA.mod.forest <- mod.cox.forest(meta.LumA.cox,'LumA')
meta.LumB.mod.forest <- mod.cox.forest(meta.LumB.cox,'LumB')
meta.Normal.mod.forest <- mod.cox.forest(meta.norm.cox,'Normal')

all.Meta.Subytpes <- rbind(basal.mod.forest,her2.mod.forest,lumA.mod.forest,LumB.mod.forest,Normal.mod.forest)

multi.forest.func(all.TCGA.Subytpes,'wwwii')

