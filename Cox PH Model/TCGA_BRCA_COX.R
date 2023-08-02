### TCGA-BRCA COX #################################################################################

### Preamble ######################################################################################
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.prognosticsignature.general);
library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(BoutrosLab.statistics.survival);
install.packages("car")
library(car)
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
cox.uni <- function(survival.time,survival.event,groups,data,reference.group,cox.model=TRUE) {
    
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
        HR[i] = round(mod$estimate[i],4)
        lb[i] = mod$estimate[i] - 1.96 * mod$std.error[i]
        ub[i] = mod$estimate[i] + 1.96 * mod$std.error[i]
        p[i] = round(mod$p.value[i],4)
        } 
    new.mod <- data.frame(
        Subtypes = subtype.name,
        feature = term,
        coef = as.numeric(HR),
        lower95.coef = as.numeric(lb),
        higher95.coef = as.numeric(ub),
        p = as.numeric(p));
    return(new.mod);
    }

#after univariate
#surv() ~ outlier + subtype + outlier * subtype

### Data Analysis #################################################################################
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
#coxph(formula = Surv(as.numeric(c(t(Overall.Survival..Months.))),as.numeric(outlier.progress.binary)) ~ Outlier.Subgroups,data = outlier.brca.clinic)
overall.outliers.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = outlier.progress.binary,
    groups = 'Outlier.Subgroups',
    data = outlier.brca.clinic,
    reference.group = '0',
    cox.model = TRUE
    );
overall.outliers.ph <- cox.zph(overall.outliers.cox)

#LumA reference group 
overall.subtypes.cox <- cox.uni(
    survival.time = 'Overall.Survival..Months.',
    survival.event = outlier.progress.binary,
    groups = 'Subtype',
    data = outlier.brca.clinic,
    reference.group = 'BRCA_LumA',
    cox.model = TRUE
    );
overall.subtypes.ph <- cox.zph(overall.subtypes.cox)

### Outlier Subgroups BY Subtypes
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
    reference.group = '0',
    cox.model = TRUE
    );
basal.ph <- cox.zph(basal.cox)

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
her2.ph <- cox.zph(her2.cox)

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
LumA.ph <- cox.zph(LumA.cox)
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
LumB.ph <- cox.zph(LumB.cox)
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
#normal.ph <- 
    cox.zph(normal.cox)
#Problem
COX.ph <- list(
    Overall.Outliers.model = overall.outliers.cox,
    Overall.Outliers.assumption = overall.outliers.ph,
    Overall.Subtypes.model = overall.subtypes.cox,
    Overall.Subtypes.assumption = overall.subtypes.ph,
    Basal.model = basal.cox,
    Basal.assumption = basal.ph,
    Her2.model = her2.cox,
    Her2.assumption = her2.ph,
    LumA.model = LumA.cox,
    LumA.assumption = LumA.ph,
    LumB.model = LumB.cox,
    LumB.assumption = LumB.ph,
    Normal.model = normal.cox
    #Normal = normal.ph
    );


### MUTLIVARIATES
#surv() ~ outlier + subtype + outlier * subtype

requested.model <- coxph(
    formula = Surv(as.numeric(c(t(Overall.Survival..Months.))),as.numeric(outlier.progress.binary)) ~ Outlier.Subgroups + Subtype + Outlier.Subgroups * Subtype, 
    data=outlier.brca.clinic)

cox.zph(requested.model)

age.as.aconfounding <- coxph(formula = Surv(as.numeric(c(t(Overall.Survival..Months.))),as.numeric(outlier.progress.binary)) ~ Outlier.Subgroups + Diagnosis.Age,data = outlier.brca.clinic)
cox.zph(age.as.aconfounding)
vif(age.as.aconfounding)



mod.cox.forest <- function(data,subgroup.name) {
    mod <- broom::tidy(data)
    term <- HR <- lb <- ub <- p <- vector()
        for (i  in 1:3) {
            term[i] = mod$term[i]
            HR[i] = round(mod$estimate[i],4)
            lb[i] = mod$estimate[i] - 1.96 * mod$std.error[i]
            ub[i] = mod$estimate[i] + 1.96 * mod$std.error[i]
            p[i] = round(mod$p.value[i],4)
            } 
    new.mod <- data.frame(Subtypes = subgroup.name,feature = term,coef = as.numeric(HR),lower95.coef = as.numeric(lb),higher95.coef = as.numeric(ub),p = as.numeric(p))
    return(new.mod)
}

basal.mod.forest <- mod.cox.forest(basal.cox,'Basal')
her2.mod.forest <- mod.cox.forest(her2.cox,'Her2')
lumA.mod.forest <- mod.cox.forest(LumA.cox,'LumA')
LumB.mod.forest <- mod.cox.forest(LumB.cox,'LumB')
Normal.mod.forest <- mod.cox.forest(normal.cox,'Normal')

all.TCGA.Subytpes <- rbind(basal.mod.forest,her2.mod.forest,lumA.mod.forest,LumB.mod.forest,Normal.mod.forest)

multi.forest.func(all.TCGA.Subytpes,'wwwii')
