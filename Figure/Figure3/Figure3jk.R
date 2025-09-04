### HISTORY #####################################################################
# This script tests whether XEGs are
# associated with survival outcomes across combined datasets.
# It computes Cox models (unadjusted and subtype-adjusted), exact Cox tests,
# log-rank tests, and visualizes significant genes with CI segment plots and KM.
# Date: 2024-08-16
################################################################################

### DESCRIPTION #################################################################
# This script processes outlier indicators for genes and tests survival
# associations using overall survival (OS). It:
# 1) Iterates through frequently outliered genes.
# 2) Fits Cox PH models and Exact Cox, derives HR/CI/p-values.
# 3) Computes log-rank p-values from KM curves.
# 4) FDR and visualizes significant genes with a segplot,
#    and an example Kaplan-Meier (KM) plot for example gene.
################################################################################


### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.survival);
library(BoutrosLab.utilities);
library(survminer);
library(ExactCox);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'outlier.patient.tag.01.brca.match.five',
    'outlier.patient.tag.01.meta.match.five',
    'outlier.symbol',
    'os.group.combine'
    ));


# Filtering criteria

outlier.patient.all.two.01 <- data.frame(
    cbind(
        outlier.patient.tag.01.brca.match.five,
        outlier.patient.tag.01.meta.match.five
        )
    );

# Set row names of the combined data frame to unique outlier symbols
rownames(outlier.patient.all.two.01) <- outlier.symbol$unique
all.cases <- outlier.patient.all.two.01[,rownames(os.group.combine)];
all.cases.gene <- list();

for(col in colnames(all.cases)) {
    all.cases.gene[[col]] <- rownames(all.cases)[all.cases[, col] %in% 1];
    }


all_cases.gene.table <- table(unlist(all.cases.gene));
frequent.gene <- all_cases.gene.table[all_cases.gene.table >= 5];

frequent.gene.cox.hr <- NULL;
frequent.gene.cox.ci <- NULL;
frequent.gene.cox.p <- NULL;
frequent.gene.cox.all <- NULL;

frequent.gene.cox.var.hr <- NULL;
frequent.gene.cox.var.ci <- NULL;
frequent.gene.cox.var.p <- NULL;
frequent.gene.cox.var.all <- NULL;

frequent.event <- NULL;
logrank.p <- NULL;
exact.cox.p <- NULL;

# Initialize containers referenced later
exact.cox.p.con <- NULL;

for (gene in names(frequent.gene)) {
    os.group.combine.gene <- os.group.combine;
    
    # Build outlier group (High: outlier present, Low: otherwise)
    outlier.patient.all.two.01.gene <- outlier.patient.all.two.01[gene, rownames(os.group.combine.gene)];
    os.group.combine.gene$out <- as.numeric(outlier.patient.all.two.01.gene);
    os.group.combine.gene$out.group <- ifelse(os.group.combine.gene$out >= 1, 'High', 'Low');
    
    # Omit missing rows
    os.group.combine.gene.na <- na.omit(os.group.combine.gene);
    
    # Count events in High group (for summary)
    frequent.event.num <- sum(os.group.combine.gene.na$status[os.group.combine.gene.na$out.group == 'High']);
    frequent.event <- c(frequent.event, frequent.event.num);
    
    # Relevel so that 'Low' is the reference
    os.group.combine.gene.na$out.group <- relevel(factor(os.group.combine.gene.na$out.group), ref = 'Low');
    
    # Cox models (unadjusted and subtype-adjusted)
    cox_fit <- coxph(Surv(os, status) ~ out.group, data = os.group.combine.gene.na);
    cox_fit.subtype <- coxph(Surv(os, status) ~ out.group + pam50, data = os.group.combine.gene.na);
    
    # Collect Cox (unadjusted)
    frequent.gene.cox.hr <- c(frequent.gene.cox.hr, summary(cox_fit)$coefficients[, 2]);
    frequent.gene.cox.ci <- rbind(frequent.gene.cox.ci, summary(cox_fit)$conf.int[, 3:4]);
    frequent.gene.cox.p <- c(frequent.gene.cox.p, summary(cox_fit)$coefficients[1, 'Pr(>|z|)']);
    frequent.gene.cox.all <- rbind(
        frequent.gene.cox.all, 
        c(
            summary(cox_fit)$coefficients[1, 2], 
            summary(cox_fit)$conf.int[1, 3:4], 
            summary(cox_fit)$coefficients[1, 'Pr(>|z|)']
        )
        );
    
    # Collect Cox (subtype-adjusted)
    frequent.gene.cox.var.hr <- c(frequent.gene.cox.var.hr, summary(cox_fit.subtype)$coefficients[1, 2]);
    frequent.gene.cox.var.ci <- rbind(frequent.gene.cox.var.ci, summary(cox_fit.subtype)$conf.int[1, 3:4]);
    frequent.gene.cox.var.p <- c(frequent.gene.cox.var.p, summary(cox_fit.subtype)$coefficients[1, 'Pr(>|z|)']);
    frequent.gene.cox.var.all <- rbind(
        frequent.gene.cox.var.all, 
        c(
            summary(cox_fit.subtype)$coefficients[1, 2], 
            summary(cox_fit.subtype)$conf.int[1, 3:4], 
            summary(cox_fit.subtype)$coefficients[1, 'Pr(>|z|)']
        )
        );
    
    # KM and log-rank p-value
    surv_obj <- Surv(os.group.combine.gene.na$os, os.group.combine.gene.na$status);
    fit <- survfit(surv_obj ~ out.group, data = os.group.combine.gene.na);
    a <- surv_pvalue(fit, data = os.group.combine.gene.na)$pval;
    logrank.p <- c(logrank.p, a);
    
    # Exact Cox
    exact.cox <- ExactCox(
        time = os.group.combine.gene.na$os,
        status = os.group.combine.gene.na$status,
        trt = os.group.combine.gene.na$out.group,
        conf.int = TRUE
        );
    exact.cox.p <- c(exact.cox.p, exact.cox$p.value);
    exact.cox.p.con <- rbind(exact.cox.p.con, as.numeric(exact.cox$conf.int));
    };

# Assemble unadjusted Cox summary table
frequent.gene.cox.all.df <- data.frame(frequent.gene.cox.all);
colnames(frequent.gene.cox.all.df) <- c('HR', 'ci.low', 'ci.max', 'p.value');
rownames(frequent.gene.cox.all.df) <- names(frequent.gene);

# Adjust p-values for unadjusted Cox (BH)
frequent.gene.cox.p.new <- p.adjust(frequent.gene.cox.p, method = 'BH');

# Store additional stats
frequent.gene.cox.all.df$new.p.new <- frequent.gene.cox.p.new;
frequent.gene.cox.all.df$event <- frequent.event;
frequent.gene.cox.all.df$exactcox <- exact.cox.p;

# Assemble subtype-adjusted Cox summary table
frequent.gene.cox.var.all.df <- data.frame(frequent.gene.cox.var.all);
colnames(frequent.gene.cox.var.all.df) <- c('var_HR', 'var_ci.low', 'var_ci.max', 'var_p.value');
rownames(frequent.gene.cox.var.all.df) <- names(frequent.gene);
frequent.gene.cox.var.all.df$var_fdr <- p.adjust(frequent.gene.cox.var.all.df$var_p.value, method = 'BH');

# Attach log-rank p-values and FDR (row alignment by gene names)
frequent.gene.cox.all.df.filter <- data.frame(log.rank = logrank.p);
rownames(frequent.gene.cox.all.df.filter) <- names(frequent.gene);
frequent.gene.cox.all.df$log.rank <- frequent.gene.cox.all.df.filter$log.rank;
frequent.gene.cox.all.df$log.rank.fdr <- p.adjust(frequent.gene.cox.all.df.filter$log.rank, method = 'BH');

# Select genes by log-rank FDR < 0.1 and order by HR
frequent.gene.cox.all.01 <- frequent.gene.cox.all.df[frequent.gene.cox.all.df$log.rank.fdr < 0.1, ];
frequent.gene.cox.all.01.order <- frequent.gene.cox.all.01[order(frequent.gene.cox.all.01$HR), ];
frequent.gene.cox.all.01.order$labels <- as.factor(rownames(frequent.gene.cox.all.01.order));

# Set dot colours by direction and significance
dot.colours <- vector(length = nrow(frequent.gene.cox.all.01.order));
dot.colours <- rep('grey70', nrow(frequent.gene.cox.all.01.order));
dot.colours[frequent.gene.cox.all.01.order$log.rank.fdr < 0.1 & frequent.gene.cox.all.01.order$HR < 1] <- 'dodgerblue2';
dot.colours[frequent.gene.cox.all.01.order$log.rank.fdr < 0.1 & frequent.gene.cox.all.01.order$HR > 1] <- 'red2';

# Segplot of HR with CI (log2 scale)
frequent.gene.os <- BoutrosLab.plotting.general::create.segplot(
    formula = labels ~ log2(ci.low) + log2(ci.max),
    data = frequent.gene.cox.all.01.order,
    centers = log2(frequent.gene.cox.all.01.order$HR),
    main.cex = 0,
    ylab.label = NULL,
    yaxis.lab = frequent.gene.cox.all.01.order$labels,
    yaxis.fontface = 1,
    xlab.cex = 1,
    xlab.label = expression('HR'),
    ylab.cex = 1,
    yaxis.cex = 1,
    xaxis.cex = 1,
    xlimits = c(-3.5, 5.3),
    xat = c(-2, 0, 2, 4),
    xaxis.lab = c(2^-2, 2^0, 2^2, 2^4),
    xaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    segments.col = dot.colours,
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -7,
    xright.rectangle = 7,
    ybottom.rectangle = seq(1.5, 63.5, 2),
    ytop.rectangle = seq(2.5, 64.5, 2),
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
    );
frequent.gene.os;


save.outlier.figure(
    frequent.gene.os,
    c('Figure3j', 'XEG_survival', 'subtype', 'segment'),
    width = 5,
    height = 3.8
    );
    

# Example: Single gene KM based solely on outlier status
gene <- 'TBX2';
os.group.combine.gene <- os.group.combine;

# Build outlier group for the example gene
outlier.patient.all.two.01.gene <- outlier.patient.all.two.01[gene, rownames(os.group.combine.gene)];
os.group.combine.gene$out <- as.numeric(outlier.patient.all.two.01.gene);
os.group.combine.gene$out.group <- ifelse(os.group.combine.gene$out >= 1, 'High', 'Low');

# Omit missing and relevel
os.group.combine.gene.na <- na.omit(os.group.combine.gene);
os.group.combine.gene.na$out.group <- relevel(factor(os.group.combine.gene.na$out.group), ref = 'Low');

# Cox model for the example gene
cox_fit <- coxph(Surv(os, status) ~ out.group, data = os.group.combine.gene.na);

# KM fit for the example gene
fit <- survfit(Surv(os, status) ~ out.group, data = os.group.combine.gene.na);

cox_summary <- summary(cox_fit);

# Log-rank test for p-value
logrank_test <- survdiff(Surv(os.group.combine.gene.na$os, os.group.combine.gene.na$status) ~ 
        as.factor(os.group.combine.gene.na$out.group))
logrank_p <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)

# Kaplan-Meier survival plot
km.os.group.gene <- create.km.plot(
    survival.object = Surv(os.group.combine.gene.na$os, os.group.combine.gene.na$status),
    main = as.expression(substitute(paste('Kaplan-Meier estimate (Combined datasets)'))),
    show.risktable = TRUE,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    xlab.label = expression('Overall survival (Months)'),
    ylab.label = expression('Estimated proportion'),
    xlimits = c(0, 280),
    xat = seq(0, 240, 80),
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    main.cex = 1.5,
    statistical.method = 'cox',
    predefined.hr = round(cox_summary$conf.int[1, "exp(coef)"], digits = 2),  # Cox HR
    predefined.hr.ci = round(cox_summary$conf.int[1, c("lower .95", "upper .95")], digits = 2),  # Cox CI
    predefined.p = logrank_p,  # Log-rank p-value
    key.stats.cex = 1.1,
    patient.groups = as.factor(os.group.combine.gene.na$out.group),
    risktable.fontsize = 11.5,
    show.key.groups = TRUE,
    risk.label.pos = -70,
    ylab.axis.padding = 2,
    risk.label.fontface = 1,
    left.padding = 5.5,
    key.groups.labels = rev(c('Outlier patients', 'Non-outlier patients')),
    key.groups.cex = 1,
    line.colours = rev(c('red3', 'dodgerblue3'))
    );
km.os.group.gene;


save.outlier.figure(
    km.os.group.gene,
    c('Figure3k', gene, 'km', 'survival'),
    width = 5,
    height = 5
    );


save.session.profile(file.path('output', 'Figure3jk.txt'));
