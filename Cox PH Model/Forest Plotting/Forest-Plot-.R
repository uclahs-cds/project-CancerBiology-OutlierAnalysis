### HR Forest Plots ###############################################################################

### Preamble ######################################################################################
library('forestplot');
library('tibble');
library('dplyr');


### forest.func ###################################################################################
#Description 
#
#Input Variables
    # (1) cox.mod
    #
    # (2) filename
    #
#Output Variable
    # (1) plot
    #
forest.func <- function(cox.mod, filename){
    mod <- broom::tidy(cox.mod)
    term <- HR <- lb <- ub <- p <- vector(length = nrow(mod))
    for (i  in 1:nrow(mod)) {
        term[i] = mod$term[i]
        HR[i] = format(x = mod$estimate[i],digits = 4)
        lb[i] = mod$estimate[i] - 1.96 * mod$std.error[i]
        ub[i] = mod$estimate[i] + 1.96 * mod$std.error[i]
        p[i] = format(x = mod$p.value[i],digits = 4)
    } 
    data <- data.frame(feature = term,coef = as.numeric(HR),lower95.coef = as.numeric(lb),higher95.coef = as.numeric(ub),p = as.numeric(p))
    forest.data <- tibble(
        mean = as.numeric(data$coef),
        lower = as.numeric(data$lower95.coef),
        upper = as.numeric(data$higher95.coef),
        Features = data$feature,
        log2HR = as.character(data$coef),
        Pvalue = as.character(data$p)
        );
    header <- tibble(
        Features = c('Features'),
        log2HR = c('log2(HR)'),
        Pvalue = c('P value')
        );
    plot.df <- bind_rows(header, forest.data);
    png(file = paste0(filename,'.png'), width = 150, height = 10 * nrow(forest.data) + 5, units='mm', res = 300);
    print(
        plot.df %>%
            forestplot(
                labeltext = c(Features, log2HR, Pvalue),
                lineheight = unit(1,'cm'),
                colgap = unit(3,'mm'),
                lwd.ci = 2,
                boxsize = 0.2,
                xlab = expression('log'[2]~'(HR)'),
                ci.vertices = TRUE,
                clip = c(-15, 15),
                graphwidth = unit(6,'cm'),
                )
        );
    dev.off();
    return();
    }

### multi.forest.func #############################################################################
#Description 
#
#Input Variables
    # (1) data
    #
    # (2) filename
    #
#Output Variable
    # (1) plot
    #
multi.forest.func <- function(data, filename){
    forest.data <- tibble(
        mean = as.numeric(data$coef),
        lower = as.numeric(data$lower95.coef),
        upper = as.numeric(data$higher95.coef),
        Subtypes = as.character(data$Subtypes),
        Features = data$feature,
        log2HR = as.character(data$coef),
        Pvalue = as.character(data$p)
    );
    header <- tibble(
        Subtypes = c('Subtypes'),
        Features = c('Features'),
        log2HR = c('log2(HR)'),
        Pvalue = c('P value')
    );
    plot.df <- bind_rows(header, forest.data);
    png(file = paste0(filename,'.png'), width = 150, height = 10 * nrow(forest.data) + 5, units='mm', res = 300);
    print(
        plot.df %>%
            forestplot(
                labeltext = c(Subtypes,Features, log2HR, Pvalue),
                lineheight = unit(1,'cm'),
                colgap = unit(3,'mm'),
                lwd.ci = 2,
                boxsize = 0.2,
                #boxfill = subgroup.colors[forest_data$subgroup],
                xlab = expression('log'[2]~'(HR)'),
                ci.vertices = TRUE,
                clip = c(-15, 15),
                graphwidth = unit(6,'cm'),
            )
    );
    dev.off();
    return();
    }

### Data Analysis #################################################################################
#General plots
forest.func(
    cox.mod = overall.outliers.cox,
    filename = 'TCGA_Outlier_COX_Model'
    );
forest.func(
    cox.mod = overall.subtypes.cox,
    filename = 'TCGA_Subtypes_COX_Model'
    );
forest.func(
    cox.mod = meta.overall.outliers.cox,
    filename = 'METABRIC_Outlier_COX_Model'
    );
forest.func(
    cox.mod = meta.overall.subtypes.cox,
    filename = 'METABRIC_Subtypes_COX_Model'
    );
#Outlier per subtype plots
multi.forest.func(
    data = all.TCGA.Subytpes,
    filename = 'TCGA_Outlier-per-Subtype_COX_Model')

