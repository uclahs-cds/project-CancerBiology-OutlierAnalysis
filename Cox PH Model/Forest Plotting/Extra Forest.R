mod <- broom::tidy(her2.cox)
term <- HR <- lb <- ub <- p <- vector()
for (i  in 1:3) {
    term[i] = mod$term[i]
    HR[i] = round(mod$estimate[i],4)
    lb[i] = mod$estimate[i] - 1.96 * mod$std.error[i]
    ub[i] = mod$estimate[i] + 1.96 * mod$std.error[i]
    p[i] = round(mod$p.value[i],4)
} 
her2.mod <- data.frame(feature = term,coef = as.numeric(HR),lower95.coef = as.numeric(lb),higher95.coef = as.numeric(ub),p = as.numeric(p))

subgroup.colors <- c('Basal' = 'grey50','Her2' = 'white','grey50','white','grey50')
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
