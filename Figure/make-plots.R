#!/usr/bin/env Rscript
library(BoutrosLab.utilities);

message('Starting up');

# Helper functions to colorize text
ansi.yellow <- function(str) {
    paste('\x1b[33;1m', str, '\x1b[0m', sep = '');
    }
ansi.green <- function(str) {
    paste('\x1b[32;1m', str, '\x1b[0m', sep = '');
    }
ansi.red <- function(str) {
    paste('\x1b[31;1m', str, '\x1b[0m', sep = '');
    }
ansi.magenta <- function(str) {
    paste('\x1b[35;1m', str, '\x1b[0m', sep = '');
    }

save.multiple.plots <- function(subfiles, datafile, output.directory = 'figures') {
    dir.create(output.directory, showWarnings = FALSE);

    message(ansi.green(paste('Loading', datafile, 'into new environment')));
    data.env <- new.env();
    load(datafile, envir = data.env);

    result <- with(list(temp.varnames = ls(data.env)), {
        list(
            datafile = datafile,
            variables = setNames(lapply(temp.varnames, function(x) digest::digest(get(x, envir = data.env))), temp.varnames),
            figures = vector(mode = 'list', length = length(subfiles))
            )
        });

    for (figure.index in seq_along(subfiles)) {
        figure.file <- subfiles[figure.index];
        message(ansi.yellow(figure.file));
        # create.histogram resets warn to 0, so continually set it back
        options(warn = 1);

        accessed.vars <- new.env(parent = emptyenv());

        plot.env <- new.env(parent = data.env);
        plot.sub.env <- new.env(parent = plot.env);

        track.access <- function(name) {
            makeActiveBinding(
                name,
                function(value) {
                    if (missing(value)) {
                        assign(name, 'ACCESSED', envir = accessed.vars)
                        get(name, envir = plot.env)
                        } else {
                        assign(name, value, envir = plot.env)
                        }
                    },
                plot.sub.env
                )
            }

        for (var in ls(data.env)) {
            track.access(var)
            }

        if (inherits(
            try(source(figure.file, local = plot.sub.env, echo = FALSE)),
            'try-error'
            )) {
            message(ansi.red(paste('Problem with', figure.file, '!')));
            }

        for (variable.name in ls(plot.env)) {
            if (exists(variable.name, data.env)) {
                if (digest::digest(get(variable.name, data.env)) == digest::digest(get(variable.name, plot.env))) {
                    assign(variable.name, 'REDUNDANT', envir = accessed.vars);
                    } else {
                    assign(variable.name, paste0('CHANGED:', digest::digest(get(variable.name, envir = plot.env))), envir = accessed.vars);
                    }
                }
            }

        for (variable.name in ls(plot.sub.env)) {
            if (!exists(variable.name, data.env)) {
                assign(variable.name, paste0('DEFINED:', digest::digest(get(variable.name, envir = plot.sub.env))), envir = accessed.vars);
                }
            }


        result$figures[[figure.index]] <- mget(ls(accessed.vars, all.names = TRUE), envir = accessed.vars);
        names(result$figures)[figure.index] <- tools::file_path_sans_ext(basename(figure.file));
        }

    write(
        x = rjson::toJSON(result, indent = 2),
        file = paste0(
            tools::file_path_sans_ext(basename(datafile)),
            '.json'
            )
        )
    }
# Use the original data on the cluster, falling back to the copied local data

# The processed / subsetted data intended for plotting
data.dir <- '/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/data';
if (!dir.exists(data.dir)) {
    data.dir <- 'untracked_data/data';
    }

# The original data
dev.dir <- '/hot/users/jyhan/TCGA/RNA-seq/outlier';
if (!dir.exists(dev.dir)) {
    dev.dir <- 'untracked_data/outlier';
    }

save.multiple.plots(
    c(
        'Figure1/Figure1b.R',
        'Figure1/Figure1c.R',
        'Figure1/Figure1d.R',
        'Figure1/Figure1e.R',
        'Figure1/Figure1fg.R',
        'Figure1/Figure1i.R'
        ),
    file.path(data.dir, '2024-09-10_Figure1.rda')
    );


save.multiple.plots(
    c(
        'Figure2/Figure2a.R',
        'Figure2/Figure2b.R',
        'Figure2/Figure2c.R',
        'Figure2/Figure2d.R'
        ),
    file.path(data.dir, '2024-08-23_Figure2a-d.rda')
    );

save.multiple.plots(
    c(
        'Figure2/Figure2e.R',
        'Figure2/Figure2f.R'
        ),
    file.path(data.dir, '2024-09-10_Figure2ef_drivergene.rda')
    );

save.multiple.plots(
    c(
        'Figure2/Figure2h.R',
        'Figure2/Figure2i.R',
        'Figure2/Figure2j.R',
        'Figure2/Figure2k.R',
        'Figure2/Figure2l.R'
        ),
    file.path(data.dir, '2024-08-26_Figure2h-l_input.rda')
    );

save.multiple.plots(
    c(
        'Figure3/Figure3a.R',
        'Figure3/Figure3b.R',
        'Figure3/Figure3c.R',
        'Figure3/Figure3d.R'
        ),
    file.path(data.dir, '2024-09-10_Figure3a-d.rda')
    );

save.multiple.plots(
    c(
        'Figure3/Figure3e.R',
        'Figure3/Figure3f.R',
        'Figure3/Figure3g.R',
        'Figure3/Figure3h.R',
        'Figure3/Figure3i.R'
        ),
    file.path(data.dir, '2024-09-11_Figure3e-i.rda')
    );

save.multiple.plots(
    c(
        'Figure4/Figure4a.R',
        'Figure4/Figure4b.R',
        'Figure4/Figure4c.R',
        'Figure4/Figure4d.R',
        'Figure4/Figure4e.R',
        'Figure4/Figure4f.R',
        'Figure4/Figure4g.R',
        'Figure4/Figure4h.R',
        'Figure4/Figure4i.R',
        'Figure4/Figure4j.R',
        'Figure4/Figure4k.R',
        'Figure4/Figure4l.R',
        'Figure4/Figure4m.R',
        'Figure4/Figure4n.R'
        ),
    file.path(data.dir, '2024-09-10_Figure4.rda')
    );
