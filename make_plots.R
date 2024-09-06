#!/usr/bin/env Rscript
invisible(loadNamespace(library(logger)));
library(BoutrosLab.utilities);

# Set up logging code
logger::log_threshold(DEBUG);
logger::log_appender(logger::appender_file('plotting.log'));
logger::log_errors();
logger::log_warnings();
logger::log_messages();

logger::log_info('Starting up');

# Upgrade all warnings to errors
options(warn = 2);

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

warn.env.duplicates <- function(base.env, plot.env) {
    message('Comparing dupes');
    base.stuff <- ls(base.env);
    local.stuff <- ls(plot.env);
    for (duplicate.object.name in intersect(base.stuff, local.stuff)) {
        if (digest::digest(as.environment(base.env)[[duplicate.object.name]]) ==
                digest::digest(as.environment(plot.env)[[duplicate.object.name]])) {
            message(ansi.green(paste('Redundant:', duplicate.object.name)));
        } else {
            message(ansi.red(paste('DIFFERENT:', duplicate.object.name)));
        }
    }
}

save.multiple.plots <- function(datafiles, subfiles, output.directory = 'figures') {
    data.env.names <- character(length(datafiles));

    dir.create(output.directory, showWarnings = FALSE);

    for (i in seq_along(datafiles)) {
        message(ansi.green(paste('Sourcing', datafiles[i])));
        # Don't fail on warnings while loading the datafile
        options(warn = 1);
        data.env.names[i] <- attr(attach(datafiles[i]), 'name');
        options(warn = 2);
    }

    for (figure.file in subfiles) {
        message(ansi.yellow(figure.file));
        # create.histogram resets warn to 0, so continually set it back
        options(warn = 1);

        accessed.vars <- new.env(parent = emptyenv());
        plot.env <- new.env();
        plot.sub.env <- new.env(parent = plot.env);

        track.access <- function(name) {
            makeActiveBinding(
                name,
                function(value) {
                    if (missing(value)) {
                        assign(name, get(name, envir = plot.env), envir = accessed.vars)
                        get(name, envir = plot.env)
                    } else {
                        assign(name, value, envir = plot.env)
                    }
                },
                plot.sub.env
            )
        }

        for (data.env.name in data.env.names) {
            for (var in ls(data.env.name)) {
                track.access(var)
            }
        }

        if (inherits(
            try(source(figure.file, local = plot.sub.env, echo = FALSE)),
            'try-error'
        )) {
            message(ansi.red(paste('Problem with', figure.file, '!')));
        }

        for (accessed.var in ls(accessed.vars)) {
            message(ansi.yellow(paste('Accessed', accessed.var)));
        }

        for (data.env.name in data.env.names) {
            warn.env.duplicates(data.env.name, plot.env);
        }
    }

    for (data.env.name in data.env.names) {
        message(ansi.green(paste('Detaching', data.env.name)));
        detach(data.env.name, character.only = TRUE);
    }
}


# Generate the full set of plots twice - once with the original dataset, once
# with the restricted dataset for plotting.
make.plots.twice <- function(figure.files, full.dataset, restricted.dataset) {
    save.multiple.plots(
        full.dataset,
        figure.files,
        'full_figures'
    );

    save.multiple.plots(
        restricted.dataset,
        figure.files,
        'restricted_figures'
    );
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

make.plots.twice(
    c(
        'Figure/Figure1/Figure1b.R',
        'Figure/Figure1/Figure1c.R',
        'Figure/Figure1/Figure1d.R',
        'Figure/Figure1/Figure1e.R',
        'Figure/Figure1/Figure1h.R',
        'Figure/Figure1/Figure1i.R'
    ),
    file.path(dev.dir, '2024-08-27_metabric_tcga_ispy_matador_icgc.RData'),
    file.path(data.dir, '2024-08-27_Figure1.rda')
);

make.plots.twice(
    c(
        'Figure/Figure2/Figure2a.R',
        'Figure/Figure2/Figure2b.R',
        'Figure/Figure2/Figure2c.R',
        'Figure/Figure2/Figure2d.R'
    ),
    file.path(dev.dir, '2024-08-27_cnv_all_brca_meta_icgc.RData'),
    file.path(data.dir, '2024-08-23_Figure2a-d.rda')
);

make.plots.twice(
    c(
        'Figure/Figure2/Figure2e.R',
        'Figure/Figure2/Figure2f.R'
    ),
    file.path(dev.dir, '2024-05-05_driver_gene.RData'),
    file.path(data.dir, '2024-08-24_Figure2ef_drivergene.rda')
);

make.plots.twice(
    c(
        'Figure/Figure2/Figure2h.R',
        'Figure/Figure2/Figure2i.R',
        'Figure/Figure2/Figure2j.R',
        'Figure/Figure2/Figure2k.R',
        'Figure/Figure2/Figure2l.R'
    ),
    file.path(dev.dir, '2024-08-26_meta_brca_methylation_merge.RData'),
    file.path(data.dir, '2024-08-26_Figure2h-l_input.rda')
);

make.plots.twice(
    c(
        'Figure/Figure3/Figure3a.R',
        'Figure/Figure3/Figure3b.R',
        'Figure/Figure3/Figure3c.R',
        'Figure/Figure3/Figure3d.R'
    ),
    file.path(dev.dir, '2024-02-20_brca.RData'),
    file.path(data.dir, '2024-08-28_Figure3a-d.rda')
);
