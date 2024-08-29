#!/usr/bin/env Rscript
invisible(loadNamespace(library(logger)));
library(BoutrosLab.utilities);

logger::log_threshold(DEBUG);

logger::log_appender(logger::appender_file('plotting.log'));
logger::log_errors();
logger::log_warnings();
logger::log_messages();

logger::log_info('Starting up');

# Upgrade all warnings to errors
options(warn = 2);

ansi.yellow <- function(str) {
    paste('\x1b[33;1m', str, '\x1b[0m', sep = '');
}

ansi.green <- function(str) {
    paste('\x1b[32;1m', str, '\x1b[0m', sep = '');
}

ansi.red <- function(str) {
    paste('\x1b[31;1m', str, '\x1b[0m', sep = '');
}

save.multiple.plots <- function(datafile, subfiles) {
    message(ansi.green(paste('Sourcing', datafile)));

    data.env.name <- attr(attach(datafile), 'name');

    for (figure.file in subfiles) {
        message(ansi.yellow(figure.file));
        # create.histogram resets warn to 0, so continually set it back
        options(warn = 2);

        if (inherits(
            try(source(figure.file, local = new.env(), echo = FALSE)),
            'try-error'
        )) {
            message(ansi.red(paste('Problem with', figure.file, '!')));
        }
        warnings();
    }

    message(ansi.green(paste('Detaching', datafile)));
    detach(data.env.name, character.only = TRUE);
}

save.multiple.plots(
    'untracked_data/outlier/2024-08-27_metabric_tcga_ispy_matador_icgc.RData',
    c(
        'Figure/Figure1/Figure1b.R',
        'Figure/Figure1/Figure1c.R',
        'Figure/Figure1/Figure1d.R',
        'Figure/Figure1/Figure1e.R',
        'Figure/Figure1/Figure1h.R',
        'Figure/Figure1/Figure1i.R'
    )
);

save.multiple.plots(
    'untracked_data/outlier/2024-08-27_cnv_all_brca_meta_icgc.RData',
    c(
        'Figure/Figure2/Figure2a.R',
        'Figure/Figure2/Figure2b.R',
        'Figure/Figure2/Figure2c.R',
        'Figure/Figure2/Figure2d.R'
    )
);

save.multiple.plots(
    'untracked_data/outlier/2024-05-05_driver_gene.RData',
    c(
        'Figure/Figure2/Figure2e.R',
        'Figure/Figure2/Figure2f.R'
    )
);

save.multiple.plots(
    'untracked_data/outlier/2024-08-26_meta_brca_methylation_merge.RData',
    c(
        'Figure/Figure2/Figure2h.R'
    )
);
