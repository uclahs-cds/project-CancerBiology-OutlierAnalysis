#!/usr/bin/env Rscript
invisible(loadNamespace(library(logger)));
library(BoutrosLab.utilities);

logger::log_threshold(DEBUG);

# logger::log_layout(logger::layout_glue_colors);
logger::log_appender(logger::appender_file('plotting.log'));
logger::log_errors();
logger::log_warnings();
logger::log_messages();

logger::log_info('Starting up');
# load('2024-08-23_Figure1.rda');
load('/Users/nwiltsie/src/project-CancerBiology-OutlierAnalysis/untracked_data/outlier/2024-08-27_metabric_tcga_ispy_matador_icgc.RData');

subfiles <- c(
    'Figure/Figure1/Figure1b.R',
    'Figure/Figure1/Figure1c.R',
    'Figure/Figure1/Figure1d.R',
    'Figure/Figure1/Figure1e.R',
    # 'Figure/Figure1/Figure1h.R', # No plotting code present
    'Figure/Figure1/Figure1i.R'
);

for (i in seq_along(subfiles)) {
    message(paste('\x1b[33;1m', subfiles[i], '\x1b[0m', sep = ''));
    source(subfiles[i], echo = TRUE);
}

logger::log_info('Shutting down');
