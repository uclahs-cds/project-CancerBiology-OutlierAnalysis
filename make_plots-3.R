#!/usr/bin/env Rscript
invisible(loadNamespace(library(logger)));
library(BoutrosLab.utilities);

logger::log_threshold(DEBUG);

logger::log_appender(logger::appender_file('plotting.log'));
logger::log_errors();
logger::log_warnings();
logger::log_messages();

logger::log_info('Starting up');
load('/Users/nwiltsie/src/project-CancerBiology-OutlierAnalysis/untracked_data/outlier/2024-05-05_driver_gene.RData');

subfiles <- c(
    # 'Figure/Figure2/Figure2e.R',
    'Figure/Figure2/Figure2f.R'
);

for (i in seq_along(subfiles)) {
    message(paste('\x1b[33;1m', subfiles[i], '\x1b[0m', sep = ''));
    source(subfiles[i], echo = TRUE);
}

logger::log_info('Shutting down');
