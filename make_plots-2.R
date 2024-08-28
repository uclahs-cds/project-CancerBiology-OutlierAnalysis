#!/usr/bin/env Rscript
invisible(loadNamespace(library(logger)));
library(BoutrosLab.utilities);

logger::log_threshold(DEBUG);

logger::log_appender(logger::appender_file('plotting.log'));
logger::log_errors();
logger::log_warnings();
logger::log_messages();

logger::log_info('Starting up');
load('/Users/nwiltsie/src/project-CancerBiology-OutlierAnalysis/untracked_data/outlier/2024-08-27_cnv_all_brca_meta_icgc.RData');

subfiles <- c(
    'Figure/Figure2/Figure2a.R',
    'Figure/Figure2/Figure2b.R',
    'Figure/Figure2/Figure2c.R',
    'Figure/Figure2/Figure2d.R',
    'Figure/Figure2/Figure2e.R',
    'Figure/Figure2/Figure2f.R',
    'Figure/Figure2/Figure2h.R',
    'Figure/Figure2/Figure2i.R',
    'Figure/Figure2/Figure2j.R',
    'Figure/Figure2/Figure2k.R',
    'Figure/Figure2/Figure2l.R'
);

for (i in seq_along(subfiles)) {
    message(paste('\x1b[33;1m', subfiles[i], '\x1b[0m', sep = ''));
    source(subfiles[i], echo = TRUE);
}

logger::log_info('Shutting down');
