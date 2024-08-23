library(logger);
library(BoutrosLab.utilities);

log_threshold(DEBUG);

log.filename <- paste0('plotting_', format(Sys.time(), '%Y%m%d_%H%M%S'), '.log');

log_appender(appender_tee(log.filename))
log_errors();
log_warnings();
log_messages();

log_info('Starting up');
load('2024-08-23_Figure1.rda');

source('Figure/Figure1/Figure1b.R');
