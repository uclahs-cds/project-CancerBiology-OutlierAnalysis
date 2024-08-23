invisible(loadNamespace(library(logger)));
library(BoutrosLab.utilities);

logger::log_threshold(DEBUG);

logger::log_appender(logger::appender_file('plotting.log'));
logger::log_errors();
logger::log_warnings();
logger::log_messages();

logger::log_info('Starting up');
load('2024-08-23_Figure1.rda');

subfiles <- c(
    'Figure/Figure1/Figure1c.R',
    'Figure/Figure1/Figure1d.R',
    'Figure/Figure1/Figure1e.R'
);

for (i in seq_along(subfiles)) {
    message(subfiles[i]);
    source(subfiles[i]);
}

logger::log_info('Shutting down');
