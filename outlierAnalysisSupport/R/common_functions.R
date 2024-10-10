### UTILITY FUNCTIONS ##############################################################################

### get.outlier.data.dir ##########################################################################
# Description:
#
# Return the path to the outlier data directory. The first existing value in
# the following list is returned:
#
#   * $OUTLIER_DATA_DIR
#   * /hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/data/
#   * ./untracked_data/data/
#
# Output variables:
#
#   * Path to the outlier data directory

get.outlier.data.dir <- function() {
    data.dirs <- c(
        Sys.getenv('OUTLIER_DATA_DIR'),
        '/hot/project/process/CancerBiology/OUTA-000164-GeneExpressionOABRCA/data/',
        './untracked_data/data/'
        );

    for (data.dir in data.dirs) {
        if (data.dir != '' && dir.exists(data.dir)) {
            return(data.dir);
            }
        }

    stop('No data directory found!');
    }

### save.outlier.figure ############################################################################
# Description:
#
# Save a BPG plot as a PDF and a PNG in a consistent manner.
#
# Input variables:
#
#   plot.object: The trellis plot object to be saved
#   name.segments: A vector of filename segments to be joined with underscores
#   width, height: The dimensions at which to save the plot
#   output.directory: The output directory, defaulting to 'output'. Will be
#     created if necessary

save.outlier.figure <- function(plot.object, name.segments, width, height, output.directory = 'output') {
    if (!dir.exists(output.directory)) {
        dir.create(output.directory);
        }

    save.basename <- file.path(output.directory, paste(name.segments, sep = '_', collapse = '_'));

    # Save a PDF
    grDevices::pdf(
        file = paste0(save.basename, '.pdf'),
        width = width,
        height = height
        );
    print(plot.object);
    grDevices::dev.off();

    # Save a PNG
    grDevices::png(
        file = paste0(save.basename, '.png'),
        width = width,
        height = height,
        unit = 'in',
        res = 1200
        );
    print(plot.object);
    grDevices::dev.off();
    return(invisible(NULL));
    }

### cache.computed.variable #######################################################################
# Description:
#
# Cache a variable for use in a later script. This deliberately only caches a
# single variable at a time. The variable can be re-loaded into the local
# environment with load.computed.variable.
#
# Input variables:
#
#   object.name: The name of the variable to be cached.
cache.computed.variable <- function(object.name) {
    output.directory <- here::here('variable-cache');

    if (!dir.exists(output.directory)) {
        dir.create(output.directory);
        }

    data.file <- file.path(output.directory, paste0(object.name, '.rda'));
    save(list = c(object.name), file = data.file);
    }


cache.multiple.computed.variables <- function(object.names) {
    for (name in object.names) {
        cache.computed.variable(name);
        }
    }

### load.computed.variable #######################################################################
# Description:
#
# Reload a single variable computed in a prior script.
#
# Input variables:
#
#   object.name: The name of the variable to be loaded
load.computed.variable <- function(object.name) {
    data.file <- file.path(here::here('variable-cache'), paste0(object.name, '.rda'));
    load(data.file);
    assign(object.name, get(object.name), envir = parent.frame());
}

load.multiple.computed.variables <- function(object.names) {
    for (name in object.names) {
        load.computed.variable(name);
        assign(name, get(name), envir = parent.frame());
        }
    }
