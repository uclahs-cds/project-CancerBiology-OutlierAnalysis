### UTILITY FUNCTIONS ##############################################################################

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
    pdf(
        file = paste0(save.basename, '.pdf'),
        width = width,
        height = height
        );
    print(plot.object);
    dev.off();

    # Save a PNG
    png(
        file = paste0(save.basename, '.png'),
        width = width,
        height = height,
        unit = 'in',
        res = 1200
        );
    print(plot.object);
    dev.off();
    }

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
