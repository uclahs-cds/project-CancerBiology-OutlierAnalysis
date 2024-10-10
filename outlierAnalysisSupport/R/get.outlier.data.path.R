#' Return the path to the outlier data file
#'
#' @param data.dir Outlier data directory
#' @param data.file Outlier data filename
#'
#' @return Full path to the outlier data file, suitable for 'attach'
#' @export
get.outlier.data.path <- function(data.dir = Sys.getenv('OUTLIER_DATA_DIR'),
                                  data.file = '2024-10-08_Figure1_2_3_4_min_input.rda') {
    full.data.file <- file.path(data.dir, data.file);

    # Check if file exists before attaching
    if (!file.exists(full.data.file)) {
        stop(paste('The file', full.data.file, 'does not exist.'));
        }

    return(full.data.file);
}