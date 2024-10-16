#' Download and cache a file locally.
#'
#' If the file extension ends with `.gz`, it will automatically be uncompressed.
#'
#' @param url The full URL to the file to be cached
#'
#' @return Path to the cached uncompressed file.
#' @export

cache.file <- function(url) {
    cache.directory <- here::here('Figure', 'output', 'download-cache');

    if (!dir.exists(cache.directory)) {
        dir.create(cache.directory, recursive = TRUE);
        }

    file.name <- basename(url);
    cached.file <- file.path(cache.directory, file.name);

    if (!file.exists(cached.file)) {
        download.file(url, destfile = cached.file, mode = 'wb');
        }

    # Check if the file is gzipped (based on extension)
    if (grepl('\\.gz$', cached.file)) {
        decompressed.file <- sub('\\.gz$', '', cached.file);

        if (!file.exists(decompressed.file)) {
            # Open the gzipped file for reading in binary mode
            gz.con <- gzfile(cached.file, 'rb')

            # Read the content
            content <- readBin(gz.con, what = 'raw', n = 1e6)
            close(gz.con)

            # Write the decompressed content to a new file
            writeBin(content, decompressed.file)
            }

        return(decompressed.file);
        }

    return(cached.file);
    }
