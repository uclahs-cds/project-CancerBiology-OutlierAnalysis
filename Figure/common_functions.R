
save.outlier.figure <- function(plot.object, name.segments, width, height, depth = 0) {
    # Find or create an output folder with the current PID and datestamp
    output.dir <- paste('figures', format(Sys.time(), '%Y%m%dT%H%M%S'), Sys.getpid(), sep = '-');
    existing.dirs <- Sys.glob(paste0('figures-*-', Sys.getpid()));

    if (length(existing.dirs) >= 2) {
        print(existing.dirs);
        stop('Multiple output directories detected');
    } else if (length(existing.dirs) == 1) {
        output.dir <- existing.dirs[1];
    } else {
        dir.create(output.dir);
    }

    # Get the basename of the file calling this function
    sourcing.filename <- tools::file_path_sans_ext(basename(parent.frame(3 + depth)$ofile));

    save.basename <- file.path(output.dir, paste(c(sourcing.filename, name.segments), sep = '_', collapse = '_'));

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
