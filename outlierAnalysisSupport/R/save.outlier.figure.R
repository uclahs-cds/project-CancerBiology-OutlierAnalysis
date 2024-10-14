#' Save a figure in standard formats.
#'
#' @param plot.object The trellis plot object to be saved
#' @param name.segments A vector of filename segments to be joined with underscores
#' @param width The output image's width
#' @param height The output image's height
#' @param output.directory The output directory, defaulting to 'output'. Will be created if necessary.
#'
#' @return The path to the PNG image.
#' @export
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
  return(paste0(save.basename, 'png'));
}
