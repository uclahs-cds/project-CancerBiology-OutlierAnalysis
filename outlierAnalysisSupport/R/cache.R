#' Cache a variable for later use.
#' 
#' Cache a variable for use in a later script. This deliberately only caches a
#' single variable at a time. The variable can be re-loaded into the local
#' environment with load.computed.variable.
#'
#' @param object.name The name of the variable to cache.
#'
#' @return `NULL`.
#' @export
cache.computed.variable <- function(object.name) {
  output.directory <- here::here('variable-cache');
  
  if (!dir.exists(output.directory)) {
    dir.create(output.directory);
  }
  
  data.file <- file.path(output.directory, paste0(object.name, '.rda'));
  save(list = c(object.name), file = data.file);
  return(invisible(NULL));
}

#' Cache multiple variables for later use.
#'
#' @param object.names A vector of variable names to cache
#'
#' @return `NULL`.
#' @export
cache.multiple.computed.variables <- function(object.names) {
  for (name in object.names) {
    cache.computed.variable(name);
  }
  return(invisible(NULL));
}

### load.computed.variable #######################################################################
# Description:
#
# Reload a single variable computed in a prior script.
#
# Input variables:
#
#   object.name: The name of the variable to be loaded
#' Reload a cache variable.
#' 
#' This will load the named variable into the global namespace.
#'
#' @param object.name The name of the variable to be loaded.
#'
#' @return `NULL`.
#' @export
load.computed.variable <- function(object.name) {
  data.file <- file.path(here::here('variable-cache'), paste0(object.name, '.rda'));
  load(data.file);
  assign(object.name, get(object.name), envir = parent.frame());
  return(invisible(NULL));
}

#' Reload multiple variables computed in a prior script.
#'
#' @param object.names A vector of variable names to be loaded.
#'
#' @return `NULL`.
#' @export
load.multiple.computed.variables <- function(object.names) {
  for (name in object.names) {
    load.computed.variable(name);
    assign(name, get(name), envir = parent.frame());
    }
  return(invisible(NULL));
}