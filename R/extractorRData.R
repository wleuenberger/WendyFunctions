# From https://stackoverflow.com/questions/65964064/programmatically-extract-an-object-from-collection-of-rdata-files
#' Extract data from an .RData file
#'
#' @param file
#' @param object
#'
#' @return
#' @export
#'
#' @examples
extractorRData <- function(file, object) {
  #' Function for extracting an object from a .RData file created by R's save() command
  #' Inputs: RData file, object name
  E <- new.env()
  load(file=file, envir=E)
  return(get(object, envir=E, inherits=F))
}
