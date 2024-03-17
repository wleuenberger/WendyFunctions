#' Calculate figure sizes for the zipkinlab front page
#'
#' @param OriginalWidth: Original height of the figure
#' @param OriginalHeight: Original width of the figure
#' @param OutputWidth: Desired width of the output figure. Set to 200 px as default, but can be changed
#'
#' @return: Width: Desired width. Currently set to 200 as default but can be changed
#' Height: Desired height based on the desired width and original dimensions
#' @export
#'
#' @examples
#'
#'
FigSize <- function(OriginalWidth, OriginalHeight, OutputWidth = 200){
  WidthRatio <- OutputWidth / OriginalWidth
  OutputHeight <- OriginalHeight * WidthRatio
  cat(paste0('Width = ', round(OutputWidth), '\nHeight = ', round(OutputHeight)))
}
