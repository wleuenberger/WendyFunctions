#' Run a post hoc regression on my butterfly response data
#'
#' @param Data: Data in the format iteration x county x year, 1000 x 262 x 32
#' @param Logged: Whether or not to log transform the response variable
#'
#' @return
#' @export
#'
#' @examples
ButterflyPostHoc <- function( Data, Logged = TRUE, YearSplit = FALSE) {
  if( any(dim(Data) == c(1000, 262, 32)) == FALSE ) {
    stop(
      paste0('Different dimensions than expected. Expected c(1000, 262, 32).
            Data had dimensions of ',
             dim(Data)[1], ', ',
             dim(Data)[2], ', ',
             dim(Data)[3])
    )
  } # Dimension check

  # Reformat data into a matrix
  # Get the dimensions and years
  J <- ncol(Data)
  n.samples <- nrow(Data)
  n.years <- dim(Data)[3]
  unique.years <- 1:dim(Data)[3]
  plot.years <- unique.years
  # Repeat data to change into 2D
  years <- rep(1:n.years, each = J)
  # Put a breakpoint halfway through
  if(YearSplit == TRUE) {
    Yearsplit <- tibble(years_raw = years)
    Yearsplit  <- Yearsplit
      mutate(yearsplit = case_when(
        years_raw <= n.years/2 ~ paste0('1:', n.years/2),
        years_raw > n.years/2 ~ paste0((n.years/2) + 1, ':', n.years)
      ),
      years = ifelse(years_raw <= n.years/2, years_raw, years_raw - n.years/2)
      )
  } # YearSplit == TRUE
  sites <- rep(1:J, times = n.years)
  # Format data for postHocLM
  y <- matrix(Data, nrow = n.samples, ncol = J * n.years)
  if(YearSplit == TRUE){
    covs <- data.frame(years = Yearsplit$years,
                       sites,
                       yearsplit = Yearsplit$yearsplit)
  } else {
    covs <- data.frame(years, sites)
  }

  if (Logged == TRUE) {
    new.data.list <- list(y = log(y),
                          covs = covs)
  } else {
    new.data.list <- list(y = y,
                          covs = covs)
  }

  if( YearSplit == TRUE) {
    out.posthoc <- postHocLM(formula = ~ years * yearsplit + (1 | sites),
                             data = new.data.list,
                             n.chains = 1,
                             verbose = TRUE)
  } else {
    out.posthoc <- postHocLM(formula = ~ years + (1 | sites),
                               data = new.data.list,
                               n.chains = 1,
                               verbose = TRUE)
  }
  return(out.posthoc)
} # ButterflyPostHoc function end
