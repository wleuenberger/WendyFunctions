#' Run a post hoc regression on my butterfly response data
#'
#' @param Data: Data in the format iteration x county x year, 1000 x 262 x 32
#' @param Logged: Whether or not to log transform the response variable
#'
#' @return
#' @export
#'
#' @examples
ButterflyPostHoc <- function( Data, Logged = TRUE) {
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
  sites <- rep(1:J, times = n.years)
  # Format data for postHocLM
  y <- matrix(Data, nrow = n.samples, ncol = J * n.years)
  covs <- data.frame(years, sites)

  if (Logged == TRUE) {
    new.data.list <- list(y = log(y),
                          covs = covs)
  } else {
    new.data.list <- list(y = y,
                          covs = covs)
  }

  out.posthoc <- postHocLM(
    formula = ~ years + (1 | sites),
    data = new.data.list,
    n.chains = 1,
    verbose = TRUE
  )
  return(list(ModelOutput = out.posthoc,
              BetaYears = out.posthoc$beta.samples[,"years"]))
} # ButterflyPostHoc function end
