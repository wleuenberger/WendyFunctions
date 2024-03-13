#' RhatHighest: Grab and plot any terms that haven't converged, or the highest rhat in a category
#'
#' @param parameter
#'
#' @return
#' @export
#'
#' @examples
# Function to grab and plot any terms that haven't converged (rhat >= 1.1) or the highest rhat in a category. Quick glimpse at potential problem terms.
RhatHighest = function(parameter){
  # Create names and indices
  # Add .samples to the end to be able to grab the mcmc chains and every iteration
  samples = paste0(parameter, '.samples')
  # Index the objects with the samples
  samplesindex = which(names(out) == samples)
  # Index the rhat values
  rhatindex = which(names(out$rhat) == parameter)
  # Test if there are terms that haven't converged
  if(any(out$rhat[[rhatindex]] > 1.1)){
    # Keep track
    print('params with rhats >= 1.1')
    # Grab the parameter names and rhat's to display nicely
    Names = dimnames(out[[samplesindex]])[[2]][out$rhat[[rhatindex]] >= 1.1]
    Rhat = out$rhat[[rhatindex]][out$rhat[[rhatindex]] >= 1.1]
    print(tibble(Names, Rhat))
    # Make the traceplot
    traceplot(out[[samplesindex]][,out$rhat[[rhatindex]] >= 1.1])
  } else {
    # If everything has converged, grab the highest rhat to still get a glimpse
    print('beta with highest rhat (< 1.1)')
    # Grab the parameter names and rhat's to display nicely
    Names = dimnames(out[[samplesindex]])[[2]][
      out$rhat[[rhatindex]] == max(out$rhat[[rhatindex]])
    ]
    Rhat = out$rhat[[rhatindex]][out$rhat[[rhatindex]] == max(out$rhat[[rhatindex]])]
    print(tibble(Names, Rhat))
    # Make the traceplot
    traceplot(out[[samplesindex]][,out$rhat[[rhatindex]] == max(out$rhat[[rhatindex]])])
  }
}
