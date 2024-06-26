#' RhatHighest: Grab and plot any terms that haven't converged, or the highest rhat in a category
#'
#' @param parameter
#'
#' @return
#' @export
#'
#' @examples
# Function to grab and plot any terms that haven't converged (rhat >= 1.1) or the highest rhat in a category. Quick glimpse at potential problem terms.
RhatHighest = function(parameter, Colors = 1:out$n.chains) {
  # Create names and indices
  # Add .samples to the end to be able to grab the mcmc chains and every iteration
  samples = paste0(parameter, '.samples')
  # Index the objects with the samples
  samplesindex = which(names(out) == samples)
  # Index the rhat values
  rhatindex = which(names(out$rhat) == parameter)
  # Test if there are terms that haven't converged
  if (any(out$rhat[[rhatindex]] > 1.1)) {
    # Keep track
    print('params with rhats >= 1.1')
    # Grab the parameter names and rhat's to display nicely
    Names = dimnames(out[[samplesindex]])[[2]][out$rhat[[rhatindex]] >= 1.1]
    Rhat = out$rhat[[rhatindex]][out$rhat[[rhatindex]] >= 1.1]
    print(tibble(Names, Rhat))
    # Make the traceplot
    # traceplot(as.mcmc(out[[samplesindex]][, out$rhat[[rhatindex]] >= 1.1]))

    for (rr in 1:length(Names)) {
      # Pull out the columns of data with rhats > 1.1
      Above1.1 <- out[[samplesindex]][ , out$rhat[[rhatindex]] > 1.1]
      if(length(Names) != 1){
        # Use different indexing if there's multiple terms with rhats > 1.1
        Above1.1 <- out[[samplesindex]][ , out$rhat[[rhatindex]] > 1.1][,rr]
      }
      for (cc in 1:out$n.chains) {
        iters <- (((cc - 1) * out$n.post) + 1):(cc * out$n.post)
        # Chain <- as.mcmc(out[[samplesindex]][iters, out$rhat[[rhatindex]] == max(out$rhat[[rhatindex]])])
        if (cc == 1) {
          traceplot(
            as.mcmc(Above1.1[iters]),
            main = paste0(Names[rr], ' \n Rhat = ', Rhat[rr] %>% round(3)))
        } # cc == 1
        if (cc != 1) {
          lines(Above1.1[iters],
                col = Colors[cc])
        } # cc != 1

      } # cc
    } # rr

  } else {

    # If everything has converged, grab the highest rhat to still get a glimpse
    print('beta with highest rhat (< 1.1)')
    # Grab the parameter names and rhat's to display nicely
    Names = dimnames(out[[samplesindex]])[[2]][out$rhat[[rhatindex]] == max(out$rhat[[rhatindex]])]
    Rhat = out$rhat[[rhatindex]][out$rhat[[rhatindex]] == max(out$rhat[[rhatindex]])]
    print(tibble(Names, Rhat))
    # Make the traceplot
    # Chains concatenated one after another
    # traceplot(out[[samplesindex]][,out$rhat[[rhatindex]] == max(out$rhat[[rhatindex]])])
    # Chains on top of one another
    for (cc in 1:out$n.chains) {
      iters <- (((cc - 1) * out$n.post) + 1):(cc * out$n.post)
      Chain <- as.mcmc(out[[samplesindex]][iters, out$rhat[[rhatindex]] == max(out$rhat[[rhatindex]])])
      if (cc == 1) {
        traceplot(Chain, main = paste0(Names, '; \nRhat = ', Rhat %>% round(3)))
      }
      if (cc != 1) {
        lines(Chain, col = Colors[cc])
      }
    }

  }
}
