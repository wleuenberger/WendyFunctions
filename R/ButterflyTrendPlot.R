# Deprecated!
#' Use output from ButterflyPostHoc to organize trend, plot data, and return values
#'
#' @param Data: Butterfly relative abundance data in the formating n.iter (1000) x n.counties x n.years
#' @param Model: Model using the above data
#' @param Logged: Whether or not to log-transform the response variable
#' @param ylabel: Ylabel on resulting figure
#' @param Name: Name on resulting figure
#' @param SaveFig: Whether or not to save the results figure
#'
#' @return
#' @export
#'
#' @examples
ButterflyTrendPlot <- function (
    Data,
    Model,
    Logged = TRUE,
    Name = 'ACHLYC',
    ylabel = 'Relative Abundance Index (log scale)',
    RunFigCode = FALSE) {
  post.hoc.df <- data.frame(trend.est = NA,
                            int.est = NA,
                            int.est.median = NA,
                            trend.prob.pos = NA,
                            trend.025 = NA,
                            trend.25 = NA,
                            trend.75 = NA,
                            trend.975 = NA,
                            Name = Name)
  # Save results and generate summary figure
  # Calculate the probability that the trend over time is positive
  prob.pos <- apply(Model$beta.samples, 2, function(a) mean(a > 0))
  # Save the trend into the post.hoc.df object for later
  post.hoc.df$trend.prob.pos <- prob.pos[2]
  # Calculate quantiles - median (0.5), as well as 2.5 and 97.5
  if (Logged == TRUE) {
    # These are quantiles of the Data (which was model output from original),
    # Not from the model output of the mean. Therefore, they are much wider than
    # the mean trend
    sum.indx.log.quants <- apply(log(Data), 3, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    # Means of the data for comparison
    means <- apply(log(Data), 2, mean)
    # Calculate quantiles of the intercept and slope terms. Median and 95% CI
    # On real scale
    beta.quants <- apply(exp(Model$beta.samples), 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    # On log scale
    beta.log.quants <- apply(Model$beta.samples, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  } else {
    # These are quantiles of the Data (which was model output from original),
    # Not from the model output of the mean. Therefore, they are much wider than
    # the mean trend
    # Keep the log name for simplicity later
    sum.indx.log.quants <- apply(Data, 3, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    # Means of the data
    means <- apply(Data, 2, mean)
    # Calculate quantiles of the intercept and slope terms. Median and 95% CI
    # On real scale
    beta.quants <- apply(Model$beta.samples, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    # Also on the real scale, but keep the log name
    beta.log.quants <- apply(Model$beta.samples, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  } # if else logged

  # Fill in the rest of the post.hoc.df summary
  post.hoc.df$trend.est <- beta.log.quants[3, 2]
  post.hoc.df$int.est <- beta.log.quants[3, 1]
  post.hoc.df$trend.025 <- beta.log.quants[1, 2]
  post.hoc.df$trend.25 <- beta.log.quants[2, 2]
  post.hoc.df$trend.75 <- beta.log.quants[4, 2]
  post.hoc.df$trend.975 <- beta.log.quants[5, 2]

  unique.years <- unique(Model$X[, 'years'])
  plot.df <- data.frame(med = sum.indx.log.quants[3, ],
                        means = means[1],
                        low = sum.indx.log.quants[1, ],
                        twentyfive = sum.indx.log.quants[2, ],
                        seventyfive = sum.indx.log.quants[4, ],
                        high = sum.indx.log.quants[5, ],
                        year = unique.years)

  if (RunFigCode == TRUE) {
    # Graphing output
    # SLOPE is fine and the same either way
    # INTERCEPT is challenging because medians are lower than means
    # Median points requires median(Intercept + Ranef) then Median again (for one point)
    #   If we don't add the ranef, then the line won't match the plots
    # Mean points do not require that transformation, can just use the intercept from the mode
    # Happening because of huge butterfly counts/overdispersion - Means are much higher than median
    # intercept needs to be adjusted by ranef
    # Intercept samples
    # (1 per iteration)
    int.samples <- Model$beta.samples[, 1]
    # Add the random effect samples for each site to the intercept samples
    # (1 per iteration x county)
    full.int.samples <- Model$beta.star.samples + int.samples
    # Take the median values of the ranef + intercept values
    # (1 per iteration, median across the rows of counties)
    overall.mean.samples <- apply(full.int.samples, 1, median)
    # Take the median overall to get a point estimate
    # (1 value)
    post.hoc.df$int.est.median <- median(overall.mean.samples)
    # For figure - median beta.star
    beta.star.median <- apply(Model$beta.star.samples, 1, median)

    unique.years <- unique(Model$X[, 'years'])
    plot.min <- min(sum.indx.log.quants)
    plot.max <- max(sum.indx.log.quants)
    plot.vals <- seq(
      from = min(unique.years),
      to = max(unique.years),
      length.out = 100
    )

    # Generate min and max y values
    for (ss in 1:999) {
      tmp <- Model$beta.samples[ss, 1] + Model$beta.samples[ss, 2] * plot.vals
      if (min(tmp) < plot.min) {
        plot.min <- min(tmp)
      } # min
      if (max(tmp) > plot.max) {
        plot.max <- max(tmp)
      } # max
    } # ss samples

    # Lines <- data.frame(year = plot.vals,
    #                     med = tmp)
    # sp.name.plot <- simpleCap(str_replace(curr.sp.name, '-', ' '))

    # ggplot structure
    my.plot <- ggplot(data = plot.df, aes(x = year, y = med))
    # 1 line per sample (takes awhile)
    for (ss in sample(1:999, 200)) {
      my.plot <- my.plot + geom_abline(
        slope = Model$beta.samples[ss, 2],
        intercept = Model$beta.samples[ss, 1] - beta.star.median,
        color = 'darkorchid4',
        alpha = 0.1
      )
    }
    # Add points, error bars, axis labels
    my.plot <- my.plot +
      geom_point(size = 3) +
      # geom_segment(aes(x = year, y = low, xend = year, yend = high)) +
      theme_classic(base_size = 18) +
      labs(x = 'Year', y = ylabel)
    # Add stats and species name
    my.plot <- my.plot +
      ggtitle(
        label = paste(
          Name,
          ": P(trend < 0) = ",
          round(1 - prob.pos[2], 2),
          "\n% Change/Year: ",
          round((beta.quants[3, 2] - 1) * 100, 2),
          " (",
          round((beta.quants[1, 2] - 1) * 100, 2),
          ", ",
          round((beta.quants[5, 2] - 1) * 100, 2),
          ")",
          sep = ''
        )
      )
    # Add trendline
    my.plot <- my.plot +
      geom_abline(
        intercept = post.hoc.df$int.est.median,
        slope = post.hoc.df$trend.est,
        linewidth = 1
      )
    # lines(plot.vals, tmp, col = 'black', lty = 1, lwd = 5)

    if (SaveFig == TRUE) {
      print(my.plot)
      pdf(
        paste('Output/trend-figures/', Name, '-trend-species.pdf', sep = ''),
        width = 5.87,
        height = 7.52
      )
      print(my.plot)
      dev.off()
    } else {
      # print(my.plot)
    }# save figure


    # save(out.posthoc, file = paste('Output/trend-figures/', ID, 'species-',
    # species[ii], '-posthoc-for-figure.rda', sep = ''))

  } # RunFigCode

  if(RunFigCode == FALSE) {
    OutputList <- list(post.hoc.df, plot.df)
  } else {
    OutputList <- list(post.hoc.df, plot.df, my.plot)
  }

  return(OutputList)

} # ButterflyTrendPlot
