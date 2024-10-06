#' ButterflyRichness calculation
#'
#' @param RichnessData
#' @param WhichLists
#'
#' @return
#' @export
#'
#' @examples
ButterflyRichness <- function (
    RichnessData = RichnessList,
    WhichLists
    ) {
  # Preparation code
  # Alphabetize which list to decrease labelling problems
  # More sophisticated methods exist to tie output to names, but this works
  # for my data
  WhichLists <- WhichLists %>% sort
  # Figure out how many groups to loop over and make a list of that size
  NGroups <- WhichLists %>% length
  ThisList <- vector('list', length = NGroups)

  # Pull out the specific parts of RichnessList
  for ( nn.groups in 1:NGroups ) {
    ThisList[[nn.groups]] <- RichnessData[[WhichLists[nn.groups]]]
    names(ThisList)[nn.groups] <- WhichLists[nn.groups]
  }

  # Dimensions
  # n.iter x n.site x n.year x n.groups
  # Sites (counties)
  J <- ncol(ThisList[[1]])
  # Number of iterations
  n.samples <- nrow(ThisList[[1]])
  # Number of years
  unique.years <- 1:dim(ThisList[[1]])[3]
  # Numbers for each year, duplicated for the number of groups
  years <- rep(1:length(unique.years), each = J)
  years %<>% rep(times = NGroups)
  # Number of years for duplicating other info
  n.years <- length(unique(unique.years))
  # Create site numbers for each year
  sites <- rep(1:J, times = n.years)
  # Duplicate for the number of groups
  sites %<>% rep(times = NGroups)
  # Duplicate group for each site/year
  groups <- rep(WhichLists, each = J * n.years)

  # Format data for postHocLM
  # Make into matrix instead of list with particular dimensions
  y <- matrix(unlist(ThisList),
              nrow = n.samples, ncol = J * n.years * NGroups)
  # Some richness values are 0, so we can't do log(y) (log[0] = -Inf)
  new.data.list <- list(y = y,
                        covs = data.frame(years, sites, groups))

  # Run model
  out.posthoc <- postHocLM(formula = ~ years * groups +
                             (1 | sites),
                           data = new.data.list,
                           n.chains = 1,
                           verbose = TRUE)

  # Post-processing
  # Pull out medians/posterior stats on trend
  # Maybe add intercept later?

  # Mean trends
  # Create vector to hold means
  # Useful for comparisons
  MeanTrend <- vector('list', length = NGroups)

  for ( nn.groups in 1:NGroups ) {
    # 'Years' is the reference level, start the data frame with those quantiles
    if( nn.groups == 1 ) {
      MeanTrend[[nn.groups]] <- out.posthoc$beta.samples[ , 'years']
      post.hoc.df <- data.frame(
        # Median
        trend.est = quantile(out.posthoc$beta.samples[, 'years'], 0.5) %>% unname,
        # Proportion of posterior above 0
        trend.prob.pos = mean(out.posthoc$beta.samples[, 'years'] > 0),
        # Quantiles
        trend.025 = quantile(out.posthoc$beta.samples[, 'years'], 0.025) %>% unname,
        trend.25 = quantile(out.posthoc$beta.samples[, 'years'], 0.25) %>% unname,
        trend.75 = quantile(out.posthoc$beta.samples[, 'years'], 0.75) %>% unname,
        trend.975 = quantile(out.posthoc$beta.samples[, 'years'], 0.975) %>% unname,
        # Which group
        Name = WhichLists[nn.groups]
      )
    } else { # nn.groups > 1
      # Create the name (customized to this model)
      CompareName <- paste0('years:groups', WhichLists[nn.groups])
      # Effects parameterization - add the intercept + group effect together
      Added <- out.posthoc$beta.samples[ , CompareName] +
        out.posthoc$beta.samples[, 'years']
      # Add to the mean trend list for comparisons
      MeanTrend[[nn.groups]] <- Added

      # Fill in the dataframe, same as above
      post.hoc.df %<>%
        add_row( Name = WhichLists[nn.groups],
                 trend.prob.pos = mean(Added > 0),
                 trend.est = quantile(Added, 0.5) %>% unname,
                 trend.025 = quantile(Added, 0.025) %>% unname,
                 trend.25 = quantile(Added, 0.25) %>% unname,
                 trend.75 = quantile(Added, 0.75) %>% unname,
                 trend.975 = quantile(Added, 0.975) %>% unname)
    }
  }

  # Run comparisons
  # Determine how many comparisons and among all levels
  Listdf <- tibble(One = 1:length(WhichLists),
                   Two = 1:length(WhichLists)) %>%
    expand.grid %>%
    filter(One > Two)

  # Loop through comparisons
  for ( nn.compare in 1:nrow(Listdf) ) {
    Compared <- MeanTrend[[ Listdf[ nn.compare, 1 ] ]] %>%
      subtract(MeanTrend[[ Listdf[ nn.compare, 2 ] ]])
    post.hoc.df %<>%
      add_row( Name = paste(WhichLists[ Listdf[ nn.compare, 1 ] ],
                            '-',
                            WhichLists[ Listdf[ nn.compare, 2 ] ]),
               trend.prob.pos = mean(Compared > 0),
               trend.est = quantile(Compared, 0.5) %>% unname,
               trend.025 = quantile(Compared, 0.025) %>% unname,
               trend.25 = quantile(Compared, 0.25) %>% unname,
               trend.75 = quantile(Compared, 0.75) %>% unname,
               trend.975 = quantile(Compared, 0.975) %>% unname)

  }

  ReturnList <- list('Data' = new.data.list,
                     'Model' = out.posthoc,
                     'SummaryStats' = post.hoc.df)
  return(ReturnList)


}
