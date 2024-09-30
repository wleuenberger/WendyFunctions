#' ButterflyComparePostHoc
#' Comparing two or more group means to each other
#' Doesn't give all pairwise comparisons for >2 groups
#'
#' @param SpeciesGroups
#' @param Betas
#' @param WhichLists
#' @param NColumns
#'
#' @return
#' @export
#'
#' @examples
ButterflyComparePostHoc <- function(
    SpeciesGroups = SpeciesLists,
    Betas = BetaYears,
    WhichLists,
    NColumns = 136 # Shouldn't need to change, but in case we're not using the full set of species
) {
  NGroups <- length(WhichLists)
  ThisList <- vector('list', length = NGroups)
  SppList <- vector('list', length = NGroups)
  GroupList <- vector('list', length = NGroups)
  for(nn.groups in 1:NGroups) {
    ThisList[[nn.groups]] <- BetaYears[
      names(BetaYears) %in% SpeciesLists[[WhichLists[nn.groups]]]
    ]
    names(ThisList)[nn.groups] <- WhichLists[nn.groups]
    SppList[[nn.groups]] <- SpeciesLists[[WhichLists[nn.groups]]]
    GroupList[[nn.groups]] <- WhichLists[nn.groups] %>%
      rep(times = length(SppList[[nn.groups]]))
  }
  # Note: we don't log these because the response variable is the trends
  # y is not the abundances (trends instead)
  y <- matrix(unlist(ThisList),
              nrow = 1000,
              ncol = NColumns)
  covs <- data.frame(
    spp = 1:NColumns,
    SppName = unlist(SppList),
    Group = unlist(GroupList)
  )
  new.data.list <- list(y = y, covs = covs)
  ComparedResults <- postHocLM(formula = ~ Group,
                               data = new.data.list, n.chains = 1,
                               verbose = TRUE)
  MeanResults <- postHocLM(formula = ~ Group - 1,
                           data = new.data.list, n.chains = 1,
                           verbose = TRUE)

  cat('\n************************************** \n\n',
      WhichLists,
      '\nEffects parameterization
      \n************************************** \n')
  cat(ComparedResults %>% summary)
  cat('\n\n************************************** \n\n',
      WhichLists,
      '\nMeans parameterization
      \n************************************** \n')
  cat(MeanResults %>% summary)

  # Pull out medians/posterior stats
  post.hoc.df <- data.frame(trend.est = NA,
                            trend.prob.pos = NA,
                            trend.025 = NA,
                            trend.25 = NA,
                            trend.75 = NA,
                            trend.975 = NA,
                            Name = WhichLists)
  # Save results and generate summary figure
  # Calculate the probability that the trend over time is positive
  prob.pos <- apply(MeanResults$beta.samples, 2, function(a) mean(a > 0))
  # Save the trend into the post.hoc.df object for later
  # Because we're using the MeanResults object, we remove the [2] from the
  # single-species function
  post.hoc.df$trend.prob.pos <- prob.pos
  # Calculate quantiles of the intercept and slope terms. Median and 95% CI
  # On real scale
  beta.quants <- apply(MeanResults$beta.samples, 2, quantile,
                       probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

  # Fill in the rest of the post.hoc.df summary
  post.hoc.df$trend.est <- beta.quants[3,]
  post.hoc.df$trend.025 <- beta.quants[1,]
  post.hoc.df$trend.25 <- beta.quants[2,]
  post.hoc.df$trend.75 <- beta.quants[4,]
  post.hoc.df$trend.975 <- beta.quants[5,]

  # Add data from the comparison model
  if (length(WhichLists) == 2) {
    beta.quants <- apply(ComparedResults$beta.samples,
                         2,
                         quantile,
                         probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    post.hoc.df %<>%
      add_row(
        Name = paste(WhichLists[1], '-', WhichLists[2]),
        trend.prob.pos = apply(ComparedResults$beta.samples, 2, function(a)
          mean(a > 0))[2],
        trend.est = beta.quants[3, 2],
        trend.025 = beta.quants[1, 2],
        trend.25 = beta.quants[2, 2],
        trend.75 = beta.quants[4, 2],
        trend.975 = beta.quants[5, 2]
      )
    # WhichLists == 2 ends
  } else {
    Listdf <- tibble(One = 1:length(WhichLists),
                     Two = 1:length(WhichLists)) %>%
      expand.grid %>%
      filter(One > Two)
    ListBetas <- MeanResults$beta.samples
    for( n.compare in 1:nrow(Listdf) ){
      Compared <- ListBetas[,Listdf[n.compare, 1]] %>%
        subtract(ListBetas[,Listdf[n.compare, 2]])
      beta.quants <- quantile(Compared,
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
      post.hoc.df %<>%
        add_row(Name = paste(WhichLists[Listdf[n.compare, 1]],
                             '-',
                             WhichLists[Listdf[n.compare, 2]]),
                trend.prob.pos = mean(Compared > 0),
                trend.est = beta.quants[3],
                trend.025 = beta.quants[1],
                trend.25 = beta.quants[2],
                trend.75 = beta.quants[4],
                trend.975 = beta.quants[5]
        )
    } # n.compare
  } # WhichLists != 2


  ReturnList <- list('Data' = new.data.list,
                     'Compared' = ComparedResults,
                     'Mean' = MeanResults,
                     'SummaryStats' = post.hoc.df)


  return(ReturnList)
}
