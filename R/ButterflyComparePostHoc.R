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
  ReturnList <- list('Data' = new.data.list,
                     'Compared' = ComparedResults,
                     'Mean' = MeanResults)
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
  return(ReturnList)
}
