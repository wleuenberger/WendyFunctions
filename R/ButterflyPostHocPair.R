#' Paired (or multiple covariates) test of butterfly post hoc data
#'
#' @param Data: List of species names in different groups
#' @param Names: Named elements in the Data list that will be tested
#' @param Files: Where the species-specific files are found
#'
#' @return: Model (out.posthoc, model results) and the combined data set (CombinedData)
#' @export
#'
#' @examples
ButterflyPostHocPair <- function(
    Data,
    Names,
    Files) {
  # State which groups we're working on
  print( Names)

  DataGroups <- tibble(Trait = character(), Species = character())
  for( nn.traits in 1:length(Names) ) {
    DataPartialGroup <- Data[names(Data) %in% Names[nn.traits]] %>% as_tibble %>%
      pivot_longer(cols = everything(),
                   names_to = 'Trait',
                   values_to = 'Species')
    DataGroups %<>% bind_rows(DataPartialGroup)
  } # nn.traits


  for( nn.traits in 1:length(Names)) {
    Composite.indx <- array(data = 0, dim = c(1000, 262, 32))
    Spp <- DataGroups %>%
      filter(Trait == Names[nn.traits])
    # Grab files
    CompositeFiles <- Files %>%
      filter(SpeciesNames %in% Spp$Species)
    Spp$Species[!Spp$Species %in% CompositeFiles$SpeciesNames]
    for ( nn.species in 1:dim(Spp)[1] ) {
      # print(CompositeFiles$SpeciesNames[nn.species])
      Species.indx <- extractorRData(CompositeFiles$FullPath[nn.species],
                                     object = 'Species.indx')
      Composite.indx <- Composite.indx + Species.indx
    } # nn.species
    # Put the trait info with the species info
    CompositeData <- ButterflyPostHocPrep(Data = Composite.indx,
                                          Logged = TRUE)
    CompositeData$covs$Trait <- Names[nn.traits]
    if( nn.traits == 1 ) {
      print('Create CompositeTraitReady')
      CompositeTraitReady <- CompositeData
    } else {
      print('Combine new composite with other traits composite')
      CompositeTraitReady$y <- cbind(CompositeTraitReady$y, CompositeData$y)
      CompositeTraitReady$covs <- rbind(CompositeTraitReady$covs, CompositeData$covs)
    } # if 1 else loop end
  } # nn.trait
  out.posthoc <- postHocLM(formula = ~ years * Trait  + (1 | sites),
                           data = CompositeTraitReady,
                           n.chains = 1,
                           verbose = TRUE)
  return(list(Model = out.posthoc, CombinedData = CompositeTraitReady))

}
