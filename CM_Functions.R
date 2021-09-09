library(purrr)

##############################################################################################
# Make functions which calculates 95th percentile of expected mistags per species
# AND 95th percentile of mistags per dual-index bin

# one function assuming uniform partitioning of mistags across dual-index bins
MaxMistagsPerBinU <- function(Samples, SpeiciesAndAbundances, MistagRate, iters) {
  SpeiciesAndAbundances <- SpeiciesAndAbundances[order(SpeiciesAndAbundances[,2]),]
  SeqsPerSpeciesList <- SpeiciesAndAbundances[,2]
  Taxa <- SpeiciesAndAbundances[,1]
  MistagLst <- qbinom(0.05, SeqsPerSpeciesList, MistagRate, lower.tail=F)
  MtsPerBinU <- c() 
  MinMtsPerBin <- 0
  for (a in 1:length(MistagLst)) {
    MaxValsU <- c()
    for (b in 1:iters) {
      BinLstU <- rep(0, Samples)
      du <- rdunif(MistagLst[a],Samples) }
      for (c in du) { 
        BinLstU[c] <- BinLstU[c] + 1 }
      MaxValsU[b] <- max(BinLstU, na.rm = T) 
    MtsPerBinU[a] <- quantile(MaxValsU,0.95, na.rm = T)
    if (MtsPerBinU[a] < MinMtsPerBin) {
      MtsPerBinU[a] <- MinMtsPerBin }
    MinMtsPerBin <- MtsPerBinU[a]
  } 
  df <- data.frame("Taxa"=Taxa,"SeqsPerSpecies"=SeqsPerSpeciesList,"Mistags95thP"=MistagLst,"MaxMTsPerBin_Uniform"=MtsPerBinU)
  return(df)
}

######################################################################################
# another function assuming clustered partitioning of mistags across dual-index bins
# severity of clustering modeled as X ~ Pois(Lambda) 
# Appropriate lambda estimated from distribution of mistags per control bin

MaxMistagsPerBinP <- function(Samples, SpeiciesAndAbundances, MistagRate, iters, lambda) {
  SpeiciesAndAbundances <- SpeiciesAndAbundances[order(SpeiciesAndAbundances[,2]),]
  SeqsPerSpeciesList <- SpeiciesAndAbundances[,2]
  Taxa <- SpeiciesAndAbundances[,1]
  MistagLst <- qbinom(0.05, SeqsPerSpeciesList, MistagRate, lower.tail=F)
  MtsPerBinP <- c() 
  MinMtsPerBin <- 0
  for (a in 1:length(MistagLst)) {
    MaxValsP <- c()
    for (b in 1:iters) {
      BinLstP <- rep(0, Samples)
      dp <- rpois(MistagLst[a],lambda = lambda) }
    for (c in dp) { 
      BinLstP[c] <- BinLstP[c] + 1 }
    MaxValsP[b] <- max(BinLstP, na.rm = T) 
    MtsPerBinP[a] <- quantile(MaxValsP,0.95, na.rm = T)
    if (MtsPerBinP[a] < MinMtsPerBin) {
      MtsPerBinP[a] <- MinMtsPerBin }
    MinMtsPerBin <- MtsPerBinP[a]
  } 
  df <- data.frame("Taxa"=Taxa,"SeqsPerSpecies"=SeqsPerSpeciesList,"Mistags95thP"=MistagLst,"MaxMTsPerBin_Poisson"=MtsPerBinP)
  return(df)
}

