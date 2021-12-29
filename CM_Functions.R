library(purrr)

##############################################################################################
# Make functions which calculates 95th percentile of expected mistags per species
# AND 95th percentile of mistags per dual-index bin

# one function assuming uniform partitioning of mistags across dual-index bins
MaxMistagsPerBinU <- function(Samples, SpeiciesAndAbundances, MistagRate, Iters, FDR) {
  FDR <- as.numeric(as.character(FDR))
  UpperPercentile <- (1 - FDR) 
  SpeiciesAndAbundances <- SpeiciesAndAbundances[order(SpeiciesAndAbundances[,2]),]
  SeqsPerSpeciesList <- SpeiciesAndAbundances[,2]
  Taxa <- SpeiciesAndAbundances[,1]
  MistagLst <- qbinom(FDR, SeqsPerSpeciesList, MistagRate, lower.tail=F)
  MtsPerBinU <- c() 
  MinMtsPerBin <- 0
  for (a in 1:length(MistagLst)) {
    MaxValsU <- c()
    for (b in 1:Iters) {
      set.seed(b)
      BinLstU <- rep(0, Samples)
      du <- rdunif(MistagLst[a],Samples) }
      for (c in du) { 
        BinLstU[c] <- BinLstU[c] + 1 }
      MaxValsU[b] <- max(BinLstU, na.rm = T) 
    MtsPerBinU[a] <- quantile(MaxValsU,UpperPercentile, na.rm = T)
    if (MtsPerBinU[a] < MinMtsPerBin) {
      MtsPerBinU[a] <- MinMtsPerBin }
    MinMtsPerBin <- MtsPerBinU[a]
  } 
  df <- data.frame("Taxa"=Taxa,"SeqsPerSpecies"=SeqsPerSpeciesList,"MistagsAtFDR"=MistagLst,"MaxMTsPerBin_Uniform"=MtsPerBinU)
  rownames(df) <- df$Taxa
  return(df)
}

######################################################################################
# another function assuming clustered partitioning of mistags across dual-index bins
# severity of clustering modeled as X ~ Pois(Lambda) 
# Appropriate lambda estimated from distribution of mistags per control bin

MaxMistagsPerBinP <- function(Samples, SpeiciesAndAbundances, MistagRate, Iters, Lambda, FDR) {
  FDR <- as.numeric(as.character(FDR))
  UpperPercentile <- (1 - FDR)
  SpeiciesAndAbundances <- SpeiciesAndAbundances[order(SpeiciesAndAbundances[,2]),]
  SeqsPerSpeciesList <- SpeiciesAndAbundances[,2]
  Taxa <- SpeiciesAndAbundances[,1]
  MistagLst <- qbinom(FDR, SeqsPerSpeciesList, MistagRate, lower.tail=F)
  MtsPerBinP <- c() 
  MinMtsPerBin <- 0
  for (a in 1:length(MistagLst)) {
    MaxValsP <- c()
    for (b in 1:Iters) {
      BinLstP <- rep(0, Samples)
      dp <- rpois(MistagLst[a],lambda = Lambda) }
    for (c in dp) { 
      BinLstP[c] <- BinLstP[c] + 1 }
    MaxValsP[b] <- max(BinLstP, na.rm = T) 
    MtsPerBinP[a] <- quantile(MaxValsP,UpperPercentile, na.rm = T)
    if (MtsPerBinP[a] < MinMtsPerBin) {
      MtsPerBinP[a] <- MinMtsPerBin }
    MinMtsPerBin <- MtsPerBinP[a]
  } 
  df <- data.frame("Taxa"=Taxa,"SeqsPerSpecies"=SeqsPerSpeciesList,"MistagsAtFDR"=MistagLst,"MaxMTsPerBin_Poisson"=MtsPerBinP)
  rownames(df) <- df$Taxa
  return(df)
}

#####################################################################################################
# Make a function which filters your data using MistagThresholds calculated using MaxMistagsPerBinU()
# Detections input table must have ################# 
FilterDetections <- function(MistagThresholds, Detections){
  Out <- as.data.frame(matrix(nrow=nrow(Detections),ncol = ncol(Detections)))
  for (row in 1:nrow(Out)) {
    for (col in 1:ncol(Out)){
      Fam <- rownames(Detections)[row] 
      Out[row,col] <- ifelse(as.numeric(as.character(Detections[row,col])) <= as.numeric(as.character(MistagThresholds[Fam,4])),0,as.numeric(as.character(Detections[row,col])))
    }
  }
  rownames(Out) <- rownames(Detections) 
  colnames(Out) <- colnames(Detections)
  return(Out)
}

