#####################################################################################################
# Try it out on real data 

# load the data
df <- read.csv("~/Projects/CriticalMistagNB/Resubmission/SimulatedDataset.csv", row.names = 1)

CMs <- colSums(df[,31:40]) # list of mistags per control bin
Samples <- 40 # total number of dual-index bins (experimental samples plus control samples)
TotalSeqs <- sum(colSums(df)) # total sequences across all samples

# make mistags per taxon per dual-index bin data long form
Abundance <- rep(rowSums(df),10)
MTs <- c()
for (i in 31:40) {
  MTs <- c(MTs,df[,i])
} 
dfLong <- data.frame("Abundance"=Abundance,"MTs"=MTs)

# find best fit Poisson
# here, we assume that spatial clustering in the dual-index matrix (i.e. forward vs reverse tags)
# may lead to specific dual-index bins receiving more than their fair share of mistags
# so we find the best-fit Poisson distribution, which tells us about how non-uniform the data is
library(fitdistrplus)
y <- CMs 
y <- round((y/max(y))*40, 0) # standardize to discrete series with max = 40 (i.e. the number of bins)
f3 <- fitdist(y,"pois") # estimate lambda assuming Poisson 
plotdist(y,"pois",para=list(lambda=f3$estimate[1]))
hist(rpois(1000,lambda = f3$estimate[1]))
hist(y, breaks = 20)

# mistag probability if we just base it off the mean
# To estimate mistag probability in this way
# find average number of mistags per dual-index combo,
# multiply by number of dual-index combos to estimate total mistags,
# then P = Mistags/TotalSequences
MistagRate <- (round(mean(CMs),0)*40)/TotalSeqs

# however, we don't really know the mean mistags per dual-index bin
# and it very reasonably could be greater then estimated with the mean
# so it's better to calculate the upper end of the 95% confidence interval of the mean
# take the distribution of critical mistags per control dual-index bin
# and calculate from t-distribution
ks.test(CMs, "pnorm", mean=mean(CMs), sd=sd(CMs)) # test if data significantly deviates from normality
CMs95CL <- mean(CMs) + qt(0.975, df = length(CMs)-1)*(sd(CMs)/length(CMs)^(1/2))
MistagRate <- (CMs95CL*Samples)/TotalSeqs #mistags/sequences


###############################################################################################
# Now calculate 95th percentile of max mistags expected per species
# and 95th percentile of max mistags expetected per species per bin
# assuming uniform distribution of mistags across bins
df$Total <- rowSums(df)
SpeiciesAndAbundances <- data.frame("Taxa"=rownames(df),"Total"=df$Total)
Samples <- Samples
MistagRate <- MistagRate
Iters <- 25000
Lambda <- f3$estimate[1]
FDR <- 0.05

# use MaxMistagsPerBinU() function from CM_Functions.R
dfP <- MaxMistagsPerBinP(Samples,SpeiciesAndAbundances,MistagRate,Iters,Lambda,FDR)

# now plot results
plot(log(dfP$SeqsPerSpecies+1), log(dfP$MistagsAtFDR+1), type='l', main="Thresholds for Simulated Data", col="blue2",
     lwd=3, xlim=c(4,9),ylim=c(0,6), ylab="log(Mistags)", xlab="log(Total sequences per taxon)")
lines(log(dfP$SeqsPerSpecies+1), log(dfP$MaxMTsPerBin_Poisson+1), col="lightblue1", lwd=3)
points(log(rowSums(df[,1:30])), log(rowSums(df[,31:40])), pch=20, cex=1, col="gray")
points(log(dfLong$Abundance+1), log(dfLong$MTs+1), pch=20, cex=0.4)

# at FDR<=0.05, we expect fewer then 5 in 100 points to be above light blue line
# is that true?
dfLong$Expected <- NA
for (i in 1:nrow(dfLong)) {
  tSeqs <- dfLong[i,"Abundance"]
  dfLong[i,"Expected"] <- dfP[which(dfP$SeqsPerSpecies==tSeqs),"MaxMTsPerBin_Poisson"] }
dfLong$ExceedsThreshold <- ifelse(dfLong$MTs > dfLong$Expected, 1, 0)
sum(dfLong$ExceedsThreshold)/nrow(dfLong)
# can do fishers exact to test if significantly greater then expected
prop.test(sum(dfLong$ExceedsThreshold), nrow(dfLong), 0.05)

#####################################################################################################
# Use FilterDetections() function to filter data 
FilteredData <- FilterDetections(dfP, df[,1:30])

df <- as.matrix(df[,1:30]) # convert detections to matrix
df[(df>0)] <- 1 # convert detections to 1
FilteredData <- as.matrix(FilteredData) # convert to matrix
d <- (df > FilteredData) # Find cells which were filtered
d <- d*1 # change from T/F to 1/0
sum(d) # total detection removed
sum(d)/sum(df) # proportion of detections removed


