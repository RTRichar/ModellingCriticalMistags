#####################################################################################################
# Try it out on real data 

# load the data
# here, "NL" and "PC" in column names indicate 'no library negative' and 'positive control,' respectively
# there's no evidence of contamination levels differing between these control types
# so we can assume they both primarily contain mistags as opposed to lab contamination
# the fist column is the total number of sequences per taxon in experimental samples
# last column is total number of mistags per taxon across all control bins
# mistags per taxon per bin are shown in columns in between
df <- read.csv("~/Projects/CriticalMistagNB/Dataset1.csv")
CMs <- colSums(df[,3:14]) # list of mistags per control bin
Samples <- 42 # total number of dual-index bins
TotalSeqs <- sum(colSums(df[,c(2,15)]))

# make mistags per species per dual-index bin data long form
Abundance <- rep((df[,2]+df[,ncol(df)]),(ncol(df)-3))
MTs <- c()
for (i in 3:(ncol(df)-1)) {
  MTs <- c(MTs,df[,i])
} 
dfLong <- data.frame("Abundance"=Abundance,"MTs"=MTs)

# mistag probability = 0.00385 if we just base it off the mean
# To estimate mistag probability in this way
# find average number of mistags per dual-index combo,
# multiply by number of dual-index combos to estimate total mistags,
# then P = Mistags/TotalSequences
MistagRate <- (round(mean(CMs),0)*42)/TotalSeqs

# however, we don't really know the mean mistags per dual-index bin
# and it very reasonably could be greater then estimated with the mean
# so it's better to calculate the upper end of the 95% confidence interval of the mean
# take the distribution of critical mistags per control dual-index bin
# and calculate from t-distribution
ks.test(CMs, "pnorm", mean=mean(CMs), sd=sd(CMs)) # small sample size, but data does not significantly deviates from normality
CMs95CL <- mean(CMs) + qt(0.975, df = length(CMs)-1)*(sd(CMs)/length(CMs)^(1/2))
MistagRate <- (CMs95CL*Samples)/TotalSeqs #mistags/sequences


###############################################################################################
# Now calculate 95th percentile of max mistags expected per species
# and 95th percentile of max mistags expetected per species per bin
# assuming uniform distribution of mistags across bins
df$Total <- (df$Dat1_Exp+df$Dat1_Cont)
SpeiciesAndAbundances <- df[,c("Taxa","Total")]
Samples <- Samples
MistagRate <- MistagRate
Iters <- 25000

# use MaxMistagsPerBinU() function from CM_Functions.R
dfU <- MaxMistagsPerBinU(Samples,SpeiciesAndAbundances,MistagRate,Iters)

# now plot results
plot(log(dfU$SeqsPerSpecies+1), log(dfU$Mistags95thP+1), type='l', main="Data Set 1", col="blue2",
     lwd=3, xlim=c(0,11.2),ylim=c(0,6.1), ylab="log(Mistags)", xlab="log(Total sequences per taxon)")
lines(log(dfU$SeqsPerSpecies+1), log(dfU$MaxMTsPerBin_Uniform+1), col="lightblue1", lwd=3)
points(log(df$Dat1_Exp+1), log(df$Dat1_Cont+1), pch=20, cex=1, col="gray")
points(log(dfLong$Abundance+1), log(dfLong$MTs+1), pch=20, cex=0.4)


# at FDR<=0.05, we expect about 5 in 100 points to be above FDR limit
# is that true?
dfLong$Expected <- NA
for (i in 1:nrow(dfLong)) {
  tSeqs <- dfLong[i,"Abundance"]
  dfLong[i,"Expected"] <- unique(dfU[which(dfU$SeqsPerSpecies==tSeqs),"MaxMTsPerBin_Uniform"]) }
#dfTMP <- dfLong # with full dataset
dfTMP <- dfLong[which(dfLong$Abundance>=2000),] # with abundant taxa, because rare taxa are unlikely to produce mistag false positives
dfTMP$ExceedsThreshold <- ifelse(dfTMP$MTs > dfTMP$Expected, 1, 0)
sum(dfTMP$ExceedsThreshold)/nrow(dfTMP)
# can test if significantly greater then expected
prop.test(sum(dfTMP$ExceedsThreshold), nrow(dfTMP), 0.05)

####################################################################################################
####################################################################################################
####################################################################################################
# now try on a second dataset
df <- read.csv("~/Projects/CriticalMistagNB/Dataset2.csv")
CMs <- colSums(df[3:9]) # list of mistags per control bin
Samples <- 303 # total number of dual-index bins
TotalSeqs <- sum(colSums(df[,c(2,10)]))

# make mistags per species per dual-index bin data long form
Abundance <- rep((df[,2]+df[,ncol(df)]),(ncol(df)-3))
MTs <- c()
for (i in 3:(ncol(df)-1)) {
  MTs <- c(MTs,df[,i])
} 
dfLong <- data.frame("Abundance"=Abundance,"MTs"=MTs)

# calculate the upper end of the 95% confidence interval of the mean
# take the distribution of critical mistags per control dual-index bin
# and calculate from t-distribution
ks.test(CMs, "pnorm", mean=mean(CMs), sd=sd(CMs)) # test if data significantly deviates from normality
CMs95CL <- mean(CMs) + qt(0.975, df = length(CMs)-1)*(sd(CMs)/length(CMs)^(1/2))
MistagRate <- (CMs95CL*Samples)/TotalSeqs #mistags/sequences


###############################################################################################
# Now calculate 95th percentile of max mistags expected per species
# and 95th percentile of max mistags expected per species per bin
# assuming uniform distribution of mistags across bins
df$Total <- (df$Dat2_Exp+df$Dat2_Cont)
SpeiciesAndAbundances <- df[,c("Taxa","Total")]
Samples <- Samples
MistagRate <- MistagRate
Iters <- 25000

# use MaxMistagsPerBinU() function from CM_Functions.R
dfU <- MaxMistagsPerBinU(Samples,SpeiciesAndAbundances,MistagRate,Iters)

# now plot results
plot(log(dfU$SeqsPerSpecies+1), log(dfU$Mistags95thP+1), type='l', main="Data Set 2", col="blue2",
     lwd=3, ylab="log(Mistags)", xlab="log(Total sequences per taxon)", xlim = c(0,16))
lines(log(dfU$SeqsPerSpecies+1), log(dfU$MaxMTsPerBin_Uniform+1), col="lightblue1", lwd=3)
points(log(df$Dat2_Exp+1), log(df$Dat2_Cont+1), pch=20, cex=1, col="gray")
points(log(dfLong$Abundance+1), log(dfLong$MTs+1), pch=20, cex=0.4)

# at FDR<=0.05, we expect fewer then 5 in 100 points to be above FDR limit
# is that true?
dfLong$Expected <- NA
for (i in 1:nrow(dfLong)) {
  tSeqs <- dfLong[i,"Abundance"]
  dfLong[i,"Expected"] <- unique(dfU[which(dfU$SeqsPerSpecies==tSeqs),"MaxMTsPerBin_Uniform"]) }
#dfTMP <- dfLong
dfTMP <- dfLong[which(dfLong$Abundance>=2000),]
dfTMP$ExceedsThreshold <- ifelse(dfTMP$MTs > dfTMP$Expected, 1, 0)
sum(dfTMP$ExceedsThreshold)/nrow(dfTMP)
# can test if significantly greater then expected
prop.test(sum(dfTMP$ExceedsThreshold), nrow(dfTMP), 0.05)
