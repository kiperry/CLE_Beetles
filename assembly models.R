#_________________________________________________________________________________
# Testing patterns of community assembly among habitats
# Vacant lots, pocket prairies, forests, old fields
#
# Cleveland, Ohio and surrounding counties
#
# Ground beetle communities: null model analyses
#
# 10 October 2025
#_________________________________________________________________________________

# Import abundance data ----
library(readxl)

df <- read_excel("2019_Ground_Beetles_Raw Data.xlsx",
                 sheet = "Raw_SpeciesCode",
                 range = "A1:CU565")

df$Trmt_Code <- as.factor(df$Trmt_Code)
df$Site <- as.factor(df$Site)

summary(df)
str(df)

# need to pool the abundance data across date and trap
library(reshape2)

df.long <- melt(df, id = c("Treatment", "Trmt_Code", "Site", "Month", "Date_Set",
                           "Date_Collected", "Pitfall", "Year"))
str(df.long)

# double check only species in the variable column
levels(df.long$variable)
names(df.long)[9] <- "Species"

# now create data frame for community analyses
df2 <- dcast(df.long, Trmt_Code + Site + Year ~ Species, sum, na.rm = TRUE)

summary(df2)

# pull out treatments of interest
levels(df2$Trmt_Code)

MP <- df2[which(df2$Trmt_Code == "MP"),]
OF <- df2[which(df2$Trmt_Code == "OF"),]
PP <- df2[which(df2$Trmt_Code == "PP"),]
VL <- df2[which(df2$Trmt_Code == "VL"),]

df3 <- rbind(MP, OF, PP, VL)
str(df3)
df3 <- droplevels(df3)
str(df3)
colSums(df3[4:94])

# copy data frame
a <- df3[4:94]

colnames(a)

# change data set to presence/absence
a[a > 0] <- 1
str(a)
rowSums(a)
colSums(a)

# create a vector with the column sums for each species
# species not collected will have a 0
sp <- colSums(a)
sp

# removes any columns (i.e., species) that were not collected
a <- a[, colSums(a != 0) > 0]
str(a)
colSums(a)

#_______________________________________________________________________________

# Import trait data ----
t <- read_excel("Ground_Beetle_Trait_Data.xlsx",
                 sheet = "Relative data with codes",
                 range = "A1:L92")

str(t)
summary(t)

# add the sp vector as a column in the trait matrix, shows which species
# were collected in Madison and which were absent (i.e., with a 0)
t$sp <- sp

# use the sp values to remove rows of species not collected in Madison
# then remove the column because we don't need it anymore
t <- t[t$sp != 0, ]
t <- t[,-13]
str(t)

# need to remove species names column
t2 <- t[,-1]

# assign species codes as row names then remove column
library(tidyverse)
t2 <- t2 %>%
  column_to_rownames(var = "code")

rownames(t2)

# trim the trait data set by identifying traits that are highly correlated
# or lack sufficient variance among species
names(t2)
str(t2)
plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

# legs and antennae are correlated, remove femur length and leg (tibia + femur) length
# keep tibia length which isn't correlated with antennae length

t3 <- t2[,-6]# femur length
t3 <- t3[,-6]# leg length (tibia + femur)
t3 <- t3[,-7]# origin (not enough variation among species)

plot(t3)
cor(t3, method = c("pearson"), use = "complete.obs") # looks good

#copy data set in case transformations are needed
t4 <- t3

str(t4)
t4$dispersal <- as.factor(t4$dispersal)
t4$activity_period <- as.factor(t4$activity_period)
str(t4)

#check traits for normality
hist(t4$bs)

hist(t4$rhw)
hist(log(t4$rhw + 1))
t4$rhw <- log(t4$rhw + 1)
hist(t4$rhw)

hist(t4$rew)
hist(log(t4$rew + 1))

hist(t4$ral)
hist(log(t4$ral + 1))

hist(t4$rtl)
hist(log(t4$rtl + 1))
t4$rtl <- log(t4$rtl + 1)

# double check that all species are present in both data sets
intersect(colnames(a), rownames(t4))
names(a)

# double check if a species is present in one data set but not the other
setdiff(colnames(a), rownames(t4))
setdiff(rownames (t4), colnames(a))

# double check all species names are in the same order
rownames(t4) == colnames(a)

# need to remove species with no trait data in both data sets
t5 <- t4[-31,]# acpu
t5 <- t5[-55,]# bemi

a2 <- a[,-31]# acpu
a2 <- a2[,-55]# bemi

rownames(t5) == colnames(a2)

# good to go!

#_______________________________________________________________________________

# Observed community metrics ----

library(FD)
library(gawdis)
library(ade4)
library(betapart)
library(ape)
library(picante)

####
# observed functional CWM
cwm.obs <- functcomp(t5, as.matrix(a2), CWM.type = "all")
cwm.obs

####
# observed taxonomic beta diversity
# create beta part object for analyses
gb.core <- betapart.core(a2)

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
gb.dist <- beta.pair(gb.core, index.family = "sorensen")
str(gb.dist)

####
# observed functional beta diversity
# weight traits on the head together to limit their total influence on the metric
tdis <- gawdis(t5, w.type = "optimized", opti.maxiter = 300,
               groups.weight = T, groups = c(1, 2, 2, 2, 3, 4, 5))
attr(tdis, "correls")
attr(tdis, "weights")

# calculate the distance matrix
pcoGB <- dudi.pco(tdis, scannf = FALSE, nf = 4)
scatter(pcoGB)

pcoGB$li
sum(pcoGB$eig[1:4]) / sum(pcoGB$eig)
sum(pcoGB$eig[1:3]) / sum(pcoGB$eig)
sum(pcoGB$eig[1:2]) / sum(pcoGB$eig)

# check correlations among axes and traits
str(t5)
cor(pcoGB$li, t5[,1:5], use = "complete.obs")

# due to number of beetle species at each site (one site has only 3 bee species),
# we can only use first two axes of PCoA for functional diversity metrics
t.ax <- as.matrix(pcoGB$li[1:2])

# returns pairwise between-site values of each functional beta-diversity component
gb.fun <- functional.beta.pair(a2, t.ax, index.family = "sorensen")
str(gb.fun)

# observed functional alpha diversity
# run rao function first
gb.rao <- Rao(sample = t(a2), dfunc = tdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
gb.falpha <- gb.rao$FD$Alpha
gb.falpha

# create a copy of the abundance matrix to randomize
ra <- a2

# Null model ----

# run the null model with 999 iterations
numberReps <- 999

# create empty matrices to store the results of each iteration of the null model:

# cwms
nbs <- nrhw <- nrew <- nral <- nrtl <- ndisp_0 <- ndisp_1 <- ndisp_2 <- nactivity_0 <- nactivity_1 <- nactivity_2  <- matrix(NA,
              nrow = nrow(a), ncol = numberReps, dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# taxonomic beta diversity
nbsim <- nbsne <- nbsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# functional alpha and beta diversity
nfalpha <- nfsim <- nfsne <- nfsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                             dimnames = list(rownames(a), paste0("n", 1:numberReps)))

#create null model for each repetition:

for(i in 1:numberReps){
  print(i) 
  
  # randomize trait matrix
  ntraits <- t5[sample(1:nrow(t5)),]
  rownames(ntraits) <- rownames(t5)
  
  # randomize abundance matrix
  # independent swap constrains by species richness and frequency
  spGB <- randomizeMatrix(samp = ra, null.model = "independentswap")
  print(rownames(ntraits) == colnames(spGB)) 
  
  # randomize trait distance matrix
  ntdis <- gawdis(ntraits, w.type = "optimized", opti.maxiter = 300,
                  groups.weight = T, groups = c(1, 2, 2, 2, 3, 4, 5))

  # CWM calculations
  cwm.null <- functcomp(x = ntraits, a = as.matrix(spGB), CWM.type = "all")
  nbs[,i] <- cwm.null$bs
  nrhw[,i] <- cwm.null$rhw
  nrew[,i] <- cwm.null$rew
  nral[,i] <- cwm.null$ral
  nrtl[,i] <- cwm.null$rtl
  ndisp_0[,i] <- cwm.null$dispersal_0
  ndisp_1[,i] <- cwm.null$dispersal_1
  ndisp_2[,i] <- cwm.null$dispersal_2
  nactivity_0[,i] <- cwm.null$activity_period_0
  nactivity_1[,i] <- cwm.null$activity_period_1
  nactivity_2[,i] <- cwm.null$activity_period_2
  
  # Functional alpha diversity
  nrao <- Rao(sample = t(spGB), dfunc = ntdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
  nfalpha[,i] <- nrao$FD$Alpha
  
  # Taxonomic beta diversity indices
  ngb.core <- betapart.core(spGB)
  ngb.dist <- beta.pair(ngb.core, index.family = "sorensen")
  nsim.dist <- as.matrix(ngb.dist$beta.sim)
  nsne.dist <- as.matrix(ngb.dist$beta.sne)
  nsor.dist <- as.matrix(ngb.dist$beta.sor)
  nbsim[,i] <- colMeans(nsim.dist)
  nbsne[,i] <- colMeans(nsne.dist)
  nbsor[,i] <- colMeans(nsor.dist)
  
  # Functional beta diversity indices
  npcoGB <- dudi.pco(sqrt(ntdis), scannf = FALSE, nf = 2)
  nt <- as.matrix(npcoGB$li)
  ngb.fun <- functional.beta.pair(spGB, nt, index.family = "sorensen")
  nfsim.dist <- as.matrix(ngb.fun$funct.beta.sim)
  nfsne.dist <- as.matrix(ngb.fun$funct.beta.sne)
  nfsor.dist <- as.matrix(ngb.fun$funct.beta.sor)
  nfsim[,i] <- colMeans(nfsim.dist)
  nfsne[,i] <- colMeans(nfsne.dist)
  nfsor[,i] <- colMeans(nfsor.dist)
}

#_______________________________________________________________________________

write.csv(nbs, file = "Null_models_2Axes/nbs.csv")
write.csv(nrhw, file = "Null_models_2Axes/nrhw.csv")
write.csv(nrew, file = "Null_models_2Axes/nrew.csv")
write.csv(nral, file = "Null_models_2Axes/nral.csv")
write.csv(nrtl, file = "Null_models_2Axes/nrtl.csv")
write.csv(ndisp_0, file = "Null_models_2Axes/ndisp_0.csv")
write.csv(ndisp_1, file = "Null_models_2Axes/ndisp_1.csv")
write.csv(ndisp_2, file = "Null_models_2Axes/ndisp_2.csv")
write.csv(nactivity_0, file = "Null_models_2Axes/nactivity_0.csv")
write.csv(nactivity_1, file = "Null_models_2Axes/nactivity_1.csv")
write.csv(nactivity_2, file = "Null_models_2Axes/nactivity_2.csv")

write.csv(nbsim, file = "Null_models_2Axes/bb_tbeta_sim.csv")
write.csv(nbsne, file = "Null_models_2Axes/bb_tbeta_sne.csv")
write.csv(nbsor, file = "Null_models_2Axes/bb_tbeta_sor.csv")

write.csv(nfalpha, file = "Null_models_2Axes/bb_falpha.csv")
write.csv(nfsim, file = "Null_models_2Axes/bb_fbeta_sim.csv")
write.csv(nfsne, file = "Null_models_2Axes/bb_fbeta_sne.csv")
write.csv(nfsor, file = "Null_models_2Axes/bb_fbeta_sor.csv")

#_______________________________________________________________________________

# load the data sets

nbs <- read.csv("Null_models_2Axes/nbs.csv", row.names=1)
nrhw <- read.csv("Null_models_2Axes/nrhw.csv", row.names=1)
nrew <- read.csv("Null_models_2Axes/nrew.csv", row.names=1)
nral <- read.csv("Null_models_2Axes/nral.csv", row.names=1)
nrtl <- read.csv("Null_models_2Axes/nrtl.csv", row.names=1)
ndisp_0 <- read.csv("Null_models_2Axes/ndisp_0.csv", row.names=1)
ndisp_1 <- read.csv("Null_models_2Axes/ndisp_1.csv", row.names=1)
ndisp_2 <- read.csv("Null_models_2Axes/ndisp_2.csv", row.names=1)
nactivity_0 <- read.csv("Null_models_2Axes/nactivity_0.csv", row.names=1)
nactivity_1 <- read.csv("Null_models_2Axes/nactivity_1.csv", row.names=1)
nactivity_2 <- read.csv("Null_models_2Axes/nactivity_2.csv", row.names=1)

nbsim <- read.csv("Null_models_2Axes/bb_tbeta_sim.csv", row.names=1)
nbsne <- read.csv("Null_models_2Axes/bb_tbeta_sne.csv", row.names=1)
nbsor <- read.csv("Null_models_2Axes/bb_tbeta_sor.csv", row.names=1)

nfalpha <- read.csv("Null_models_2Axes/bb_falpha.csv", row.names=1)
nfsim <- read.csv("Null_models_2Axes/bb_fbeta_sim.csv", row.names=1)
nfsne <- read.csv("Null_models_2Axes/bb_fbeta_sne.csv", row.names=1)
nfsor <- read.csv("Null_models_2Axes/bb_fbeta_sor.csv", row.names=1)

#_______________________________________________________________________________

# SES Calculations ----

#calculate standardized effect sizes (SES) for each trait and index
#the effect size is the difference between the observed value and the expected one
#then divide the effect size by the standard deviation of the null distribution to get the standardized effect size
#allows comparison among sites with different numbers of species

## calculate SES values for each metric

## community weighted means ----

## body size
SES_bs <- (cwm.obs$bs - apply(nbs, MARGIN = 1, mean)) / apply(nbs, MARGIN = 1, sd, na.rm=T)
SES_bs

## head width
SES_hw <- (cwm.obs$rhw - apply(nrhw, MARGIN = 1, mean)) / apply(nrhw, MARGIN = 1, sd, na.rm=T)
SES_hw

## eye width
SES_ew <- (cwm.obs$rew - apply(nrew, MARGIN = 1, mean)) / apply(nrew, MARGIN = 1, sd, na.rm=T)
SES_ew

## antennae length
SES_al <- (cwm.obs$ral - apply(nral, MARGIN = 1, mean)) / apply(nral, MARGIN = 1, sd, na.rm=T)
SES_al

## tibia length
SES_tl <- (cwm.obs$rtl - apply(nrtl, MARGIN = 1, mean)) / apply(nrtl, MARGIN = 1, sd, na.rm=T)
SES_tl

## dispersal 0 = brachypterous
SES_disp_0 <- (cwm.obs$dispersal_0 - apply(ndisp_0, MARGIN = 1, mean)) / apply(ndisp_0, MARGIN = 1, sd, na.rm=T)
SES_disp_0

## dispersal 1 = dimorphic
SES_disp_1 <- (cwm.obs$dispersal_1 - apply(ndisp_1, MARGIN = 1, mean)) / apply(ndisp_1, MARGIN = 1, sd, na.rm=T)
SES_disp_1

## dispersal 2 = macropterous
SES_disp_2 <- (cwm.obs$dispersal_2 - apply(ndisp_2, MARGIN = 1, mean)) / apply(ndisp_2, MARGIN = 1, sd, na.rm=T)
SES_disp_2

## activity 0 = nocturnal
SES_activity_0 <- (cwm.obs$activity_period_0 - apply(nactivity_0, MARGIN = 1, mean)) / apply(nactivity_0, MARGIN = 1, sd, na.rm=T)
SES_activity_0

## activity 1 = both
SES_activity_1 <- (cwm.obs$activity_period_1 - apply(nactivity_1, MARGIN = 1, mean)) / apply(nactivity_1, MARGIN = 1, sd, na.rm=T)
SES_activity_1

## activity 2 = diurnal
SES_activity_2 <- (cwm.obs$activity_period_2 - apply(nactivity_2, MARGIN = 1, mean)) / apply(nactivity_2, MARGIN = 1, sd, na.rm=T)
SES_activity_2

## taxonomic beta diversity ----

beta.sor <- as.matrix(gb.dist$beta.sor)
beta.sor <- colMeans(beta.sor)

beta.sim <- as.matrix(gb.dist$beta.sim)
beta.sim <- colMeans(beta.sim)

beta.sne <- as.matrix(gb.dist$beta.sne)
beta.sne <- colMeans(beta.sne)

beta.t <- data.frame(beta.sor, beta.sim, beta.sne)

## taxonomic diveristy - beta sor
SES_bsor <- (beta.t$beta.sor - apply(nbsor, MARGIN = 1, mean)) / apply(nbsor, MARGIN = 1, sd, na.rm=T)
SES_bsor

## taxonomic diveristy - beta sim
SES_bsim <- (beta.t$beta.sim - apply(nbsim, MARGIN = 1, mean)) / apply(nbsim, MARGIN = 1, sd, na.rm=T)
SES_bsim

## taxonomic diveristy - beta sne
SES_bsne <- (beta.t$beta.sne - apply(nbsne, MARGIN = 1, mean)) / apply(nbsne, MARGIN = 1, sd, na.rm=T)
SES_bsne

## functional beta diversity ----
falpha <- as.matrix(gb.falpha)

fbeta.sor <- as.matrix(gb.fun$funct.beta.sor)
fbeta.sor <- colMeans(fbeta.sor)

fbeta.sim <- as.matrix(gb.fun$funct.beta.sim)
fbeta.sim <- colMeans(fbeta.sim)

fbeta.sne <- as.matrix(gb.fun$funct.beta.sne)
fbeta.sne <- colMeans(fbeta.sne)

beta.f <- data.frame(fbeta.sor, fbeta.sim, fbeta.sne, falpha)

## functional alpha diversity
SES_falpha <- (beta.f$falpha - apply(nfalpha, MARGIN = 1, mean)) / apply(nfalpha, MARGIN = 1, sd, na.rm=T)
SES_falpha

## functional diveristy - beta sor
SES_fbsor <- (beta.f$fbeta.sor - apply(nfsor, MARGIN = 1, mean)) / apply(nfsor, MARGIN = 1, sd, na.rm=T)
SES_fbsor

## functional diveristy - beta sim
SES_fbsim <- (beta.f$fbeta.sim - apply(nfsim, MARGIN = 1, mean)) / apply(nfsim, MARGIN = 1, sd, na.rm=T)
SES_fbsim

## functional diveristy - beta sne
SES_fbsne <- (beta.f$fbeta.sne - apply(nfsne, MARGIN = 1, mean)) / apply(nfsne, MARGIN = 1, sd, na.rm=T)
SES_fbsne


## combine all indices ----
SES <- cbind(SES_bs, SES_hw, SES_ew, SES_al, SES_tl, SES_disp_0, SES_disp_1, SES_disp_2,
             SES_activity_0, SES_activity_1, SES_activity_2, SES_bsor, SES_bsim, SES_bsne,
             SES_fbsor, SES_fbsim, SES_fbsne, SES_falpha)
SES <- as.data.frame(SES)
str(SES)

## add in treatment categories
SES$trmt <- df3$Trmt_Code
SES$site <- df3$Site
str(SES)

write.csv(SES, file = "SES_CLE.csv")

## load the SES data set ----
SES <- read.csv("SES_CLE.csv", row.names=1)

# Comparison: Null Communities ----

## pull out data for each treatment
vl <- SES[which(SES$trmt == "VL"),]
vl <- droplevels(vl)
str(vl)

pp <- SES[which(SES$trmt == "PP"),]
pp <- droplevels(pp)
str(pp)

of <- SES[which(SES$trmt == "OF"),]
of <- droplevels(of)
str(of)

mp <- SES[which(SES$trmt == "MP"),]
mp <- droplevels(mp)
str(mp)


## taxonomic diveristy - beta sor
hist(SES$SES_bsor)
plot(SES$SES_bsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
wilcox.test(vl$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(pp$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(of$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(mp$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)


## taxonomic diversity - beta sim
hist(SES$SES_bsim)
plot(SES$SES_bsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
wilcox.test(vl$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(pp$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(of$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(mp$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)


## taxonomic diversity - beta sne
hist(SES$SES_bsne)
plot(SES$SES_bsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
wilcox.test(vl$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(pp$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(of$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(mp$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)


## functional alpha diversity
hist(SES$SES_falpha)
plot(SES$SES_falpha, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_falpha ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
wilcox.test(vl$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(pp$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(of$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(mp$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)


## functional beta diversity - beta sor
hist(SES$SES_fbsor)
plot(SES$SES_fbsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
wilcox.test(vl$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(pp$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(of$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(mp$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)


## functional beta diversity - beta sim
hist(SES$SES_fbsim)
plot(SES$SES_fbsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
wilcox.test(vl$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(pp$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(of$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(mp$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)


## functional beta diversity - beta sne
hist(SES$SES_fbsne)
plot(SES$SES_fbsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
wilcox.test(vl$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(pp$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(of$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
wilcox.test(mp$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)


# Comparison: Among Treatments ----
library(emmeans)

# body size
hist(SES$SES_bs)
dotchart(SES$SES_bs, group = SES$trmt, pch = 19)
plot(SES$SES_bs, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bs ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_bs <- lm(SES_bs ~ trmt, data = SES)
summary(mod_bs)
anova(mod_bs)
emmeans(mod_bs, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_bs))
par(mfrow = c(1, 1))

# head width
hist(SES$SES_hw)
dotchart(SES$SES_hw, group = SES$trmt, pch = 19)
plot(SES$SES_hw, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_hw ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_hw <- lm(SES_hw ~ trmt, data = SES)
summary(mod_hw)
anova(mod_hw)
emmeans(mod_hw, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_hw))
par(mfrow = c(1, 1))

# eye width
hist(SES$SES_ew)
dotchart(SES$SES_ew, group = SES$trmt, pch = 19)
plot(SES$SES_ew, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_ew ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_ew <- lm(SES_ew ~ trmt, data = SES)
summary(mod_ew)
anova(mod_ew)
emmeans(mod_ew, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_ew))
par(mfrow = c(1, 1))

# antennae length
hist(SES$SES_al)
dotchart(SES$SES_al, group = SES$trmt, pch = 19)
plot(SES$SES_al, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_al ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_al <- lm(SES_al ~ trmt, data = SES)
summary(mod_al)
anova(mod_al)
emmeans(mod_al, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_al))
par(mfrow = c(1, 1))

# tibia length
hist(SES$SES_tl)
dotchart(SES$SES_tl, group = SES$trmt, pch = 19)
plot(SES$SES_tl, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_tl ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_tl <- lm(SES_tl ~ trmt, data = SES)
summary(mod_tl)
anova(mod_tl)
emmeans(mod_tl, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_tl))
par(mfrow = c(1, 1))

# brachypterous
hist(SES$SES_disp_0)
dotchart(SES$SES_disp_0, group = SES$trmt, pch = 19)
plot(SES$SES_disp_0, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_disp_0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_disp_0 <- lm(SES_disp_0 ~ trmt, data = SES)
summary(mod_disp_0)
anova(mod_disp_0)
emmeans(mod_disp_0, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_disp_0))
par(mfrow = c(1, 1))

# macropterous
hist(SES$SES_disp_2)
dotchart(SES$SES_disp_2, group = SES$trmt, pch = 19)
plot(SES$SES_disp_2, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_disp_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_disp_2 <- lm(SES_disp_2 ~ trmt, data = SES)
summary(mod_disp_2)
anova(mod_disp_2)
emmeans(mod_disp_2, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_disp_2))
par(mfrow = c(1, 1))

# nocturnal
hist(SES$SES_activity_0)
dotchart(SES$SES_activity_0, group = SES$trmt, pch = 19)
plot(SES$SES_activity_0, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_activity_0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_activity_0 <- lm(SES_activity_0 ~ trmt, data = SES)
summary(mod_activity_0)
anova(mod_activity_0)
emmeans(mod_activity_0, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_activity_0))
par(mfrow = c(1, 1))

# diurnal
hist(SES$SES_activity_2)
dotchart(SES$SES_activity_2, group = SES$trmt, pch = 19)
plot(SES$SES_activity_2, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_activity_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_activity_2 <- lm(SES_activity_2 ~ trmt, data = SES)
summary(mod_activity_2)
anova(mod_activity_2)
emmeans(mod_activity_2, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_activity_2))
par(mfrow = c(1, 1))

# taxonomic beta-diversity - beta sor
hist(SES$SES_bsor)
dotchart(SES$SES_bsor, group = SES$trmt, pch = 19)
plot(SES$SES_bsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_bsor <- lm(SES_bsor ~ trmt, data = SES)
summary(mod_bsor)
anova(mod_bsor)
emmeans(mod_bsor, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_bsor))
par(mfrow = c(1, 1))

# functional alpha diversity
hist(SES$SES_falpha)
dotchart(SES$SES_falpha, group = SES$trmt, pch = 19)
plot(SES$SES_falpha, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_falpha ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_falpha <- lm(SES_falpha ~ trmt, data = SES)
summary(mod_falpha)
anova(mod_falpha)
emmeans(mod_falpha, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_falpha))
par(mfrow = c(1, 1))

# functional beta diversity - beta sor
hist(SES$SES_fbsor)
dotchart(SES$SES_fbsor, group = SES$trmt, pch = 19)
plot(SES$SES_fbsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

mod_fbsor <- lm(SES_fbsor ~ trmt, data = SES)
summary(mod_fbsor)
anova(mod_fbsor)
emmeans(mod_fbsor, pairwise ~ trmt)
par(mfrow = c(2, 2))
(plot(mod_fbsor))
par(mfrow = c(1, 1))

## Figures----

png("Figures/SES Beta diversity.png", width = 2800, height = 1000, pointsize = 30)

par(mfrow=c(1,2)) # arrange figures in one row and four columns
par(mar=c(5,8,4,2))

boxplot(SES_bsor ~ trmt, data = SES, col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"),
        ylim = c(-2.5,3.5), ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2, cex.lab = 1.8, cex.axis = 1.8)
stripchart(SES_bsor ~ trmt, data = SES, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES_fbsor ~ trmt, data = SES, col = c("seagreen4", "seagreen2", "lightgoldenrod2", "goldenrod2" ),
        ylim = c(-2.5,3.5), ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        main = "", cex.lab = 1.8, cex.axis = 1.8)
stripchart(SES_fbsor ~ trmt, data = SES, pch = 19, cex = 2, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

dev.off()


png("Figures/SES Dispersal.png", width = 2800, height = 1000, pointsize = 30)

par(mfrow=c(1,2)) # arrange figures in one row and four columns
par(mar=c(5,8,4,2))

boxplot(SES_disp_0 ~ trmt, data = SES, col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"),
        ylim = c(-2.5,3.5), ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2, cex.lab = 1.8, cex.axis = 1.8)
stripchart(SES_disp_0 ~ trmt, data = SES, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES_disp_2 ~ trmt, data = SES, col = c("seagreen4", "seagreen2", "lightgoldenrod2", "goldenrod2" ),
        ylim = c(-2.5,3.5), ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        main = "", cex.lab = 1.8, cex.axis = 1.8)
stripchart(SES_disp_2 ~ trmt, data = SES, pch = 19, cex = 2, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

dev.off()


png("Figures/SES Activity.png", width = 2800, height = 1000, pointsize = 30)

par(mfrow=c(1,2)) # arrange figures in one row and four columns
par(mar=c(5,8,4,2))

boxplot(SES_activity_0 ~ trmt, data = SES, col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"),
        ylim = c(-2.5,3.5), ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2, cex.lab = 1.8, cex.axis = 1.8)
stripchart(SES_activity_0 ~ trmt, data = SES, pch = 19, cex = 2, add = TRUE, 
           vertical = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES_activity_2 ~ trmt, data = SES, col = c("seagreen4", "seagreen2", "lightgoldenrod2", "goldenrod2" ),
        ylim = c(-2.5,3.5), ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        main = "", cex.lab = 1.8, cex.axis = 1.8)
stripchart(SES_activity_2 ~ trmt, data = SES, pch = 19, cex = 2, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

dev.off()

# Community composition ----

# taxonomic
# run the nonmetric multidimensional scaling model
nmds.tax <- metaMDS(gb.dist$beta.sor, trymax = 500, autotransform = TRUE, k = 2)
nmds.tax # stress is quality of fit
stressplot(nmds.tax)
plot(nmds.tax) # basic plot with no treatment distinctions

# plot the NMDS model
ordiplot(nmds.tax, disp = "sites", type = "n", xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
points(nmds.tax, dis = "sites", select = which(df3$Trmt_Code=="MP"), pch = 17, cex = 2, col = "seagreen4")
points(nmds.tax, dis = "sites", select = which(df3$Trmt_Code=="OF"), pch = 18, cex = 2, col = "lightgoldenrod2")
points(nmds.tax, dis = "sites", select = which(df3$Trmt_Code=="PP"), pch = 15, cex = 2, col = "pink1")
points(nmds.tax, dis = "sites", select = which(df3$Trmt_Code=="VL"), pch = 16, cex = 2, col = "seagreen2")
ordiellipse(nmds.tax, df3$Trmt_Code, draw = "lines", col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)

legend("topleft", legend = c("Forest", "Old Field", "Pocket Prarie", "Vacant Lot"),
       pch = c(17, 18, 15, 16), cex = 1.5, bty = "n", col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"))

## Test for differences in predator composition among treatments

# PERMANOVA tests whether the group centroid of communities differs among groups
# in multivariate space (e.g. different community composition)
adonis2(gb.dist$beta.sor ~ df3$Trmt_Code, permutations = 999)

library(pairwiseAdonis)
pairwise.adonis(gb.dist$beta.sor, df3$Trmt_Code)

# BETADISPER tests whether the dispersion of a group from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogeneity of variances
tax.beta <- betadisper(gb.dist$beta.sor, df3$Trmt_Code, type = c("median"))
anova(tax.beta)
plot(tax.beta)
boxplot(tax.beta, ylab = "Distance to median")
TukeyHSD(tax.beta, which = "group", conf.level = 0.95)


# functional
# run the nonmetric multidimensional scaling model
nmds.fun <- metaMDS(gb.fun$funct.beta.sor, trymax = 500, autotransform = TRUE, k = 2)
nmds.fun # stress is quality of fit
stressplot(nmds.fun)
plot(nmds.fun) # basic plot with no treatment distinctions

# plot the NMDS model
ordiplot(nmds.fun, disp = "sites", type = "n", xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
points(nmds.fun, dis = "sites", select = which(df3$Trmt_Code=="MP"), pch = 17, cex = 2, col = "seagreen4")
points(nmds.fun, dis = "sites", select = which(df3$Trmt_Code=="OF"), pch = 18, cex = 2, col = "lightgoldenrod2")
points(nmds.fun, dis = "sites", select = which(df3$Trmt_Code=="PP"), pch = 15, cex = 2, col = "pink1")
points(nmds.fun, dis = "sites", select = which(df3$Trmt_Code=="VL"), pch = 16, cex = 2, col = "seagreen2")
ordiellipse(nmds.fun, df3$Trmt_Code, draw = "lines", col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)

legend("topleft", legend = c("Forest", "Old Field", "Pocket Prarie", "Vacant Lot"),
       pch = c(17, 18, 15, 16), cex = 1.5, bty = "n", col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"))

## Test for differences in predator composition among treatments

# PERMANOVA tests whether the group centroid of communities differs among groups
# in multivariate space (e.g. different community composition)
adonis2(gb.fun$funct.beta.sor ~ df3$Trmt_Code, permutations = 999)

library(pairwiseAdonis)
pairwise.adonis(gb.fun$funct.beta.sor, df3$Trmt_Code)

# BETADISPER tests whether the dispersion of a group from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogeneity of variances
fun.beta <- betadisper(gb.fun$funct.beta.sor, df3$Trmt_Code, type = c("median"))
anova(fun.beta)
plot(fun.beta)
boxplot(fun.beta, ylab = "Distance to median")
TukeyHSD(fun.beta, which = "group", conf.level = 0.95)


# PLSCA - Landscape variables & traits

