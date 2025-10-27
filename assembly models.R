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

df <- read_excel("2019_Ground_Beetles_Raw Data.xlsx",
                 sheet = "Raw_SpeciesCode",
                 range = "A1:CV565")

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
colSums(df3[4:95])

# copy data frame
a <- df3[4:95]

# create a vector with the column sums for each species
# species not collected will have a 0
sp <- colSums(a)
sp

# removes any columns (i.e. species) that were not collected
a <- a[, colSums(a != 0) > 0]
str(a)

#_______________________________________________________________________________
# DO THIS ONCE I ADD TRAIT DATA

# Import trait data ----

# add the sp vector as a column in the trait matrix, shows which species
# were collected in Madison and which were absent (i.e. with a 0)
t$sp <- sp
t

# use the sp values to remove rows of species not collected in Madison
# then remove the column because we don't need it anymore
t <- t[t$sp != 0, ]
t <- t[,-28]
str(t)

# trim the trait dataset by identifying traits that are highly correlated
# or lack sufficient variance among species
names(t)
str(t)
plot(t)
cor(t, method = c("pearson"), use = "complete.obs")

# lecty, nest construction, sociality, activity, parasitism, and pollen transport lack sufficient
# variance among bumble bee species - remove these traits
t2 <- t[,-8]#lecty
t2 <- t2[,-9]#nest construction
t2 <- t2[,-9]#sociality
t2 <- t2[,-9]#activity
t2 <- t2[,-9]#parasitism
t2 <- t2[,-9]#pollen transport

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

t2 <- t2[,-21]#wt

t2 <- t2[,-20]#corbicula width

t2 <- t2[,-10]#wing marginal cell length (keeping inter-tegular distance)
t2 <- t2[,-10]#wing width

t2 <- t2[,-12]#eye width

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

# remove the body size variables from the literature except body length variances
# these are correlated with each other and several other traits
t2 <- t2[,-1]#average queen body length
t2 <- t2[,-2]#average male body length
t2 <- t2[,-3]#average worker body length

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

t2 <- t2[,-10] #scape length
# will group other traits on head to limit their influence with functional diversity

#copy dataset in case transformataions are needed
t3 <- t2

str(t3)
t3$nestl <- as.factor(t3$nestl)
t3$tl <- as.factor(t3$tl)
str(t3)

#check traits for normality
hist(t2$qbl_var)
hist(log(t2$qbl_var))#
t3$qbl_var <- log(t3$qbl_var + 1)

hist(t2$mbl_var)
hist(log(t2$mbl_var))#
t3$mbl_var <- log(t3$mbl_var + 1)

hist(t2$wbl_var)
hist(log(t2$wbl_var))

hist(t2$it)
hist(log(t2$it))#
t3$it <- log(t3$it + 1)

# double check that all species are present in both data sets
intersect(colnames(a), rownames(t3))
intersect(colnames(a), p2$tip.label)
intersect(p2$tip.label, rownames(t3))
names(a)

# double check if a species is present in one data set but not the other
setdiff(colnames(a), rownames(t3))
setdiff(rownames (t3), colnames(a))
setdiff(rownames (t3), p2$tip.label)
setdiff(colnames(a), p2$tip.label)

# double check all species names are in the same order
rownames(t3) == colnames(a) 
colnames(a) == p2$tip.label
rownames(t3) == p2$tip.label

#_______________________________________________________________________________

# Observed community metrics ----

library(FD)
library(gawdis)

####
# observed functional CWM
cwm.obs <- functcomp(t3, as.matrix(a), CWM.type = "all")
cwm.obs

####
# observed taxonomic beta diversity
dis.matrix.t <- vegdist(a, method = "bray")

####
# observed functional beta diversity
# weight traits on the head together to limit their total influence on the metric
tdis <- gawdis(t3, w.type = "optimized", opti.maxiter = 300,
               groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
attr(tdis, "correls")
attr(tdis, "weights")

# calculate the distance matrix
pcoBB <- dudi.pco(sqrt(tdis), scannf = FALSE, nf = 4)#select four axes
scatter(pcoBB)

pcoBB$li
sum(pcoBB$eig[1:4]) / sum(pcoBB$eig)#0.74
sum(pcoBB$eig[1:3]) / sum(pcoBB$eig)#0.61
sum(pcoBB$eig[1:2]) / sum(pcoBB$eig)#0.43

# check correlations among axes and traits
str(t2)
cor(pcoBB$li, t2, use = "complete.obs")

# due to number of bee species at each site (one site has only 3 bee species),
# we can only use first two axes of PCoA for functional diversity metrics
t.ax <- as.matrix(pcoBB$li[1:2])

# returns pairwise between-site values of each functional beta-diversity component
bb.fun <- functional.beta.pair(a, t.ax, index.family = "sorensen")
str(bb.fun)

# observed functional alpha diversity
# run rao function first
bb.rao <- Rao(sample = t(a), dfunc = tdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
bb.falpha <- bb.rao$FD$Alpha
bb.falpha

library(ape)
library(picante)

# create a copy of the abundance matrix to randomize
ra <- a

# Null model ----

# run the null model with 999 iterations
numberReps <- 999

# create empty matrices to store the results of each iteration of the null model:
# cwms
#nbs <- nrhw <- nrew <- nral <- nrtl <- nrfl <- nrll <- ndisp_0 <- ndisp_1 <- ndisp_2 <- norigin_1 <- norigin_2  <- matrix(NA,
#              nrow = nrow(a), ncol = numberReps, dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# taxonomic beta diversity
ntmatrix <- matrix(NA, nrow = nrow(ra), ncol = numberReps, 
              dimnames = list(rownames(ra), paste0("n", 1:numberReps)))

# functional alpha and beta diversity
#nfalpha <- nfmatrix <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
#              dimnames = list(rownames(a), paste0("n", 1:numberReps)))

#create null model for each repetition:

for(i in 1:numberReps){
  print(i) 
  
  # randomize trait matrix
  #ntraits <- t3[sample(1:nrow(t3)),]
  #rownames(ntraits) <- rownames(t3)
  
  # randomize abundance matrix
  # independent swap constrains by species richness and frequency
  spGB <- randomizeMatrix(samp = ra, null.model = "independentswap")
  #print(rownames(ntraits) == colnames(spBB)) 
  
  # randomize trait distance matrix
  #ntdis <- gawdis(ntraits, w.type = "optimized", opti.maxiter = 300,
  #                groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
  
  # CWM calculations
  #cwm.null <- functcomp(x = ntraits, a = as.matrix(spBB), CWM.type = "all")
  #nqbl_var[,i] <- cwm.null$qbl_var
  #nmbl_var[,i] <- cwm.null$mbl_var
  #nwbl_var[,i] <- cwm.null$wbl_var
  #ntl_0[,i] <- cwm.null$tl_0
  #ntl_1[,i] <- cwm.null$tl_1
  #ntl_2[,i] <- cwm.null$tl_2
  #nnestl_0[,i] <- cwm.null$nestl_0
  #nnestl_2[,i] <- cwm.null$nestl_2
  #nit[,i] <- cwm.null$it
  #nwingl[,i] <- cwm.null$rwingl
  #nheadw[,i] <- cwm.null$rheadw
  #neyel[,i] <- cwm.null$reyel
  #nthairl[,i] <- cwm.null$thairl
  #nbsetael[,i] <- cwm.null$bsetael
  #nbtl[,i] <- cwm.null$rbtl
  
  # Functional alpha diversity
  #nrao <- Rao(sample = t(spBB), dfunc = ntdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
  #nfalpha[,i] <- nrao$FD$Alpha
  
  # Taxonomic beta diversity indices
  nt.dist <- vegdist(spGB, method = "bray")
  nt.dist.matrix <- as.matrix(nt.dist)
  ntmatrix[,i] <- colMeans(nt.dist.matrix)
  
  # Functional beta diversity indices
  #npcoBB <- dudi.pco(sqrt(ntdis), scannf = FALSE, nf = 2)
  #nt <- as.matrix(npcoBB$li)
  #nbb.fun <- functional.beta.pair(spBB, nt, index.family = "sorensen")
  #nfsim.dist <- as.matrix(nbb.fun$funct.beta.sim)
  #nfsne.dist <- as.matrix(nbb.fun$funct.beta.sne)
  #nfsor.dist <- as.matrix(nbb.fun$funct.beta.sor)
  #nfsim[,i] <- colMeans(nfsim.dist)
  #nfsne[,i] <- colMeans(nfsne.dist)
  #nfsor[,i] <- colMeans(nfsor.dist)
}

#_______________________________________________________________________________

write.csv(ntmatrix, file = "tax_betadiv.csv")

write.csv(nqbl_var, file = "Urban_to_Local_2Axes/nqbl_var.u.csv")
write.csv(nmbl_var, file = "Urban_to_Local_2Axes/nmbl_var.u.csv")
write.csv(nwbl_var, file = "Urban_to_Local_2Axes/nwbl_var.u.csv")
write.csv(ntl_0, file = "Urban_to_Local_2Axes/ntl_0.u.csv")
write.csv(ntl_1, file = "Urban_to_Local_2Axes/ntl_1.u.csv")
write.csv(ntl_2, file = "Urban_to_Local_2Axes/ntl_2.u.csv")
write.csv(nnestl_0, file = "Urban_to_Local_2Axes/nnestl_0.u.csv")
write.csv(nnestl_2, file = "Urban_to_Local_2Axes/nnestl_2.u.csv")
write.csv(nit, file = "Urban_to_Local_2Axes/nit.u.csv")
write.csv(nwingl, file = "Urban_to_Local_2Axes/nwingl.u.csv")
write.csv(nheadw, file = "Urban_to_Local_2Axes/nheadw.u.csv")
write.csv(neyel, file = "Urban_to_Local_2Axes/neyel.u.csv")
write.csv(nthairl, file = "Urban_to_Local_2Axes/nthairl.u.csv")
write.csv(nbsetael, file = "Urban_to_Local_2Axes/nbsetael.u.csv")
write.csv(nbtl, file = "Urban_to_Local_2Axes/nbtl.u.csv")

write.csv(nfalpha, file = "Urban_to_Local_2Axes/bb_falpha.u.csv")
write.csv(nfsim, file = "Urban_to_Local_2Axes/bb_fbeta_sim.u.csv")
write.csv(nfsne, file = "Urban_to_Local_2Axes/bb_fbeta_sne.u.csv")
write.csv(nfsor, file = "Urban_to_Local_2Axes/bb_fbeta_sor.u.csv")

write.csv(npsim, file = "Urban_to_Local_2Axes/bb_pbeta_sim.csv")
write.csv(npsne, file = "Urban_to_Local_2Axes/bb_pbeta_sne.csv")
write.csv(npsor, file = "Urban_to_Local_2Axes/bb_pbeta_sor.csv")

#_______________________________________________________________________________
# load the datasets

ntmatrix <- read.csv("tax_betadiv.csv", row.names = 1)

nqbl_var <- read.csv("Urban_to_Local_2Axes/nqbl_var.u.csv", row.names=1)
nmbl_var <- read.csv("Urban_to_Local_2Axes/nmbl_var.u.csv", row.names=1)
nwbl_var <- read.csv("Urban_to_Local_2Axes/nwbl_var.u.csv", row.names=1)
ntl_0 <- read.csv("Urban_to_Local_2Axes/ntl_0.u.csv", row.names=1)
ntl_1 <- read.csv("Urban_to_Local_2Axes/ntl_1.u.csv", row.names=1)
ntl_2 <- read.csv("Urban_to_Local_2Axes/ntl_2.u.csv", row.names=1)
nnestl_0 <- read.csv("Urban_to_Local_2Axes/nnestl_0.u.csv", row.names=1)
nnestl_2 <- read.csv("Urban_to_Local_2Axes/nnestl_2.u.csv", row.names=1)
nit <- read.csv("Urban_to_Local_2Axes/nit.u.csv", row.names=1)
nwingl <- read.csv("Urban_to_Local_2Axes/nwingl.u.csv", row.names=1)
nheadw <- read.csv("Urban_to_Local_2Axes/nheadw.u.csv", row.names=1)
neyel <- read.csv("Urban_to_Local_2Axes/neyel.u.csv", row.names=1)
nthairl <- read.csv("Urban_to_Local_2Axes/nthairl.u.csv", row.names=1)
nbsetael <- read.csv("Urban_to_Local_2Axes/nbsetael.u.csv", row.names=1)
nbtl <- read.csv("Urban_to_Local_2Axes/nbtl.u.csv", row.names=1)

nbsim <- read.csv("Urban_to_Local_2Axes/bb_tbeta_sim.u.csv", row.names=1)
nbsne <- read.csv("Urban_to_Local_2Axes/bb_tbeta_sne.u.csv", row.names=1)
nbsor <- read.csv("Urban_to_Local_2Axes/bb_tbeta_sor.u.csv", row.names=1)

nfalpha <- read.csv("Urban_to_Local_2Axes/bb_falpha.u.csv", row.names=1)
nfsim <- read.csv("Urban_to_Local_2Axes/bb_fbeta_sim.u.csv", row.names=1)
nfsne <- read.csv("Urban_to_Local_2Axes/bb_fbeta_sne.u.csv", row.names=1)
nfsor <- read.csv("Urban_to_Local_2Axes/bb_fbeta_sor.u.csv", row.names=1)

npsim<- read.csv(file = "Urban_to_Local_2Axes/bb_pbeta_sim.csv", row.names=1)
npsne<- read.csv(file = "Urban_to_Local_2Axes/bb_pbeta_sne.csv", row.names=1)
npsor<- read.csv(file = "Urban_to_Local_2Axes/bb_pbeta_sor.csv", row.names=1)

#_______________________________________________________________________________

# SES Calculations ----

#calculate standardized effect sizes (SES) for each trait and index
#the effect size is the difference between the observed value and the expected one
#then divide the effect size by the standard deviation of the null distribution to get the standardized effect size
#allows comparison among sites with different numbers of species

## calculate SES values for each metric

## community weighted means ----
## queen body length variance
SES_qblv <- (cwm.obs$qbl_var - apply(nqbl_var, MARGIN = 1, mean)) / apply(nqbl_var, MARGIN = 1, sd, na.rm=T)
SES_qblv

## taxonomic beta diversity ----

beta.bray <- as.matrix(dis.matrix.t)
beta.bray <- colMeans(beta.bray)
beta.t <- data.frame(beta.bray)

# beta.div <- data.frame(beta.sor, beta.sim, beta.sne)
# may want to combine all beta diversity into one matrix?

## taxonomic diversity - beta bray
SES_bdiv <- (beta.t$beta.bray - apply(ntmatrix, MARGIN = 1, mean)) / apply(ntmatrix, MARGIN = 1, sd, na.rm = T)
SES_bdiv <- as.data.frame(SES_bdiv)
colnames(SES_bdiv) <- c("tbeta")

SES_bdiv$trmt <- as.factor(df3$Trmt_Code)
SES_bdiv
str(SES_bdiv)

boxplot(tbeta ~ trmt, data = SES_bdiv,
        col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        cex.lab = 1.2, cex.axis = 1, ylim = c(-3,3))
stripchart(tbeta ~ trmt, data = SES_bdiv, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

dotchart(SES_bdiv$tbeta, group = SES_bdiv$trmt)
hist(SES_bdiv$tbeta)
boxplot(SES_bdiv$tbeta ~ SES_bdiv$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
plot(SES_bdiv$tbeta, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

tbeta.mod <- lm(tbeta ~ trmt, data = SES_bdiv)
summary(tbeta.mod)
plot(tbeta.mod)
qqnorm(resid(tbeta.mod))
qqline(resid(tbeta.mod))
Anova(tbeta.mod, type = "II")
emmeans(tbeta.mod, pairwise ~ trmt)

## pull out data for each treatment
mp <- SES_bdiv[which(SES_bdiv$trmt == "MP"),]
of <- SES_bdiv[which(SES_bdiv$trmt == "OF"),]
pp <- SES_bdiv[which(SES_bdiv$trmt == "PP"),]
vl <- SES_bdiv[which(SES_bdiv$trmt == "VL"),]

## compare to null expectations by treatment
tbeta.mp <- wilcox.test(mp$tbeta, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
tbeta.mp

tbeta.of <- wilcox.test(of$tbeta, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
tbeta.of

tbeta.pp <- wilcox.test(pp$tbeta, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
tbeta.pp

tbeta.vl <- wilcox.test(vl$tbeta, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
tbeta.vl



