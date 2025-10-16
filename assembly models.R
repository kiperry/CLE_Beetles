############################################################################
# Testing patterns of community assemlby among habitats
# Vacant lots, pocket prairies, forests, old fields
#
# Cleveland, Ohio and surrounding counties
#
# Ground beetle communities
# Ground beetle communities as indicators
#
# 10 October 2025
############################################################################

# import the data

library(readxl)
df <- read_excel("2019_Ground_Beetles_Raw Data.xlsx",
                 sheet = "Raw",
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

# need to remove columns with zero values
# removes species not collected in these sites during this year
library(tidyverse)
df4 <- df3 %>%
  select(where(~any(. != 0)))
colSums(df4[4:62])


## Taxonomic diversity ----

# load libraries
library(hillR)
library(chemodiv)
library(betapart)
library(vegan)
library(car)
library(emmeans)

# calculate diversity metrics

# Hill numbers
# q = 0 (default) for species richness
# q = 2 for inverse Simpson

# margin 1 is the default, indicates that sites are rows
# if sites are columns, then set margin = 2 (this format is less common)

# total abundance
df4$abund <- rowSums(df4[4:62])
boxplot(df4$abund ~ df4$Trmt_Code)

dotchart(df4$abund, group = df4$Trmt_Code)
hist(df4$abund)
boxplot(df4$abund ~ df4$Trmt_Code)

abund.mod <- glm(abund ~ Trmt_Code, family = poisson, data = df4)
summary(abund.mod)
plot(abund.mod)
qqnorm(resid(abund.mod))
qqline(resid(abund.mod))
Anova(abund.mod, type = "II")
emmeans(abund.mod, pairwise ~ Trmt_Code)

boxplot(abund ~ Trmt_Code, data = df4,
        col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"))
stripchart(abund ~ Trmt_Code, data = df4, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

# species richness
df4$rich <- hill_taxa(df4[4:62], q = 0, MARGIN = 1)
boxplot(df4$rich ~ df4$Trmt_Code)

dotchart(df4$rich, group = df4$Trmt_Code)
hist(df4$rich)
boxplot(df4$rich ~ df4$Trmt_Code)

rich.mod <- glm(rich ~ Trmt_Code, family = poisson, data = df4)
summary(rich.mod)
plot(rich.mod)
qqnorm(resid(rich.mod))
qqline(resid(rich.mod))
Anova(rich.mod, type = "II")
emmeans(rich.mod, pairwise ~ Trmt_Code)

boxplot(rich ~ Trmt_Code, data = df4,
        col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"))
stripchart(rich ~ Trmt_Code, data = df4, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

# species diversity
df4$div <- hill_taxa(df4[4:62], q = 2, MARGIN = 1)
boxplot(df4$div ~ df4$Trmt_Code)

dotchart(df4$div, group = df4$Trmt_Code)
hist(df4$div)
boxplot(df4$div ~ df4$Trmt_Code)

div.mod <- lm(log(div) ~ Trmt_Code, data = df4)
summary(div.mod)
plot(div.mod)
qqnorm(resid(div.mod))
qqline(resid(div.mod))
Anova(div.mod, type = "II")
emmeans(div.mod, pairwise ~ Trmt_Code)

boxplot(div ~ Trmt_Code, data = df4,
        col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"))
stripchart(div ~ Trmt_Code, data = df4, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

# species evenness
df4$eve <- calcDiv(df4[4:62], type = "HillEven", q = 2)
boxplot(df4$eve$HillEven ~ df4$Trmt_Code)

dotchart(df4$eve$HillEven, group = df4$Trmt_Code)
hist(df4$eve$HillEven)
boxplot(df4$eve$HillEven ~ df4$Trmt_Code)

eve.mod <- lm(eve$HillEven ~ Trmt_Code, data = df4)
summary(eve.mod)
plot(eve.mod)
qqnorm(resid(eve.mod))
qqline(resid(eve.mod))
Anova(eve.mod, type = "II")
emmeans(eve.mod, pairwise ~ Trmt_Code)

boxplot(eve$HillEven ~ Trmt_Code, data = df4,
        col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"))
stripchart(eve$HillEven ~ Trmt_Code, data = df4, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

# compute a dissimilarity matrix
# the method option let's you indicate which dissimilarity metric to calculate
# here, we are using abundance data and calculating the bray-curtis dissimilarity matrix
dis.matrix.t <- vegdist(df4[4:62], method = "bray")

# run the nonmetric multidimensional scaling model
nmds.tax <- metaMDS(dis.matrix.t, trymax = 500, autotransform = TRUE, k = 2)
nmds.tax # stress is quality of fit
stressplot(nmds.tax)
plot(nmds.tax) # basic plot with no treatment distinctions

# plot the NMDS model
levels(df4$Trmt_Code)
ordiplot(nmds.tax, disp = "sites", type = "n", xlim = c(-3, 3), ylim = c(-1.5, 1.5))
points(nmds.tax, dis = "sites", select = which(df4$Trmt_Code=="MP"), pch = 17, cex = 2, col = "seagreen4")
points(nmds.tax, dis = "sites", select = which(df4$Trmt_Code=="OF"), pch = 18, cex = 2, col = "lightgoldenrod2")
points(nmds.tax, dis = "sites", select = which(df4$Trmt_Code=="PP"), pch = 15, cex = 2, col = "pink1")
points(nmds.tax, dis = "sites", select = which(df4$Trmt_Code=="VL"), pch = 16, cex = 2, col = "seagreen2")
ordiellipse(nmds.tax, df4$Trmt_Code, draw = "lines", col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)

legend("topleft", legend = c("Forest", "Old Field", "Pocket Prarie", "Vacant Lot"),
       pch = c(17, 18, 15, 16), cex = 1.5, bty = "n", col = c("seagreen4", "lightgoldenrod2", "pink1", "seagreen2"))

## Test for differences in predator composition among treatments

# PERMANOVA tests whether the group centroid of communities differs among groups
# in multivariate space (e.g. different community composition)
adonis2(dis.matrix.t ~ df4$Trmt_Code, permutations = 999)

library(pairwiseAdonis)
pairwise.adonis(dis.matrix.t, df4$Trmt_Code)

# BETADISPER tests whether the dispersion of a group from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogeneity of variances
tax.beta <- betadisper(dis.matrix.t, df4$Trmt_Code, type = c("median"))
anova(tax.beta)
plot(tax.beta)
boxplot(tax.beta, ylab = "Distance to median")
TukeyHSD(tax.beta, which = "group", conf.level = 0.95)

################################################################################
