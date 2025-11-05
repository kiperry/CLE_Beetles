#_________________________________________________________________________________
# Testing patterns of community assembly among habitats
# Vacant lots, pocket prairies, forests, old fields
#
# Cleveland, Ohio and surrounding counties
#
# Ground beetle communities: partial least squares canonical analysis (PLSCA)
# using standardized effect sizes (SES) generated from null models
#
# 5 November 2025
#_________________________________________________________________________________

# Import SES and landscape data ----
library(readxl)

SES <- read.csv("SES_CLE.csv", row.names=1)

# load the landscape data set
land <- read_excel("Landscape_data_2000m.xlsx",
                   sheet = "2000_m",
                   range = "A1:M33")
str(land)


# Check for correlated variables ----
## pull out variables to include in the analysis

# 500 m
plot(land[, c(5:13)], pch = 19)

library(ggplot2)
library(GGally)

cp <- ggpairs(land[, c(5:13)], upper = list(continuous = wrap("cor", size = 5, color = "black")))
cp + theme(strip.text.x = element_text(size = 18), strip.text.y = element_text(size = 10))
cor(land[, c(5:13)], method = c("pearson"), use = "complete.obs")

plot(land$PLAND_GS, pch = 19)
plot(land$LPI_GS, pch = 19)
plot(land$ED_GS, pch = 19)
plot(land$PLAND_FOR, pch = 19)
plot(land$ED_FOR, pch = 19)
plot(land$PLAND_IMP, pch = 19)
plot(land$SIDI, pch = 19)
plot(land$ENN, pch = 19)

# PLSCA for landscape (2000 m) data ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics")
library(mixOmics)


trait.ses <- pls(land[, c(5:13)], SES[, c(1:12,15,18)], mode = c("canonical"), 
                 ncomp = 2, scale = TRUE, max.iter = 100)
trait.ses

trait.ses$prop_expl_var
trait.ses$loadings # identify a cutoff level

plotVar(trait.ses)
plotLoadings(trait.ses)

## Reduced PLSCA ----
## with variables that made the initial cutoff (0.3)

trait.ses.red <- pls(land[, c(5,8,11,12)], SES[, c(2:4,6,8,9,11:12,15)], mode = c("canonical"), 
                     ncomp = 2, scale = TRUE, max.iter = 100)
trait.ses.red

trait.ses.red$loadings
trait.ses.red$prop_expl_var

plotVar(trait.ses.red)
plotLoadings(trait.ses.red)

nw.trait.ses.red <- network(trait.ses.red, cutoff = 0.5, color.edge = color.spectral(2), 
                            lty.edge = c("solid", "dashed"), lwd.edge = 2)
nw.trait.ses.red

# Figure ----

# save the loadings as a new data frame to make a better figure
trait.ses.red.pred <- as.data.frame(trait.ses.red$loadings$X)
trait.ses.red.resp <- as.data.frame(trait.ses.red$loadings$Y)

png("Figures/PLSCA.png", width = 1800, height = 1500, pointsize = 30)

par(mfrow=c(1,1))

plot(trait.ses.red.pred$comp2 ~ trait.ses.red.pred$comp1, pch = 19,
     lwd = 1.9, ylim = c(-1, 1), xlim = c(-1, 1), col = "gray60",
     xlab = "PLS Axis 1", ylab = "PLS Axis 2", main = "", cex = 2)
points(trait.ses.red.resp$comp2 ~ trait.ses.red.resp$comp1, pch = 15, cex = 2)
abline(h = 0.0, v = 0.0, col = "black", lwd = 1, lty=1)

# add text for predictor variables (landscape)
text(-0.58, -0.20, "%Forest", pos = 2, font = 1, cex = 1)
text(0.30, -0.91, "%Greenspace", pos = 4, font = 1, cex = 1)
text(0.59, 0.34, "%Impervious", pos = 4, font = 1, cex = 1)
text(0.49, -0.1, "Landscape Diversity", pos = 4, font = 1, cex = 1)

# add text for response variables (traits)
text(0.36, 0.04, "Head width", pos = 2, font = 2, cex = 1)
text(-0.16, -0.687, "Eye width", pos = 2, font = 2, cex = 1)
text(-0.33, -0.33, "Antennae length", pos = 2, font = 2, cex = 1)
text(-0.37, -0.12, "Brachypterous", pos = 2, font = 2, cex = 1)
text(0.34, 0.10, "Macropterous", pos = 2, font = 2, cex = 1)
text(-0.38, 0.42, "Nocturnal", pos = 2, font = 2, cex = 1)
text(0.32, -0.40, "Diurnal", pos = 2, font = 2, cex = 1)
text(-0.34, 0.185, "T Beta", pos = 2, font = 2, cex = 1)
text(-0.36, -0.06, "F Beta", pos = 2, font = 2, cex = 1)


dev.off()




