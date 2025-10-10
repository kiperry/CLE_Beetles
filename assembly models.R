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
colSums(df3[4:94])

# need to remove columns with zero values
# removes species not collected in these sites during this year
library(tidyverse)
df4 <- df3 %>%
  select(where(~any(. != 0)))
colSums(df4[4:62])

