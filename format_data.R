library(tidyverse)
library(sp)
library(spdep)

load("covid.RData")

States <- raster::getData("GADM", country = "United States", level = 1)
States <- States[States$NAME_1 != "Alaska" & States$NAME_1 != "Hawaii",]

# plot(States)

States$Covid.2020 <- dat$positive
States$AirPol <- dat$AirPol
States$Obesity <- dat$Obesity
States$Smoking <- dat$Smoking
States$Excessive_Drinking <- dat$Excessive_Drinking
States$Popdensity <- dat$Popdensity
States$Uninsured <- dat$Uninsured

neighborhood <- poly2nb(pl = States)

