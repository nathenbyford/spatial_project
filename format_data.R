library(tidyverse)
library(sp)
library(spdep)

load("covid.RData")

States <- raster::getData("GADM", country = "United States", level = 1)
States <- States[States$NAME_1 != "Alaska" & States$NAME_1 != "Hawaii",]

# plot(States)

unemployment <- read_csv("unemployment.csv")

state_names <- c(States$NAME_1[8], States$NAME_1[-8])

unemployment <- unemployment |> filter(Year == 2022, Period == "M04") |> 
  mutate(
    State = state_names
  ) |> 
  arrange(State)

dat <- arrange(dat, State)

dat$State == States$NAME_1

# dat <- rename(dat, "NAME_1" = "State")
# 
# States <- left_join(States, dat)

States$Covid.2020 <- dat$positive
States$AirPol <- dat$AirPol
States$Obesity <- dat$Obesity
States$Smoking <- dat$Smoking
States$Excessive_Drinking <- dat$Excessive_Drinking
States$Popdensity <- dat$Popdensity
States$Uninsured <- dat$Uninsured
States$Population <- dat$Pop
States$Unemployment <- unemployment$Value

neighborhood <- poly2nb(pl = States)

save(States, neighborhood, file = "spatial.Rdata")
