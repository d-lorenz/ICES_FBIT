
library(rje)
library(ggplot2)
library(RColorBrewer)
library(rworldmap)
library(rworldxtra)
library(broom)


load("dataset/Coefficients_Bdata.RData")  

coef_int <- modcoeff[1]
coef_ll  <- modcoeff[2]
coef_OCM <- modcoeff[3]
coef_OCS <- modcoeff[4]


# only select region for which we had sampling data in step 2
barClip <- subset(bargrid,bargrid@data$MSFDhab %in% c("Shallow sublittoral mud",
                                                      "Shelf sublittoral mud",
                                                      "Shallow sublittoral sand",
                                                      "Shelf sublittoral sand"))

# predict longevity at a certain location
OCM <- ifelse(barClip@data$MSFDhab %in% c("Shallow sublittoral mud",
                                          "Shelf sublittoral mud"),1,0)
OCS <- ifelse(barClip@data$MSFDhab %in% c("Shallow sublittoral sand",
                                          "Shelf sublittoral sand"),1,0)

medLong <- exp((logit(0.5)- coef_int - coef_OCM*OCM - coef_OCS*OCS) / coef_ll)
barClip@data$medLong <- medLong

source("map_plot.R") # warnings are okay
outmap <- map_plot(barClip,"medLong",bluegreen)

# estimate the slope and intercept for each gridcell
slope <- rep(coef_ll,nrow(barClip@data))
intercept <- coef_int + coef_OCM*OCM + coef_OCS*OCS

barClip@data$intercept <- intercept
barClip@data$slope     <- slope



# save barClip
save(barClip, file="region_grid.RData")

