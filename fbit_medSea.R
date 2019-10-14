
# Adaptation of the workflow @ https://github.com/ices-eg/FBIT/tree/dev/TAF%20-%20ICES%20tutorial

library(rgdal)
library(rje)
library(ggplot2)
library(RColorBrewer)
library(rworldmap)
library(rworldxtra)
library(broom)
source("aux/map_plot.R") # warnings are okay?
source("aux/RBS.R")

#### Step 1 ####

bedShp <- readOGR(dsn = "seabed ita.shp")
bb <- bbox(bedShp)
bb[,1] <- c(7, 35)
cs <- c(0.05, 0.05)
cc <- bb[, 1] + (cs/2)
cd <- ceiling(diff(t(bb))/cs)
grd <- sp::GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
SpP_grd <- as.SpatialPolygons.GridTopology(grd)
tmp_over2 <- SpP_grd[bedShp,]

bargrid <- SpatialPolygonsDataFrame(tmp_over2)
coord <- SpatialPoints(coords = coordinates(bargrid), proj4string = CRS(proj4string(bedShp)))
tr <- over(coord, bedShp)
bargrid@data$MSFDhab <- tr$MSFD_predo 



#### Step 2 ####

load(file = "data/cooBio.rData")

# get longevity categories seperate for each station 
ID        <-rep(cooBio@data$H.Num,3)
MSFD      <-rep(cooBio@data$MSFDhab,3)
Cumb      <-c(cooBio@data$L1,(cooBio@data$L1+cooBio@data$L1_3),(cooBio@data$L1+cooBio@data$L1_3+cooBio@data$L3_10))
Longevity <-c(rep(1,nrow(cooBio@data)),rep(3,nrow(cooBio@data)),rep(10,nrow(cooBio@data)))  
Depth <-rep(cooBio@data$Depth,3)

fulldat   <-data.frame(ID,MSFD,Cumb,Longevity, Depth) 
fulldat$ll <-log(fulldat$Longevity)

# add a small number to values very close to 0 and 1 
for (i in 1:(nrow(fulldat))){
  if (fulldat$Cumb[i] < 1e-3){ fulldat$Cumb[i] <- 1e-3}
  if (fulldat$Cumb[i] > 0.999){fulldat$Cumb[i] <- 0.999}
}

fulldat <- fulldat[-which(is.na(fulldat$MSFD)),]

# fit a linear mixed model with sampling station as random factor and MSFD habitats as exploratory variable
mod1   <-  glmer(Cumb ~ ll + MSFD*ll + Depth + (1 | ID), data=fulldat, family=binomial)
mod2   <-  glmer(Cumb ~ ll + MSFD + Depth + (1 | ID), data=fulldat, family=binomial)
mod3   <-  glmer(Cumb ~ ll + (1 | ID), data=fulldat, family=binomial)
mod4   <-  glmer(Cumb ~ ll + Depth + (1 | ID), data=fulldat, family=binomial)
mod5   <-  glmer(Cumb ~ ll + Depth*ll + (1 | ID), data=fulldat, family=binomial)
mod6   <-  glmer(Cumb ~ ll + MSFD + (1 | ID), data=fulldat, family=binomial)
mod7   <-  glmer(Cumb ~ ll + MSFD*Depth + (1 | ID), data=fulldat, family=binomial)
AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7)

modGlm1   <-  glm(Cumb ~ ll + Depth, data=fulldat, family=binomial)
modGlm2   <-  glm(Cumb ~ ll, data=fulldat, family=binomial)
modGlm3   <-  glm(Cumb ~ ll + MSFD, data=fulldat, family=binomial)
modGlm4   <-  glm(Cumb ~ ll + MSFD + Depth, data=fulldat, family=binomial)

AIC(modGlm1, modGlm2, modGlm3, modGlm4)
# models give a singular fit --> the random effect is very small (but you can argue that it, in principle, has to be included)

modcoeff  <-  fixef(mod4)

save(modcoeff,file="Coefficients_Bdata.RData") 



#### Step 3 ####

load("data/Coefficients_Bdata.RData")  

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

outmap <- map_plot(barClip,"medLong",bluegreen)

# estimate the slope and intercept for each gridcell
slope <- rep(coef_ll,nrow(barClip@data))
intercept <- coef_int + coef_OCM*OCM + coef_OCS*OCS

barClip@data$intercept <- intercept
barClip@data$slope     <- slope



#### Step 4 ####


Depl_TBB  <- 0.14 * barClip@data$TBB_SurfSAR   ### data from Hiddink et al. PNAS 2017 Table S4
Depl_OT   <- 0.06 * barClip@data$OT_SurfSAR    ### data from Hiddink et al. PNAS 2017 Table S4

Depl <- cbind(Depl_TBB, Depl_OT)
Depl_tot <- rowSums(Depl, na.rm = TRUE)
barClip@data$Depl_tot <- Depl_tot

state <- c()
for(j in 1:nrow(barClip)){
  state[j] <- RBS(a = barClip@data$slope[j],
                  b = barClip@data$intercept[j],
                  Fd = barClip@data$Depl_tot[j])
}

barClip@data$state <- state

# plot state as a function of total depletion
plot(barClip@data$state~barClip@data$Depl_tot, 
     ylab="State (PD model)", xlab="Depletion (SAR*d)",las=1)

map_plot(barClip,"state",purples)

aggregate(TBB_SurfSAR ~ MSFDhab,
          data = barClip@data, mean)

# estimate state per MSFD habitat
aggregate(state ~ MSFDhab,
          data = barClip@data, mean)


