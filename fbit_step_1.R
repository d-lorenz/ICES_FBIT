

library(rgdal)
source(file = "CSquare.R")
library(marmap)

bedShp <- readOGR(dsn = "seabed ita.shp")
bb <- bbox(bedShp)
bb[,1] <- c(7, 35)
cs <- c(0.05, 0.05)
cc <- bb[, 1] + (cs/2)
cd <- ceiling(diff(t(bb))/cs)
grd <- sp::GridTopology(cellcentre.offset = cc,
                        cellsize = cs,
                        cells.dim = cd)
SpP_grd <- as.SpatialPolygons.GridTopology(grd)
tmp_over2 <- SpP_grd[bedShp,]

# assign c-squares
coord <- coordinates(tmp_over2)
squares <- CSquare(coord[,1], coord[,2], 0.05)
rnames <- sapply(slot(tmp_over2, "polygons"), function(x) slot(x, "ID"))
LOCUNI <- as.data.frame(seq(1, length(tmp_over2)))
rownames(LOCUNI) <- rnames
bargrid <- SpatialPolygonsDataFrame(tmp_over2, LOCUNI)
bargrid@data$csquares <- squares
coord <- SpatialPoints(coords = coordinates(bargrid),
                       proj4string = CRS(proj4string(bedShp)))
tr <- over(coord, bedShp)
bargrid@data$MSFDhab<-tr$MSFD_predo 


# save bargrid
save(bargrid, file="region_grid.RData")


