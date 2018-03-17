library(ncdf4)
library(sp)
library(rgdal)
library(maptools)
library(gstat)
library(raster)
library(gdal)


big.data<-nc_open("abiesalba_h1x1.cdf", write=FALSE, readunlim=TRUE, verbose=FALSE, auto_GMT=TRUE, suppress_dimvals=FALSE )

lon <- ncvar_get(big.data,"LONGITUDE")
lat <- ncvar_get(big.data,"LATITUDE")
time <- ncvar_get(big.data,"TIME")
temp <- ncvar_get(big.data,"Temp")
area <- ncvar_get(big.data,"area")
lu.area <- ncvar_get(big.data,"lu_area")


fpc.grid<-ncvar_get(big.data,"fpc_grid")

#lu.area is landuse: this tells how big percentage of the grid is occupied by seas.

tree.num <- lu.area
for (t in 1:2204) {
  tree.num[,,t]<-(1-lu.area[,,t])*area*fpc.grid[,,10,t]/40
}
tree.num[is.na(tree.num)] <- 0


# tree.num is the number of trees on each grid cell at each time point (70x55x2204)


