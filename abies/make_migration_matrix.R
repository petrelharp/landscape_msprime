library(ncdf4)
library(landsim)
library(raster)

big.data<-nc_open("abiesalba_h1x1.cdf", write=FALSE, readunlim=TRUE, verbose=FALSE, auto_GMT=TRUE, suppress_dimvals=FALSE )

lon <- ncvar_get(big.data,"LONGITUDE")
lat <- ncvar_get(big.data,"LATITUDE")
area <- ncvar_get(big.data,"area")
lu.area <- ncvar_get(big.data,"lu_area")
fpc.grid<-ncvar_get(big.data,"fpc_grid")

#lu.area is landuse: this tells how big percentage of the grid is NOT occupied by seas.

tree.num <- lu.area
for (t in 1:2204) {
  tree.num[,,t] <- lu.area[,,t]*area*fpc.grid[,,10,t]/40
}
tree.num[is.na(tree.num)] <- 0


ll_map <- raster(extent(min(lon)-0.5, max(lon)+0.5, min(lat)-0.5, max(lat)+0.5), res=1,
                  crs=crs("+proj=longlat"))
values(ll_map) <- as.vector(rowMeans(tree.num, dims=2)[70:1,55:1])

## produce migration matrices
# for forwards-time
#
# This will produce a migration matrix
# between all adjacent cells that have nonzero population
# at *any* point in history

# ROUGHLY a degree of lat/long is 1e5m
M <- migration_matrix(ll_map,
                      kern="cauchy",
                      discretize=TRUE, disc.fact=50,
                      sigma=5e2 * 1e-5, radius=sqrt(2)*max(res(ll_map)), 
                      normalize=1)

M[M < 1e-4] <- 0.0
M@x <- M@x/(rowSums(M)[M@i+1L])

Mijx <- data.frame(i=M@i,
                   j=landsim:::p.to.j(M@p),
                   x=M@x)

write.table(Mijx, file="abiesalba_longlat.migr.tsv", row.names=FALSE)




## Really we should work with an equal-area projection:

# from https://epsg.io/3035
europe_proj4 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"  

the_map <- projectRaster(ll_map, crs=crs(europe_proj4))


M <- migration_matrix(the_map,
                      kern="cauchy",
                      discretize=TRUE, disc.fact=50,
                      sigma=5e2, radius=sqrt(2)*max(res(the_map)), 
                      normalize=1)

M[M < 1e-4] <- 0.0
M@x <- M@x/(rowSums(M)[M@i+1L])

Mijx <- data.frame(i=M@i,
                   j=landsim:::p.to.j(M@p),
                   x=M@x)

write.table(Mijx, file="abiesalba_epsg3035.migr.tsv", row.names=FALSE)

