#####
# Pre-compute the discrete migration matrices
# *and* assignments of samples to grid squares
#####

library(ncdf4)
library(landsim)
library(raster)

big.data <- nc_open("abiesalba_h1x1.cdf", write=FALSE, readunlim=TRUE, verbose=FALSE, auto_GMT=TRUE, suppress_dimvals=FALSE )

lon <- ncvar_get(big.data,"LONGITUDE")
lat <- ncvar_get(big.data,"LATITUDE")
area <- ncvar_get(big.data,"area")
lu.area <- ncvar_get(big.data,"lu_area")
fpc.grid <- ncvar_get(big.data,"fpc_grid")

#lu.area is landuse: this tells how big percentage of the grid is NOT occupied by seas.

tree.num <- lu.area
for (t in 1:2204) {
  tree.num[,,t] <- lu.area[,,t]*area*fpc.grid[,,10,t]/40
}
tree.num[is.na(tree.num)] <- 0

# Here is python code to read the NETCDF in, recorded for posterity
# import netCDF4
# 
# def read_tree_nums_orig():
#     # this is the original projection
#     trees = netCDF4.Dataset("abiesalba_h1x1.cdf", "r", format="NETCDF4")
#     tree_num = np.zeros(shape=(2204, 55, 70))
#     for t in range(2204):
#         # recall python is zero-indexed so '9' = 10th variable in fpc_grid
#         tree_num[t,:,:] = trees['lu_area'][t,0,:,:] * trees['area'] * trees['fpc_grid'][t,9,:,:]/40.
#     # the result has NAs as -99999
#     tree_num[tree_num < 0] = 0.0
#     return tree_num


ll_map <- raster(extent(min(lon)-0.5, max(lon)+0.5, min(lat)-0.5, max(lat)+0.5), res=1,
                  crs=crs("+proj=longlat"))
values(ll_map) <- as.vector(t(rowMeans(tree.num, dims=2)[70:1,55:1]))

########
## First move to an equal-area projection

# from https://epsg.io/3035
europe_proj4 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"  

## Convert values to density before projecting, then back again
dens_map <- ll_map / area(ll_map)
the_map <- projectRaster(dens_map, crs=crs(europe_proj4))
the_map <- the_map * prod(res(the_map))
non_na <- !is.na(values(the_map))

M <- migration_matrix(the_map, accessible=non_na,
                      kern="cauchy",
                      discretize=TRUE, disc.fact=50,
                      sigma=5e2, radius=sqrt(2)*max(res(the_map)), 
                      normalize=1)

M[M < 1e-4] <- 0.0
M@x <- M@x/(rowSums(M)[M@i+1L])

# the resulting matrix will be zero-indexed (like python)
Mijx <- data.frame(i=M@i,
                   j=landsim:::p.to.j(M@p) - 1L,
                   x=M@x)

stopifnot(nrow(M) == ncol(M) && nrow(M) == sum(non_na))

write.table(Mijx, file="abiesalba_epsg3035.migr.tsv", row.names=FALSE)

# Sample assignments to grid squares
# ... here "cell" will be the index of the cell in the *nonempty* squares, i.e., the rows/columns of M

coords <- read.csv("coordsXY.csv", header=TRUE)
samples <- spTransform(SpatialPoints(coords[,c("X","Y")], proj4string=crs("+proj=longlat")), proj4string(the_map))
coords[,c("X","Y")] <- coordinates(samples)
stopifnot(sum(!is.na(values(the_map))) == nrow(M))
coords$cell <- match(cellFromXY(the_map, samples), which(non_na))
# reorder to match what we'll get from msprime
coords <- coords[order(coords$cell),]
coords$msp_id <- (1:nrow(coords))-1
write.table(coords, "coords_epsg3035.tsv", row.names=FALSE)

## translate the tree numbers
ll_area <- values(area(ll_map))

trans_map <- function (x) {
    values(ll_map) <- t(x[70:1, 55:1]) / ll_area
    out_map <- projectRaster(ll_map, to=the_map)
    return(prod(res(the_map)) * values(out_map))
}

# only write out those locations corresponding to nonempty portions of the map
# (i.e., the rows/columns of M)
trans.tree.num <- apply(tree.num, 3, trans_map)
write.table( t(trans.tree.num[non_na,]), file="tree_num_epsg3035.tsv", col.names=FALSE, row.names=FALSE, na="nan") 


#### UNUSED BELOW HERE
# Here is the same thing in the original projection
if (FALSE) {

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

}
