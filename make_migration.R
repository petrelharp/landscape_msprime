library(raster)
library(landsim)
library(Matrix)

full_habitat <- raster("nussear.grd")
habitat <- aggregate(full_habitat, fact=8)

M <- migration_matrix(habitat,
                      kern="gaussian",
                      discretize=TRUE, disc.fact=50,
                      sigma=5e3, radius=sqrt(2)*max(res(habitat)), 
                      normalize=1)


