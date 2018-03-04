library(raster)
library(landsim)
library(Matrix)

full_habitat <- raster("nussear.grd")
habitat <- aggregate(full_habitat, fact=8)

# a small grid

M <- migration_matrix(habitat,
                      kern="gaussian",
                      discretize=TRUE, disc.fact=50,
                      sigma=5e3, radius=sqrt(2)*max(res(habitat)), 
                      normalize=1)
M[M < 1e-3] <- 0.0
M@x <- M@x/(rowSums(M)[M@i+1L])

Mijx <- data.frame(i=M@i,
                   j=landsim:::p.to.j(M@p),
                   x=M@x)

write.table(Mijx, file="nussear_8x.migr.tsv", row.names=FALSE)

# a big one

M <- migration_matrix(full_habitat,
                      kern="gaussian",
                      discretize=TRUE, disc.fact=50,
                      sigma=5e3, radius=sqrt(2)*max(res(habitat)), 
                      normalize=1)
M[M < 1e-3] <- 0.0
M@x <- M@x/(rowSums(M)[M@i+1L])

Mijx <- data.frame(i=M@i,
                   j=landsim:::p.to.j(M@p),
                   x=M@x)

write.table(Mijx, file="nussear.migr.tsv", row.names=FALSE)
