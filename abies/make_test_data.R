#!/usr/bin/env Rscript

# Make a small, fake dataset for testing

set.seed(23)
library(Matrix)

spatial_dims <- c(5, 3)
ntimes <- 7

# array of tree population sizes
ntrees <- array(exp(rnorm(prod(spatial_dims)*ntimes, mean=4, sd=1)),
                dim=c(spatial_dims, ntimes))
for (k in 2:ntimes) {
    ntrees[,,k] <- 0.9 * ntrees[,,k-1] + 0.1 * ntrees[,,k]
}

write.table(t(matrix(ntrees,ncol=ntimes)), file="test_tree_num.tsv", 
            col.names=FALSE, row.names=FALSE, na="nan")


# migration matrix: nearest-neighbor, no diagonals
ij <- matrix(1:prod(spatial_dims), nrow=spatial_dims[1])
adj_rows <- outer(1:prod(spatial_dims), 1:prod(spatial_dims), 
                  function (i, j) abs(row(ij)[i] - row(ij)[j]) == 1)
adj_cols <- outer(1:prod(spatial_dims), 1:prod(spatial_dims), 
                  function (i, j) abs(col(ij)[i] - col(ij)[j]) == 1)
adj <- (adj_rows | adj_cols)
# we need these to be zero-indexed!
Mijx <- data.frame(i=row(adj)[adj]-1L, j=col(adj)[adj]-1L)
Mijx$x <- 1/tapply(Mijx$j, Mijx$i, length)[Mijx$i+1L]
stopifnot(all(tapply(Mijx$x, Mijx$i, sum) == 1))

write.table(Mijx, file="test.migr.tsv", row.names=FALSE)

# sample coordinates: pretend this is [0,1]^2
npops <- 5
popsizes <- 1 + rpois(npops, 3)
pop_x <- spatial_dims[1] * runif(npops)
pop_y <- spatial_dims[2] * runif(npops)
coords <- data.frame(pop=rep(1:npops, popsizes),
                     id=paste("pop", rep(1:npops, popsizes),
                              unlist(lapply(popsizes, seq_len)), sep='_'),
                     msp_id=seq_len(sum(popsizes)))
coords$X <- pop_x[coords$pop]
coords$Y <- pop_y[coords$pop]
coords$cell <- ij[cbind(ceiling(coords$X), ceiling(coords$Y))]
coords <- coords[,c("X", "Y", "pop", "id", "cell", "msp_id")]

write.table(coords, file="coords_test.tsv", row.names=FALSE)
