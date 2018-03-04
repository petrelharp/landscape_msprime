layout(matrix(1:6, nrow=2))
for (a in c("f4", "f3", "f2", "Y3", "Y2")) {
    x <- read.table(paste0("nussear_8x.",a,".tsv"), header=TRUE); 
    plot(x, main=a); 
    abline(0,1e-8, col='red') 
}

