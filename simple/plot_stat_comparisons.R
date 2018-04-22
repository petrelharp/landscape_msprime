png('sites_vs_branches.png', width=5, height=6, units='in', res=300)

par(tck=-.01, bty="l")
layout(matrix(1:4, nrow=2))
for (a in c("f4", "Y3")) {
    x <- read.table(paste0("simple.",a,".tsv"), header=TRUE);
    xlow <- read.table(paste0("simple-lowmut.",a,".tsv"), header=TRUE);
    plot(x, main=paste(gsub('f', 'F', a),': high mutation'), ylim=range(c(x[[2]], xlow[[2]])), ylab='site', xlab='branch') #expression(F[4]^(site)));
    abline(0,1e-8, col='red')
    # points(xlow, main=a, pch=19, col='darkgrey');
    # abline(0,1e-10, col='blue')
    plot(xlow, main=paste(gsub('f', 'F', a),': low mutation'), pch=19, col='darkgrey', ylab='site', xlab='branch');
    abline(0,1e-10, col='blue')
}

dev.off()

