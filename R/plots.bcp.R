"plots.bcp" <-
function(x, burnin=50, mcmc=500, cbs=FALSE) {
  x$rhos[,ncol(x$rhos)] <- 0
  blocksize.dist <- table(x$blocks[burnin:(burnin+mcmc)])
  changepoint.freq <- apply(x$rhos[burnin:(burnin+mcmc),1:dim(x$results)[2]],2,mean)
  if (cbs) {
    require("DNAcopy")
    cbs <- segment(CNA(x$data, rep(1, length(x$data)), 1:length(x$data)), verbose = 0)
    cbs.ests <- rep(unlist(cbs$output[6]), unlist(cbs$output[5]))
    par(mfrow=c(3,1),col.lab="black",col.main="black")
    plot(1:length(x$posterior.mean), x$posterior.mean, type="l",
         xlab="Location", ylab="Posterior Mean", main="Posterior Means")
    lines(cbs.ests, col="red")
    points(x$data)
    plot(1:length(x$posterior.mean), changepoint.freq, type="l",
         xlab="Location", ylab="Posterior Probability", main="Change Point Locations")
      for(i in 1:(dim(cbs$output)[1]-1)){
        abline(v=cbs$output$loc.end[i], col="red")
      }
    plot(as.numeric(names(blocksize.dist )), blocksize.dist/nrow(x$results),
         xlab="Number of Blocks", ylab="Posterior Probability", main="Distribution of Number of Blocks")
    abline(v=dim(cbs$output)[1], col="red")
  } else {
    par(mfrow=c(3,1),col.lab="black",col.main="black")
    plot(1:length(x$data), x$data, xlab="Location", ylab="Posterior Mean", main="Posterior Means")
    lines(x$posterior.mean)
    plot(1:length(x$posterior.mean), changepoint.freq, type="l",
         xlab="Location", ylab="Posterior Probability", main="Change Point Locations")
    plot(as.numeric(names(blocksize.dist )), blocksize.dist/nrow(x$results),
         xlab="Count", ylab="Posterior Probability", main="Distribution of Number of Blocks")
    }
}

