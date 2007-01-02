\name{RealInt}
\alias{RealInt}
\title{US Ex-post Real Interest Rate data, 1961(1):1986(3)}
\description{
US ex-post real interest rate: the three-month treasury bill
deflated by the CPI inflation rate.
}
\usage{data("RealInt")}
\format{
A quarterly time series from 1961(1) to 1986(3).
}

\source{The data is available online in the data archive of the
\emph{Journal of Applied Econometrics}. url: 
\url{http://qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.}

\references{
Bai J., Perron P. (2003), Computation and Analysis of Multiple Structural Change Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.

Achim Zeileis, Friedrich Leisch, Bruce Hansen, Kurt Hornik, Christian Kleiber (2006), \emph{The strucchange Package, version 1.3-1}, CRAN: The Comprehensive R Network. 
}

\examples{
  ##### RealInt #####
  data("RealInt")
  bcp.ri <- bcp(as.vector(RealInt))
  summary.bcp(bcp.ri)
  plot.bcp(bcp.ri)
  
  # to see bcp and Bai and Perron results run:
  if(require("strucchange")) {
  bcp.ri$rhos[,ncol(bcp.ri$rhos)] <- 0
  changepoint.freq <- apply(bcp.ri$rhos[bcp.ri$burnin:(bcp.ri$burnin+bcp.ri$burnin),1:dim(bcp.ri$results)[2]],2,mean)
  bp <- breakpoints(RealInt ~ 1, h = 2)$breakpoints
  rho <- rep(0, length(RealInt))
  rho[bp] <- 1
  b.num<-1 + c(0,cumsum(rho[1:(length(rho)-1)]))
  bp.mean <- unlist(lapply(split(RealInt,b.num),mean))
  bp.ri <- rep(0,length(RealInt))
  for(i in 1:length(bp.ri)) bp.ri[i] <- bp.mean[b.num[i]]
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  xax <- seq(1961,1987, length=103)
  plot(xax, bcp.ri$posterior.mean, type="l", xlab="Time", ylab="Mean", main="Posterior Means")
  lines(xax, bp.ri, col="blue")
  points(RealInt)
  plot(xax, changepoint.freq, type="l", xlab="Time", ylab="Posterior Probability", main="Posterior Probability of a Change")
  for(i in 1:length(bp.ri)) abline(v=xax[bp[i]], col="blue")
  par(op)
  } else {
  cat("strucchange is not loaded")
  }

}

\keyword{datasets}
