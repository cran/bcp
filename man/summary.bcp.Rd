% --- Source file: summary.bcp.Rd ---
\name{summary.bcp}

\alias{summary.bcp}

\title{Summarizing Bayesian change point results}

\description{
 summary method for class ``bcp''.
}

\usage{
 summary.bcp(object, digits = max(3, .Options$digits - 3), ...)
}

\arguments{
  \item{object}{the result of a call to `bcp'.}
  \item{digits}{the number of digits displayed in the summary statistics.}
  \item{...}{additional arguments.}
  
  
}

\details{
The function returns the posterior probability of a change point for each position, the posterior means and standard deviations.  These results are modeled after the `summary.mcmc' output of the coda package (Plummer \emph{et al.}, 2006).  If return.mcmc=TRUE (i.e., if full mcmc results are returned), `bcp' objects can be converted into `mcmc' objects to view `mcmc' summaries -- see examples below.
}

\author{Chandra Erdman}


\seealso{`bcp', `print.bcp', and `plot.bcp'.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(20), rnorm(10, 5, 1), rnorm(20))
  bcp.0 <- bcp(testdata)
  summary.bcp(bcp.0)
  plot.bcp(bcp.0)
  
  ##### An MCMC summary from the ``coda'' package #####
   if(require("coda")) {
     if(is.na(bcp.0$mcmc.means)==FALSE) {
       bcp.mcmc <- as.mcmc(bcp.0$mcmc.means)
       summary(bcp.mcmc)
    }
  }
  
  ##### Coriell chromosome 11 #####
  data(coriell)
  chrom11 <- as.vector(na.omit(coriell$Coriell.05296[coriell$Chromosome==11]))
  bcp.11 <- bcp(chrom11)
  summary.bcp(bcp.11)
  plot.bcp(bcp.11)
  
  # to see bcp and Circular Binary Segmentation results run:
  if(require("DNAcopy")) {
  bcp.11$posterior.prob[length(bcp.11$posterior.brob)] <- 0
  n <- length(chrom11)
  cbs <- segment(CNA(chrom11, rep(1, n), 1:n), verbose = 0)
  cbs.ests <- rep(unlist(cbs$output[6]), unlist(cbs$output[5]))
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  plot(1:n, bcp.11$posterior.mean, type="l", xlab="Location", ylab="Posterior Mean", main="Posterior Means")
  lines(cbs.ests, col="red")
  points(chrom11)
  plot(1:n, bcp.11$posterior.prob, type="l", ylim=c(0,1), xlab="Location", ylab="Posterior Probability of a Change", main="Change Point Locations")
  for(i in 1:(dim(cbs$output)[1]-1)) abline(v=cbs$output$loc.end[i], col="red")
  par(op)
  } else {
	cat("DNAcopy is not loaded")
  }
  
  ##### RealInt #####
  data("RealInt")
  bcp.ri <- bcp(as.vector(RealInt))
  summary.bcp(bcp.ri)
  plot.bcp(bcp.ri)
  
  # to see bcp and Bai and Perron results run:
  if(require("strucchange")) {
  bcp.ri$posterior.prob[length(bcp.ri$posterior.brob)] <- 0
  bp <- breakpoints(RealInt ~ 1, h = 2)$breakpoints
  rho <- rep(0, length(RealInt))
  rho[bp] <- 1
  b.num<-1 + c(0,cumsum(rho[1:(length(rho)-1)]))
  bp.mean <- unlist(lapply(split(RealInt,b.num),mean))
  bp.ri <- rep(0,length(RealInt))
  for(i in 1:length(bp.ri)) bp.ri[i] <- bp.mean[b.num[i]]
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  xax <- seq(1961, 1987, length=103)
  plot(xax, bcp.ri$posterior.mean, type="l", xlab="Time", ylab="Mean", main="Posterior Means")
  lines(xax, bp.ri, col="blue")
  points(RealInt)
  plot(xax, bcp.ri$posterior.prob, type="l", ylim=c(0,1), xlab="Time", ylab="Posterior Probability", main="Posterior Probability of a Change")
  for(i in 1:length(bp.ri)) abline(v=xax[bp[i]], col="blue")
  par(op)
  } else {
  cat("strucchange is not loaded")
  }

}

\references{
Bai J., Perron P. (2003), Computation and Analysis of Multiple Structural Change Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22. url: \url{http:
//qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.

Daniel Barry and J. A. Hartigan (1993), A Bayesian Analysis for Change Point Problems, \emph{Journal of The American Statistical Association}, \bold{88}, 309-19.

Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004), Circular binary segmentation for the analysis of array-based DNA copy number data, \emph{Biostatistics}, \bold{5}, 557-572. url: \url{http://www.bioconductor.org/repository/release1.5/package/html/DNAcopy.html}.

Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines (2006), \emph{The coda Package, version 0.10-7}, CRAN: The Comprehensive R Network. 

Snijders \emph{et al.} (2001), Assembly of microarrays for genome-wide measurement of DNA copy number, \emph{Nature Genetics}, \bold{29}, 263-264. 

Achim Zeileis, Friedrich Leisch, Bruce Hansen, Kurt Hornik, Christian Kleiber (2006), \emph{The strucchange Package, version 1.3-1}, CRAN: The Comprehensive R Network. 
}

\keyword{datasets}
