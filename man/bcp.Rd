% --- Source file: bcp.Rd ---
\name{bcp}

\alias{bcp}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{A Package for Performing a Bayesian Analysis of Change Point Problems}

\description{
 \code{bcp()} implements the Barry and Hartigan (1993) product partition model for the standard change point problem using Markov Chain Monte Carlo. Versions < 2.0 are quadratic in speed, and perform the default 550 iterations in approximately 0.75 seconds for a sequence of length 100.  Versions => 2.0 are linear in speed and partition a sequence of length 10,000 in approximately 45 seconds (compared with 45 minutes for versions < 2.0).  These times were computed on a PC with Windows XP, a Pentium D Processor (2.99 GHz) and 3.50GB of RAM. 
}

\usage{
 bcp(x, w0 = 0.2, p0 = 0.2, burnin = 50, mcmc = 500,
     return.mcmc = FALSE, nwssleigh = NULL)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{a vector of numerical data. }
  \item{w0}{an optional numeric value specifying the prior variance. If  
   no value is specified, the default value of 0.2 is used, as    
   recommended by Barry and Hartigan.  }
  \item{p0}{an optional numeric value specifying the prior probability of a change point at each location in the sequence. If no value is specified, the default value of 0.2 is used, as recommended by Barry and Hartigan. }
  \item{burnin}{the number of burnin iterations.}
  \item{mcmc}{the number of iterations after burnin. }
  \item{return.mcmc}{if \code{return.mcmc=TRUE} the posterior means and the partitions after each iteration are returned.}
  \item{nwssleigh}{the NWS sleigh (see details).}
}

\details{
This algorithm is used when there exists an unknown partition of a sequence into 
contiguous blocks such that the mean is constant within blocks.

When NetWorkSpaces (package \code{nws}) is used for parallel Markov chains, there is network
overhead as well as the overhead associated with burning in each chain.  A special note
is needed about the values returned (with respect to \code{blocks} and, if 
\code{return.mcmc=TRUE}, about \code{mcmc.means} and \code{mcmc.rhos}): the burnins
are collected at the beginning of the vector (or matrix, as the case may be) of results,
followed by blocks of the mcmc results.  So, for example, positions \code{1:burnin}
contain the burnins from the first worker.  And if there are \code{k} workers, then
positions \code{k*burnin+1} through \code{k*burnin+mcmc/k} will contain the results
of the first chain (ideally \code{mcmc} should be divisible by \code{k}).  As a result, convergence diagnostics (for example, the \code{heidel.diag()} function of the coda package) should be applied to parallel \code{bcp} objects with care.


If you are unfamiliar with NetWorkSpaces, it requires both the \code{nws} package
(available from CRAN) and a NetWorkSpaces server (freely available but somewhat more
work to get started).  We will keep a page updated with links and information to help
you get started: \url{http://www.stat.yale.edu/~jay/nws/}.
}

\value{
  \code{bcp()} returns a list containing the following components:
  \item{data}{a vector copy of the data.}
  \item{return.mcmc}{if \code{return.mcmc=TRUE}, the posterior means and the partitions after each iteration are returned.}
  \item{mcmc.means}{if \code{return.mcmc=TRUE}, \code{mcmc.means} contains the posterior means at for each iteration conditional on the state of the partition.}
  \item{mcmc.rhos}{if \code{return.mcmc=TRUE}, \code{mcmc.rhos} contains the partitions after each iteration. A value of 1 indicates the end of a block.}
  \item{blocks}{a vector of the number of blocks after each iteration.}
  \item{posterior.mean}{a vector of the posterior means over the iterations, excuding burnin.}
  \item{posterior.var}{a vector of the posterior variances over the iterations, excuding burnin.}
  \item{posterior.prob}{a vector of the posterior probability of changes over the iterations, excuding burnin.}
  \item{burnin}{the number of burnin iterations.}
  \item{mcmc}{the number of iterations after burnin.}
  \item{w0}{an optional numeric value specifying the prior variance. If no value is specified, the default value of 0.2 is used, as recommended by Barry and Hartigan.}
  \item{p0}{an optional numeric value specifying the prior probability of a change point at each location in the sequence. If no value is specified, the default value of 0.2 is used, as recommended by Barry and Hartigan.}

  The generic accessor functions \code{data}, \code{blocks}, and \code{posterior.mean}, \code{posterior.prob}, \code{burnin}, \code{mcmc}, \code{w0}, and \code{p0} extract various useful components of the list returned by \code{bcp()}.  The functions \code{summary.bcp()}, \code{print.bcp()}, and \code{plot.bcp()} are used to obtain summaries of the results.
}

\author{Chandra Erdman and John W. Emerson}

\note{The vector of data must contain no missing or \code{NA} values.  The user may set \code{.Random.seed} to control the MCMC if desired.}

\seealso{\code{\link{plot.bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for summaries of the results.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
  bcp.0 <- bcp(testdata)
  plot.bcp(bcp.0)

  \dontrun{
    ##### An example using NetWorkSpaces for parallel MCMC; note that
    ##### you must have package nws as well as an NWS server available,
    ##### specified in the nwsHost argument (see Details, above).
    if (require(nws)) {
      s <- sleigh(workerCount=4, nwsHost='hostname.XXX.YYY.ZZZ')
      bcp.parallel <- bcp(testdata, nwssleigh=s)
      stopSleigh(s)
    }
  }
  
  ##### Coriell chromosome 11 #####
  data(coriell)
  chrom11 <- as.vector(na.omit(coriell$Coriell.05296[coriell$Chromosome==11]))
  bcp.11 <- bcp(chrom11)
  plot.bcp(bcp.11)
  
  # to see bcp and Circular Binary Segmentation results run:
  if(require("DNAcopy")) {
  n <- length(chrom11)
  cbs <- segment(CNA(chrom11, rep(1, n), 1:n), verbose = 0)
  cbs.ests <- rep(unlist(cbs$output[6]), unlist(cbs$output[5]))
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  op2 <- par(mar=c(0,4,4,2),xaxt="n", cex.axis=0.75)
  plot(1:n, bcp.11$data, col="grey", pch=20, xlab="Location", ylab="Posterior Mean", main="Posterior Means and Probabilities of a Change")
  lines(cbs.ests, col="red")
  lines(bcp.11$posterior.mean, lwd=2)
  par(op2)
  op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
  plot(1:n, bcp.11$posterior.prob, type="l", ylim=c(0,1), xlab="Location", ylab="Posterior Probability", main="")
  for(i in 1:(dim(cbs$output)[1]-1)) abline(v=cbs$output$loc.end[i], col="red")
  par(op3)
  par(op)
  } else {
        cat("DNAcopy is not loaded")
  }
  
  ##### RealInt #####
  data("RealInt")
  bcp.ri <- bcp(as.vector(RealInt))
  plot.bcp(bcp.ri)
  
  # to see bcp and Bai and Perron results run:
  if(require("strucchange")) {
  bp <- breakpoints(RealInt ~ 1, h = 2)$breakpoints
  rho <- rep(0, length(RealInt))
  rho[bp] <- 1
  b.num<-1 + c(0,cumsum(rho[1:(length(rho)-1)]))
  bp.mean <- unlist(lapply(split(RealInt,b.num),mean))
  bp.ri <- rep(0,length(RealInt))
  for(i in 1:length(bp.ri)) bp.ri[i] <- bp.mean[b.num[i]]
  xax <- seq(1961, 1987, length=103)
  op<-par(mfrow=c(2,1),col.lab="black",col.main="black")
  op2 <- par(mar=c(0,4,4,2),xaxt="n", cex.axis=0.75)
  plot(1:length(bcp.ri$data), bcp.ri$data, col="grey", pch=20, xlab="", ylab="Posterior Mean", main="U.S. Ex-Post Interest Rate")
  lines(bcp.ri$posterior.mean, lwd=2)
  lines(bp.ri, col="blue")
  par(op2)
  op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
  plot(xax, bcp.ri$posterior.prob, yaxt="n", type="l", ylim=c(0,1), xlab="Year", ylab="Posterior Probability", main="")
  for(i in 1:length(bp.ri)) abline(v=xax[bp[i]], col="blue")
  axis(2, yaxp=c(0, 0.9, 3))
  par(op3)
  par(op)
  } else {
    cat("strucchange is not loaded")
  }
}

\references{
J. Bai and P. Perron (2003), Computation and Analysis of Multiple Structural Change Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22. \url{http://qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.

Daniel Barry and J. A. Hartigan (1993), A Bayesian Analysis for Change Point Problems, \emph{Journal of The American Statistical Association}, \bold{88}, 309-19.

Chandra Erdman and John W. Emerson (2008), A Fast Bayesian Change Point Analysis for the Segmentation of Microarray Data, \emph{Bioinformatics}, 24(19), 2143--2148. \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btn404}.

Chandra Erdman and John W. Emerson (2007), bcp: An R Package for Performing a Bayesian Analysis of Change Point Problems. \emph{Journal of Statistical Software}, 23(3), 1--13. \url{http://www.jstatsoft.org/v23/i03/}.

A. B. Olshen, E. S. Venkatraman, R. Lucito, M. Wigler (2004), Circular binary segmentation for the analysis of array-based DNA copy number data, \emph{Biostatistics}, \bold{5}, 557-572.  \url{http://www.bioconductor.org/repository/release1.5/package/html/DNAcopy.html}.

Snijders \emph{et al.} (2001), Assembly of microarrays for genome-wide measurement of DNA copy number, \emph{Nature Genetics}, \bold{29}, 263-264. 

Achim Zeileis, Friedrich Leisch, Kurt Hornik, Christian Kleiber (2002), strucchange: An R Package for Testing for Structural Change in Linear Regression Models, \emph{Journal of Statistical Software}, \bold{7}(2), 1--38. \url{http://www.jstatsoft.org/v07/i02/}. 
}

\keyword{datasets}% at least one, from doc/KEYWORDS
