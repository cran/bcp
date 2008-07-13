% --- Source file: print.bcp.Rd ---
\name{print.bcp}

\alias{print.bcp}

\title{Summarizing Bayesian change point results}

\description{
 print method for class \code{bcp}.
}

\usage{
 print.bcp(x, digits = max(3, .Options$digits - 3), ...)
}

\arguments{
  \item{x}{the result of a call to \code{bcp()}.}
  \item{digits}{an optional specification of the number of digits displayed in the summary statistics.}
  \item{...}{additional arguments.}
  
}

\details{
The function returns the posterior probability of a change point for each position, the posterior means and standard deviations. These results are modeled after the summary method of the coda package (Plummer \emph{et al.}, 2006). If \code{return.mcmc=TRUE} (i.e., if full MCMC results are returned), \code{bcp} objects can be converted into \code{mcmc} objects to view \code{mcmc} summaries -- see examples below.
}

\author{Chandra Erdman and John W. Emerson}


\seealso{\code{\link{bcp}}, \code{\link{summary.bcp}}, and \code{\link{plot.bcp}}.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
  bcp.0 <- bcp(testdata)
  print.bcp(bcp.0)
  plot.bcp(bcp.0)
  
  ##### An MCMC summary from the ``coda'' package #####
   if(require("coda")) {
       bcp.0 <- bcp(testdata, return.mcmc=TRUE)
       BURNIN <- bcp.0$burnin
       M <- bcp.0$burnin + bcp.0$mcmc
       bcp.mcmc <- as.mcmc(bcp.0$mcmc.means)
       summary(bcp.mcmc)
       heidel.diag(bcp.mcmc) # an example convergence diagnostic from the coda package.
   }
  
}


\keyword{datasets}
