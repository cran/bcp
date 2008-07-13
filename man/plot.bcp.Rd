% --- Source file: plot.bcp.Rd ---
\name{plot.bcp}

\alias{plot.bcp}

\title{ Plotting Bayesian change point results}

\description{
  \code{plot.bcp()} produces summary plots of the results of \code{bcp()}.
}

\usage{
 plot.bcp(x, ...)
}

\arguments{
  \item{x}{the result of a call to \code{bcp()}.}
  \item{...}{additional arguments.}  
}

\details{
\code{plot.bcp()} produces the following plots:
     
Posterior Means: location in the sequence versus the posterior mean over the iterations. 

Posterior Probability of a Change: location in the sequence versus the relative frequency of iterations which resulted in a change point. 

}

\author{Chandra Erdman and John W. Emerson}


\seealso{\code{\link{bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for complete results and summary statistics.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
  bcp.0 <- bcp(testdata)
  plot.bcp(bcp.0)

}

\keyword{datasets}
