% --- Source file: plot.bcp.legacy.Rd ---
\name{plot.bcp.legacy}

\alias{legacyplot}

\title{Plotting univariate Bayesian change point results}

\description{
  \code{legacyplot()} produces summary plots of the results of \code{bcp()} when used for univariate analysis; it was the default method prior to package version 3.0.0.
}

\usage{
 legacyplot(x, ...)
}

\arguments{
  \item{x}{the result of a call to \code{bcp()}.}
  \item{...}{additional arguments.}  
}

\details{
\code{legacyplot()} produces the following plots using \code{base} graphics:
     
Posterior Means: location in the sequence versus the posterior mean over the iterations. 

Posterior Probability of a Change: location in the sequence versus the relative frequency of iterations which resulted in a change point. 

}

\author{Chandra Erdman and John W. Emerson}


\seealso{\code{\link{plot.bcp}}, \code{\link{bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for complete results and summary statistics.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
  bcp.0 <- bcp(testdata)
  legacyplot(bcp.0)

}

\keyword{datasets}
