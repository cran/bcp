% --- Source file: interval.prob.Rd ---
\name{interval.prob}

\alias{interval.prob}

\title{The probability of a change point in a specified interval.}

\description{
From a \code{bcp} object (with \code{return.mcmc=TRUE}), \code{interval.prob()} calculates the probability of a change point in a specified interval. 
}

\usage{
 interval.prob(object, start, end)
}

\arguments{
  \item{object}{the result of a call to \code{bcp()}.}
  \item{start}{the starting index of the interval.}
  \item{end}{the ending index of the interval.}
  
}

\details{
The function returns the posterior probability of at least one change point in the specified interval.
}

\author{Chandra Erdman and John W. Emerson}

\note{\code{return.mcmc} must be \code{TRUE}}


\seealso{\code{\link{bcp}} and \code{\link{plot.bcp}}.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
  bcp.0 <- bcp(testdata, return.mcmc=TRUE)
  plot.bcp(bcp.0)
  interval.prob(bcp.0, 45, 55)
  
}

\keyword{datasets}
