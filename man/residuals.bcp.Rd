\name{residuals.bcp}

\alias{residuals.bcp}

\title{Extract model residuals}

\description{
 residuals method for class \code{bcp}.
}

\usage{
 residuals.bcp(object, ...)
}

\arguments{
  \item{object}{the result of a call to \code{bcp()}.}
  \item{...}{additional arguments.}
  
}

\value{
Residuals extracted from the \code{bcp} object.
}

\author{Chandra Erdman and John W. Emerson}


\seealso{\code{\link{bcp}} and \code{\link{plot.bcp}}.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
  bcp.0 <- bcp(testdata)
  plot.bcp(bcp.0)
  residuals.bcp(bcp.0)
  
}

\keyword{datasets}
