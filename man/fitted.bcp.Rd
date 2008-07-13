\name{fitted.bcp}

\alias{fitted.bcp}

\title{Extract model fitted values}

\description{
 fitted method for class \code{bcp}.
}

\usage{
 fitted.bcp(object, ...)
}

\arguments{
  \item{object}{the result of a call to \code{bcp()}.}
  \item{...}{additional arguments.}
  
}

\value{
Fitted values extracted from the \code{bcp} object.
}

\author{Chandra Erdman and John W. Emerson}


\seealso{\code{\link{plot.bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for summaries of the results.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
  bcp.0 <- bcp(testdata)
  plot.bcp(bcp.0)
  fitted.bcp(bcp.0)
  
}

\keyword{datasets}
