% --- Source file: summary.bcp.Rd ---
\name{summary.bcp}
\alias{print.bcp}
\alias{summary.bcp}

\title{Summarizing Bayesian change point analysis results}

\description{
 Summary and print methods for class \code{bcp}.
}

\usage{
 \method{summary}{bcp}(object, digits = max(3, .Options$digits - 3), ...)
 \method{print}{bcp}(x, digits = max(3, .Options$digits - 3), ...)
}

\arguments{
  \item{x}{the result of a call to \code{bcp()}.}
  \item{object}{the result of a call to \code{bcp()}.}
  \item{digits}{the number of digits displayed in the summary statistics.}
  \item{...}{additional arguments.}
}

\value{The matrix of results is returned invisibly.}

\details{
The functions print (and return invisibly) the estimated posterior probability of a change point for each position and the estimated posterior means. These results are modeled after the summary method of the \code{coda} package (Plummer \emph{et al.}, 2006). If \code{return.mcmc=TRUE} (i.e., if full MCMC results are returned), \code{bcp} objects can be converted into \code{mcmc} objects to view \code{mcmc} summaries -- see examples below.
}

\author{Chandra Erdman and John W. Emerson}


\seealso{\code{\link{bcp}} and \code{\link{plot.bcp}}.}

\examples{

  ##### A random sample from a few normal distributions #####
  testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
  bcp.0 <- bcp(testdata)
  summary.bcp(bcp.0)
  plot.bcp(bcp.0)
  
  ##### An MCMC summary from the ``coda'' package #####
  if (require("coda")) {
    bcp.0 <- bcp(testdata, return.mcmc=TRUE)
    bcp.mcmc <- as.mcmc(t(bcp.0$mcmc.means))
    summary(bcp.mcmc)
    heidel.diag(bcp.mcmc) # an example convergence diagnostic
                          # from the coda package.
  }
}

\references{
Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines (2006), \emph{The coda Package, version 0.10-7}, CRAN: The Comprehensive R Network. 
}

\keyword{datasets}
