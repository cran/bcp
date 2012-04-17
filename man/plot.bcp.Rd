% --- Source file: plot.bcp.Rd ---
\name{plot.bcp}

\alias{plot.bcp}

\title{ Plotting Bayesian change point results}

\description{
  \code{plot.bcp()} produces summary plots of the results of \code{bcp()}.
}

\usage{
 \method{plot}{bcp}(x, separated = FALSE, 
   outer.margins = list(left=unit(4, "lines"),
     bottom=unit(3, "lines"),
     right=unit(2, "lines"), 
     top=unit(2, "lines")),
   lower.area = unit(0.33, "npc"),
   size.points = unit(0.25, "char"),
   pch.points = 20,
   colors = NULL,
   main = NULL,      
   cex.axes = list(cex.xaxis = 0.75,
     cex.yaxis.lower = 0.75,
     cex.yaxis.upper.default = 0.75,
     cex.yaxis.upper.separated = 0.5),
    ...)
}

\arguments{
  \item{x}{the result of a call to \code{bcp()}.}
  \item{separated}{if \code{separated=TRUE} and the data is multivariate, each series is plotted separately.}
  \item{outer.margins}{an optional list of units specifying the left, bottom, right and top margins.  For more information on units, see the documentation for \code{grid}.}
  \item{lower.area}{an optional unit specifying the proportion of the plot occupied by the posterior probabilities of change points.}
  \item{size.points}{an optional unit specifying the size of the data points.}
  \item{pch.points}{an optional unit specifying the style of the data points.}
  \item{colors}{an optional vector specifying the colors in which to plot each data series.}
  \item{main}{an optional plot title. Use \code{""} for no title.}
  \item{cex.axes}{an optional list specifying the sizes of the axes labels. \code{cex.xaxis} specifies the size of the x-axis label, \code{cex.yaxis.lower} specifies the size of the y-axis label of the posterior probability plot, \code{cex.yaxis.upper.default} specifies the size of the y-axis labels of the posterior means plot when the series are displayed in a single plot, and \code{cex.yaxis.upper.separated} specifies the size of the y-axis labels of the posterior means plots when each series is plotted separately.}
  \item{...}{additional arguments, ignored.}  
}

\details{
\code{plot.bcp()} produces the following plots using \code{grid} graphics instead of \code{base}:
     
Posterior Means: location in the sequence versus the posterior means over the iterations. 

Posterior Probability of a Change: location in the sequence versus the relative frequency of iterations which resulted in a change point. 

}

\author{Chandra Erdman and John W. Emerson}

\note{The original plot method for \code{bcp} from package versions prior to 3.0.0 is available in \code{\link{legacyplot}}; this plot method uses \code{base} graphics.}


\seealso{\code{\link{legacyplot}}, \code{\link{bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for complete results and summary statistics.}

\examples{

  testdata <- cbind( c(rnorm(50), rnorm(50, -5, 1), rnorm(50)),
                     c(rnorm(50), rnorm(50, 10.8, 1), rnorm(50, -3, 1)) )
  bcp.0 <- bcp(testdata)
  plot.bcp(bcp.0)
  plot.bcp(bcp.0, separated=TRUE)

}

\keyword{datasets}
