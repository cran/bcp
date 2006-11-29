% --- Source file: bcp.Rd ---
\name{bcp}

\alias{bcp}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Bayesian Change Point }

\description{
 'bcp' implements the Barry and Hartigan (1993) product partition model for the standard change point problem using Markov Chain Monte Carlo.  
}

\usage{
 bcp(x, w0 = 0.2, p0 = 0.2, burnin = 50, mcmc = 500)
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
}

\details{
 This algorithm is used when there exists an unknown partition of a sequence into contiguous blocks such that the mean is constant within blocks.}

\value{
  'bcp(x)' returns a list containing the following components:
  \item{data}{a vector copy of the data.}
  \item{rhos}{a matrix containing the partitions after each pass through the data.}
  \item{blocks}{a vector of the number of blocks after each pass.}
  \item{posterior.means}{a vector of the posterior means over the M passes.}
  The generic accessor functions 'rhos', 'posterior.means', and 'blocks' extract various useful components of the list returned by 'bcp'.  The function 'summary.bcp' is used to obtain and print plots of the results.
}

\author{Chandra Erdman}

\note{The vector of data must contain no missing or "NA" values.}

\seealso{'plots.bcp' for plots of the results. }

\examples{

# A random sample from a few normal distributions.
  testdata <- c(rnorm(10), rnorm(5, 5, 1), rnorm(10))
  bcp.0 <- bcp(testdata)
  plots.bcp(bcp.0)

# A couple real-life examples.
#  data(coriell)
#  chrom10 <- na.omit(coriell$Coriell.05296[coriell$Chromosome==10])
#  chrom11 <- na.omit(coriell$Coriell.05296[coriell$Chromosome==11])
#  bcp.10 <- bcp(chrom10[1:length(chrom10)])
#  bcp.11 <- bcp(chrom11[1:length(chrom11)])
#  par(mfrow=c(2,1))
#  plot(chrom10,ylab="Log2 of Copy Number", main="Coriell Chromosome 10")
#  lines(bcp.10$posterior.means, col=2)
#  plot(chrom11,ylab="Log2 of Copy Number", main="Coriell Chromosome 11")
#  lines(bcp.11$posterior.means, col=2)
   
}

\references{
Daniel Barry and J. A. Hartigan (1993) A Bayesian Analysis for Change Point Problems. \emph{Journal of The American Statistical Association}, 88, 309-19.

Snijders et al. (2001) Assembly of microarrays for genome-wide measurement of DNA copy number, \emph{Nature Genetics}.
}

\keyword{datasets}% at least one, from doc/KEYWORDS
