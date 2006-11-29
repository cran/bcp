% --- Source file: plots.bcp.Rd ---
\name{plots.bcp}

\alias{plots.bcp}

\title{ Bayesian Change Point Summary }

\description{
  `plots.bcp' produces summary plots of the results of the `bcp' function.
}

\usage{
 plots.bcp(x, burnin=50, mcmc=500, cbs=FALSE)
}

\arguments{
  \item{x}{the result of a call to `bcp'.}
  \item{burnin}{the number of burnin iterations.}
  \item{mcmc}{the number of iterations after burnin.}
  \item{cbs}{if TRUE, the results of the Circular Binary Segemntation algorithm will be provided in red.}  
}

\details{
`plots.bcp' produces the following plots:
     
Posterior Means: location in the sequence versus the posterior mean over the iterations. 

Change Point Locations: location in the sequence versus the relative frequency of iterations which resulted in a change point. 

Distribution of Number of Blocks: number of blocks in the partition versus the number of iterations resulting in that number of blocks.

}

\author{Chandra Erdman}


\seealso{`bcp' for complete results.}

\examples{

  # A random sample from a few normal distributions.
  testdata <- c(rnorm(10), rnorm(20, 5, 1), rnorm(10))
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
  #  plots.bcp(bcp.10)
  #  plots.bcp(bcp.11, cbs=TRUE)
  
}

\references{
Daniel Barry and J. A. Hartigan (1993) A Bayesian Analysis for Change Point Problems. \emph{Journal of The American Statistical Association}, 88, 309-19.

Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004) Circular binary segmentation for the analysis of array-based DNA copy number data. \emph{Biostatistics} *5*, 557-572. url: http://www.bioconductor.org/repository/release1.5/package/html/DNAcopy.html

Snijders et al. (2001) Assembly of microarrays for genome-wide measurement of DNA copy number, \emph{Nature Genetics}. 
}

\keyword{datasets}
