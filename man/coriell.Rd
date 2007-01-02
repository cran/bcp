% --- Source file: coriell.Rd ---
\name{coriell}

\alias{coriell}

\docType{data}

\title{ Array CGH data set of Coriell cell lines }

\description{
 These are two data array CGH studies of Corriel cell lines taken from the reference below.
}
\usage{data(coriell)}

\format{ A data frame containing five variables: first is clone name, second is clone chromosome, third is   clone position, fourth and fifth  are log2ratio for two cell lines.
}
  
\source{ \url{htp://www.nature.com/ng/journal/v29/n3/suppinfo/ng754\_S1.html} }

\references{Snijders \emph{et al.} (2001), Assembly of microarrays for genome-wide measurement of DNA copy number, \emph{Nature Genetics}, \bold{29}, 263-264. 

}

\examples{
  ##### Coriell chromosome 11 #####
  data(coriell)
  chrom11 <- na.omit(coriell$Coriell.05296[coriell$Chromosome==11])
  n <- length(chrom11)
  bcp.11 <- bcp(chrom11[1:n])
  summary.bcp(bcp.11)
  plot.bcp(bcp.11)
  
  # to see bcp and Circular Binary Segmentation results run:
  if(require("DNAcopy")) {
  bcp.11$rhos[,ncol(bcp.11$rhos)] <- 0
  changepoint.freq <- apply(bcp.11$rhos[bcp.11$burnin:(bcp.11$burnin+bcp.11$burnin),1:dim(bcp.11$results)[2]],2,mean)
  cbs <- segment(CNA(chrom11[1:n], rep(1, n), 1:n), verbose = 0)
  cbs.ests <- rep(unlist(cbs$output[6]), unlist(cbs$output[5]))
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  plot(1:n, bcp.11$posterior.mean, type="l", xlab="Location", ylab="Posterior Mean", main="Posterior Means")
  lines(cbs.ests, col="red")
  points(chrom11)
  plot(1:n, changepoint.freq, type="l", xlab="Location", ylab="Posterior Probability of a Change", main="Change Point Locations")
  for(i in 1:(dim(cbs$output)[1]-1)) abline(v=cbs$output$loc.end[i], col="red")
  par(op)
  } else {
	cat("DNAcopy is not loaded")
  }

}
  
\keyword{datasets} 

