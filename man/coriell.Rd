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

\references{Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004), Circular binary segmentation for the analysis of array-based DNA copy number data, \emph{Biostatistics}, \bold{5}, 557-572. url: \url{http://www.bioconductor.org/repository/release1.5/package/html/DNAcopy.html}.

Snijders \emph{et al.} (2001), Assembly of microarrays for genome-wide measurement of DNA copy number, \emph{Nature Genetics}, \bold{29}, 263-264. 

}

\examples{
  ##### Coriell chromosome 11 #####
  data(coriell)
  chrom11 <- na.omit(coriell$Coriell.05296[coriell$Chromosome==11])
  n <- length(chrom11)
  bcp.11 <- bcp(chrom11[1:n])
  plot.bcp(bcp.11)
  
  ##### Coriell chromosome 11 #####
  data(coriell)
  chrom11 <- as.vector(na.omit(coriell$Coriell.05296[coriell$Chromosome==11]))
  bcp.11 <- bcp(chrom11)
  plot.bcp(bcp.11)
  
  # to see bcp and Circular Binary Segmentation results run:
  if(require("DNAcopy")) {
  n <- length(chrom11)
  cbs <- segment(CNA(chrom11, rep(1, n), 1:n), verbose = 0)
  cbs.ests <- rep(unlist(cbs$output[6]), unlist(cbs$output[5]))  
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  op2 <- par(mar=c(0,4,4,2),xaxt="n", cex.axis=0.75)
  plot(1:n, bcp.11$data, col="grey", pch=20, xlab="Location",
       ylab="Posterior Mean",
       main="Posterior Means and Probabilities of a Change")
  lines(cbs.ests, col="red")
  lines(bcp.11$posterior.mean, lwd=2)
  par(op2)
  op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
  plot(1:n, bcp.11$posterior.prob, type="l", ylim=c(0,1),
       xlab="Location", ylab="Posterior Probability", main="")
  for(i in 1:(dim(cbs$output)[1]-1)) {
    abline(v=cbs$output$loc.end[i], col="red")
  }
  par(op3)
  par(op)
  } else {
        cat("DNAcopy is not loaded")
  }
    
}
  
\keyword{datasets} 

