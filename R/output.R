summary.bcp <-
  function (object, digits = max(3, .Options$digits - 3), ...) {
  
  if (!is.matrix(object$data)) {
    object$data <- matrix(object$data, ncol=1)
    names(object$data) <- "Mean" 
  }
  
  if (is.null(colnames(object$data))) {
    colnames(object$data) <- paste(rep("X", ncol(object$data)), 
                                   1:ncol(object$data), sep="")
  }

  statnames <- c("Probability", colnames(object$data))
  varstats <- matrix(NA, nrow = nrow(object$data), 
                     ncol = length(statnames), 
                     dimnames = list(1:nrow(object$data), statnames))
						   
  varstats[,1] <- object$posterior.prob
  varstats[,2:ncol(varstats)] <- object$posterior.mean

  cat("\nBayesian Change Point (bcp) summary:\n\n")
  cat("\nProbability of a change in mean and Posterior Means:\n\n")
  print(varstats, digits=digits)
  cat("\n")

  invisible(varstats)
}

print.bcp <- function(x, digits = max(3, .Options$digits - 3), ...) {
  summary(x, digits=digits, ...)
}
	
interval.prob <- function(object, start, end) {
  if (!object$return.mcmc) stop("bcp must be run with return.mcmc=TRUE")
  return( sum(apply(object$mcmc.rhos[start:end,-c(1:object$burnin)], 2, sum) > 0) /
          object$mcmc)
}

plot.bcp.legacy <- function(x, ...) { 	
  if (is.matrix(x$data) && ncol(x$data)>1)
    stop("Legacy bcp plot invalid for multivariate bcp object.")
  posterior.prob <- x$posterior.prob
  posterior.prob[length(posterior.prob)] <- 0
		
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  op2 <- par(mar=c(0,4,4,2),xaxt="n", cex.axis=0.75)
  plot(1:length(x$data), x$data, col="grey", 
       pch=20, xlab="", ylab="Posterior Mean", 
       main="Posterior Means and Probabilities of a Change", ...)
  lines(x$posterior.mean, lwd=2)
  par(op2)
  op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
  plot(1:length(x$posterior.mean), posterior.prob, 
       yaxt="n", type="l", ylim=c(0,1),
       xlab="Location", ylab="Posterior Probability", main="")
  axis(2, yaxp=c(0, 0.9, 3))
  par(op3)
  par(op)
}

plot.bcp <- function(x, separated = FALSE, 
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
                     ...) {
  if (!require(grid)) stop("library(grid) is required and unavailable.\n\n")
  grid.newpage()

  thisjust <- c("left", "bottom")
  n <- length(x$posterior.prob)
  if (!is.matrix(x$data)) {
    x$data <- matrix(x$data, ncol=1) 
    x$posterior.mean <- matrix(x$posterior.mean, ncol=1)
  }
  m <- ncol(x$data)
  if (is.null(colors)) {
    colors <- 2:(m+1)
  } else {
    if (length(colors)!=m) stop("length(colors) must equal number of series.")
  }
  if (is.null(main)) main <- "Posterior Means and Probabilities of a Change"

  vp.main <- viewport(x=outer.margins$left, y=outer.margins$bottom,
                      w=unit(1, "npc")-outer.margins$right-outer.margins$left,
                      h=unit(1, "npc")-outer.margins$top-outer.margins$bottom,
                      just=thisjust, name="main", clip="off")
  pushViewport(vp.main)
  grid.rect()
  grid.text(main, 0.5, unit(1, "npc")+unit(1,"lines"), 
            gp=gpar(fontface="bold"))

  # Upper plot
  if (!separated) {
    pushViewport(viewport(x=unit(0, "npc"), y=lower.area,
                          w=unit(1, "npc"), h=unit(1, "npc") - lower.area,
                          just=thisjust, name="upper", clip="off",
                          default.units="native",
                          xscale=c(-1, n+2),
                          yscale=range(pretty(x$data, n=10))))
    grid.rect()
    grid.yaxis(gp=gpar(cex=cex.axes$cex.yaxis.upper.default), main=FALSE)
    grid.text("Posterior Means", unit(-3, "lines"), 0.5, rot=90)
    for (i in 1:m) {
      grid.points(1:n, x$data[,i], size=size.points, gp=gpar(col=colors[i]),
                  pch=pch.points)
      grid.lines(1:n, x$posterior.mean[,i], gp=gpar(col=colors[i]),
                 default.units="native")
    }
    popViewport(1)
  } else {
    pushViewport(viewport(x=unit(0, "npc"), y=lower.area,
                          w=unit(1, "npc"), h=unit(1, "npc") - lower.area,
                          just=thisjust, name="upper", clip="off"))
    grid.text("Posterior Means", unit(-3, "lines"), 0.5, rot=90)
    yloc <- FALSE
    for (i in 1:m) {
      pushViewport(viewport(x=unit(0, "npc"),
                            y=unit((i-1)/m, "npc"),
                            w=unit(1, "npc"),
                            h=unit(1/m, "npc"),
                            just=thisjust, name="upper", clip="off",
                            default.units="native",
                            xscale=c(-1, n+2),
                            yscale=range(pretty(x$data[,i]))))
      grid.rect()
      grid.yaxis(gp=gpar(cex=cex.axes$cex.yaxis.upper.separated), main=yloc)
      grid.points(1:n, x$data[,i], size=size.points, gp=gpar(col=colors[i]),
                  pch=pch.points)
      grid.lines(1:n, x$posterior.mean[,i], gp=gpar(col=colors[i]),
                    default.units="native")
      popViewport(1)
      yloc <- !yloc
    }
    popViewport(1)
  }

  # Lower plot
  pushViewport(viewport(x=unit(0, "npc"), y=unit(0, "npc"),
                        w=unit(1, "npc"), h=lower.area,
                        just=thisjust, name="lower", clip="off",
                        default.units="native",
                        xscale=c(-1, n+2),
                        yscale=c(-0.05, 1.05)))
  grid.rect()
  grid.yaxis(gp=gpar(cex=cex.axes$cex.yaxis.lower))
  grid.xaxis(gp=gpar(cex=cex.axes$cex.xaxis))
  grid.text("Posterior Probability", unit(-3, "lines"), 0.5, rot=90)
  grid.text("Location", 0.5, unit(-2, "lines"))
  grid.lines(1:n, x$posterior.prob, default.units="native")
  popViewport(2) # vp.main and vp.lower

}

residuals.bcp <- function(object, ...) { 
  return(object$data - object$posterior.mean)
}

fitted.bcp <- function(object, ...) {
  return(object$posterior.mean)
}


