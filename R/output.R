summary.bcp <-
	function (object, digits = max(3, .Options$digits - 3), ...) {
	
	statnames <- c("Probability", "Mean", "SD")
	varstats <- matrix(nrow = length(object$data), ncol = length(statnames), dimnames = list(1:length(object$data), statnames))
	xprob <- object$posterior.prob
	xmean <- object$posterior.mean
	xvar <- object$posterior.var
	varstats[, 1] <- xprob
	varstats[, 2] <- xmean
	varstats[, 3] <- sqrt(xvar)
	varstats <- drop(varstats)

	# print everything
	out <- list(statistics = varstats)
	cat("\nBayesian Change Point (bcp) output:\n\n")
	cat("\nProbability of a change in mean, Posterior Mean,\n")
	cat("and SD for each position:\n\n")
        print(out$statistics, digits=digits)
	cat("\n")
}
	

##############################################################################
print.bcp <- function(x, digits = max(3, .Options$digits - 3), ...) {
	
	statnames <- c("Probability", "Mean", "SD")
	varstats <- matrix(nrow = length(x$data), ncol = length(statnames), dimnames = list(1:length(x$data), statnames))
	xprob <- x$posterior.prob
	xmean <- x$posterior.mean
	xvar <- x$posterior.var
	varstats[, 1] <- xprob
	varstats[, 2] <- xmean
	varstats[, 3] <- sqrt(xvar)
	varstats <- drop(varstats)
  
	# print everything
	out <- list(statistics = varstats)
	cat("\nProbability of a change in mean, Posterior Mean,\n")
	cat("and SD for each position:\n\n")
	print(out$statistics, digits=digits)
	cat("\n")
}

########################################################################
plot.bcp <-
	function(x, ...) { 	
	posterior.prob <- x$posterior.prob
	posterior.prob[length(posterior.prob)] <- 0
		
	op<-par(mfrow=c(2,1),col.lab="black",col.main="black")
		op2 <- par(mar=c(0,4,4,2),xaxt="n", cex.axis=0.75)
			plot(1:length(x$data), x$data, col="grey", pch=20, xlab="", ylab="Posterior Mean", main="Posterior Means and Probabilities of a Change", ...)
			lines(x$posterior.mean, lwd=2)
		par(op2)
		op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
			plot(1:length(x$posterior.mean), posterior.prob, yaxt="n", type="l", ylim=c(0,1),
			xlab="Location", ylab="Posterior Probability", main="")
			axis(2, yaxp=c(0, 0.9, 3))
		par(op3)
	par(op)
	
}

########################################################################
residuals.bcp <-
	function(object, ...) { 
        residuals <- object$posterior.mean-object$data	
	print(residuals)	
}

########################################################################
fitted.bcp <-
	function(object, ...) {
        fitted <- object$posterior.mean	
	print(fitted)	
}


