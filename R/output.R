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
	plot(1:length(x$data), x$data, xlab="Location", ylab="Posterior Mean", main="Posterior Means")
		lines(x$posterior.mean)
	plot(1:length(x$posterior.mean), posterior.prob, type="l", ylim=c(0,1),
		xlab="Location", ylab="Posterior Probability", main="Posterior Probability of a Change")
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


