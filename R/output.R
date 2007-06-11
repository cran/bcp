# the code for spectrum0, do.spectrum0, safespec0, and summary.bcp was taken from the coda package 
# (Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines, 2006).

####################################################
#spectrum0 <- function(x, max.freq=0.5, order=1, max.length=200)
#{
#  x <- as.matrix(x$results)
#  if (!is.null(max.length) && nrow(x) > max.length) {
#    batch.size <- ceiling(nrow(x)/max.length)
#    x <- aggregate(ts(x, frequency=batch.size), nfreq = 1, FUN=mean)
#  }
#  else {
#    batch.size <- 1
#  }
#  out <- do.spectrum0(x, max.freq=max.freq, order=order)
#  out$spec <- out$spec * batch.size
#  return(out)
#}

####################################################
#do.spectrum0 <- function(x, max.freq=0.5, order=1)
#{
#  ## Estimate spectral density of time series x at frequency 0.
#  ## spectrum0(x)/length(x) estimates the variance of mean(x)
#  
#  fmla <- switch(order+1,
#                 spec ~ one,
#                 spec ~ f1,
#		   spec ~ f1 + f2)
#  if(is.null(fmla))
#    stop("invalid order")
#
#  N <- nrow(x)
#  Nfreq <- floor(N/2)
#  freq <- seq(from = 1/N, by = 1/N, length = Nfreq)
#  f1 <- sqrt(3) * (4 * freq - 1)
#  f2 <- sqrt(5) * (24 * freq^2 - 12 * freq + 1)
#  v0 <- numeric(ncol(x))
#  for(i in 1:ncol(x)) {
#    y <- x[,i]
#    if (var(y) == 0) {
#      v0[i] <- 0
#    }
#    else {
#      yfft <- fft(y)
#      spec <- Re(yfft * Conj(yfft))/ N
#      spec.data <- data.frame(one = rep(1, Nfreq), f1=f1, f2=f2,
#                              spec = spec[1 + (1:Nfreq)],
#                              inset = I(freq<=max.freq))
#      
#      glm.out <- glm(fmla, family=Gamma(link="log"), data=spec.data,
#                     subset=inset)
#      v0[i] <- predict(glm.out, type="response",
#                       newdata=data.frame(spec=0,one=1,f1=-sqrt(3),f2=sqrt(5)))
#    }
#  }
#  return(list(spec=v0))
#}

###############################
#safespec0 <-
#  function (x) {
#  result <- try(spectrum0(x)$spec)
#  ## R
#  if (class(result) == "try-error") result <- NA
#  ## S-Plus
#  if (class(result) == "try") result <- NA
#  result
#}

##########################################################
#summary.bcp <-
#  function (object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = max(3, .Options$digits - 3), ...) 
#{
#	results <- as.matrix(object$results)
#	rhos <- as.matrix(object$rhos)
#  
#	# calculations from summary.mcmc
#	statnames <- c("Probability", "Mean", "SD", "Naive SE", "Time-series SE")
#	varstats <- matrix(nrow = ncol(results), ncol = length(statnames), dimnames = list(1:ncol(results), statnames))
#	xprob <- apply(rhos[object$burnin:nrow(results),], 2, mean)
#	xmean <- apply(results[object$burnin:nrow(results),], 2, mean)
#	xvar <- apply(results[object$burnin:nrow(results),], 2, var)
#	xtsvar <- safespec0(object)
#	varquant <- t(apply(results[object$burnin:nrow(results),], 2, quantile, quantiles))
#	varstats[, 1] <- xprob
#	varstats[, 2] <- xmean
#	varstats[, 3] <- sqrt(xvar)
#	varstats[, 4] <- sqrt(xvar/object$mcmc)
#	varstats[, 5] <- sqrt(xtsvar/object$mcmc)
#	varstats <- drop(varstats)
#	varquant <- drop(varquant)
#  
#	# print everything
#	out <- list(statistics = varstats, quantiles = varquant)
#	cat("\nBayesian Change Point (bcp) output:\n\n")
#	cat("\n1. Posterior probability of a change in mean and Posterior mean for each position,")
#	cat("\n   plus standard deviation and standard error of the mean:\n\n")
#      print(out$statistics, digits=digits)
#	cat("\n2. Quantiles for each position:\n\n")
#	print(out$quantiles, digits = digits)
#	cat("\n")
#}

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


