"bcp" <-
function(x, w0=0.2, p0=0.2, burnin=50, mcmc=500, return.mcmc=FALSE) {

	if (class(try(.Random.seed, silent=TRUE))=="try-error") junk <- runif(1)

	# INITIALIZATION
	n <- length(x)		 # n = sample size.             
	M <- burnin + mcmc	            
	rho <- rep(0,n)		 # rho = vector of 0/1, specifying a partition.
        rho[n] <- 1
	blocks <- rep(0,M)	 # blocks = vector of number of blocks after each iteration.
	if (return.mcmc==TRUE){
		rhos <- matrix(0,M,n)	 # rhos = matrix of rho[m] for m in 1:M.
		results <- matrix(0,M,n) # results = matrix of posterior means.
		mcmcreturn <- 1		# for use in C.
	} else {
		rhos <- 0
		results <- 0
		mcmcreturn <- 0
	}
	pmean <- rep(0,n)
	pvar <- rep(0,n)
	pchange <- rep(0,n)
	
	# LOAD C SCRIPT 
	out <- .C("Cbcp", 
		PACKAGE="bcp", 
		data = as.double(x), 
		mcmcreturn = as.integer(mcmcreturn),
		n = as.integer(n), 
		burnin = as.integer(burnin), 
		mcmc = as.integer(mcmc),     
		rho = as.integer(rho),
		rhos = as.integer(rhos),
		blocks = as.integer(blocks),
		results = as.double(results),
        	a = as.double(p0),
		c = as.double(w0),
		pmean = as.double(pmean),
		pvar = as.double(pvar),
		pchange = as.double(pchange)
      	)

	if (return.mcmc==TRUE) {
		# STUFF LONG VECTORS FROM C INTO MATRICES IN R
		start <- rep(0,M)
		end <- rep(0,M)
		for (m in 1:M) {
			start[m] <- (m-1)*n + 1
			end[m] <- m*n
			rhos[m,] <- out$rhos[start[m]:end[m]]
			results[m,] <- out$results[start[m]:end[m]]
		}
	} else {
		rhos <- NA
		results <- NA
	}
	
	# RETURN RESULTS
	y <- (list(data=x,
		return.mcmc=return.mcmc,
		mcmc.means=results,
		mcmc.rhos=rhos,
		blocks=out$blocks,
		posterior.mean=out$pmean,
		posterior.var=out$pvar,
		posterior.prob=out$pchange,
		burnin=burnin,  
		mcmc=mcmc,     
		p0=p0,		 
                w0=w0)		 
                )
	class(y) <- "bcp"
	return(y)
}
