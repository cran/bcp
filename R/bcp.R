"bcp" <-
function(x, w0=0.2, p0=0.2, burnin=50, mcmc=500) {

	# INITIALIZATION
	n <- length(x)		# n = sample size.             
	M <- burnin + mcmc	# REMOVE THIS LINE.                 
	rho <- rep(0,n)		# rho = vector of 0/1, specifying a partition.
      rho[n] <- 1
	rhos <- matrix(0,M,n)	# rhos = matrix of rho[m] for m in 1:M.
	blocks <- rep(0,M)	# blocks = vector of number of blocks after each iteration.
	results <- matrix(0,M,n)	# results = matrix of posterior means.

	# LOAD C SCRIPT 
	out <- .C("Cbcp", 
		PACKAGE="bcp", 
		data = as.double(x), 
		n = as.integer(n), 
		# M = as.integer(M), 	     # REMOVE THIS LINE 
		burnin = as.integer(burnin), # ADDED THIS LINE
		mcmc = as.integer(mcmc),     # ADDED THIS LINE
		rho = as.integer(rho),
		rhos = as.integer(rhos),
		blocks = as.integer(blocks),
		results = as.double(results),
        	a = as.double(p0),
		c = as.double(w0)
      	)

	# STUFF LONG VECTORS FROM C INTO MATRICES IN R
	start <- rep(0,M)
	end <- rep(0,M)
	for (m in 1:M) {
  		start[m] <- (m-1)*n + 1
  		end[m] <- m*n
  		rhos[m,] <- out$rhos[start[m]:end[m]]
  		results[m,] <- out$results[start[m]:end[m]]
	}

	# RETURN RESULTS
	return(list(data=x,
		   results=results,
               rhos=rhos,
               blocks=out$blocks,
               posterior.mean=apply(results[burnin:M,1:n],2,mean),
               burnin=burnin,  # ADDED THIS LINE
		   mcmc=mcmc,      # ADDED THIS LINE
	         p0=p0,		 # ADDED THIS LINE
               w0=w0)		 # ADDED THIS LINE
              )
}
