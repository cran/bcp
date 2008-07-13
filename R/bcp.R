"worker.bcp" <- function(mcmc, x, w0, p0, burnin, return.mcmc) {

  require(bcp)
  return(bcp(x, w0, p0, burnin, mcmc, return.mcmc))

}

"bcp" <-
function(x, w0=0.2, p0=0.2, burnin=50, mcmc=500, return.mcmc=FALSE, nwssleigh = NULL) {

  if (!is.null(nwssleigh)) {

    if (!require(nws)) stop("The nws package is not available.")

    numworkers <- length(nwssleigh@nodeList)
    
    if (numworkers < 2) stop("Parallel computing requires at least 2 workers.")
    
    mcmc.w <- rep(floor(mcmc/numworkers), numworkers)
    mcmc.w[1] <- mcmc - sum(mcmc.w[2:numworkers])

    ans <- eachElem(nwssleigh, worker.bcp,
                    elementArgs=mcmc.w,
                    fixedArgs=list(x, w0, p0, burnin, return.mcmc))

    # Here, process the result and return! 

    gmeans <- ans[[1]]$posterior.mean
    gprob <- ans[[1]]$posterior.prob
    gvar <- ans[[1]]$posterior.var * (mcmc.w[1]-1)
    for (i in 2:numworkers) {
      gmeans <- gmeans + ans[[i]]$posterior.mean
      gprob <- gprob + ans[[i]]$posterior.prob
      gvar <- gvar + ans[[i]]$posterior.var * (mcmc.w[i]-1)
    }
    gmeans <- gmeans / numworkers
    gprob <- gprob / numworkers
    # gvar at this point is the within-group SS (WSS)
    bss <- mcmc.w[1] * (ans[[1]]$posterior.mean - gmeans)^2
    for (i in 2:numworkers) {
      bss <- bss + mcmc.w[i] * (ans[[i]]$posterior.mean - gmeans)^2
    }
    gvar <- gvar + bss
    gvar <- gvar / (mcmc-1)

    gbcp <- ans[[1]]
    gbcp$posterior.mean <- gmeans
    gbcp$posterior.prob <- gprob
    gbcp$posterior.var <- gvar
    gbcp$mcmc <- mcmc
    gbcp$burnin <- burnin * numworkers
    bblocks <- gbcp$blocks[1:burnin]
    allblocks <- gbcp$blocks[(burnin+1):length(gbcp$blocks)]
    for (i in 2:numworkers) {
      bblocks <- c(bblocks, ans[[i]]$blocks[1:burnin])
      allblocks <- c(allblocks, ans[[i]]$blocks[(burnin+1):length(ans[[i]]$blocks)])
    }
    gbcp$blocks <- c(bblocks, allblocks)
    if (return.mcmc) {
      lastpos <- mcmc.w[1] + burnin
      mcmc.mburn <- gbcp$mcmc.means[1:burnin,]
      mcmc.means <- gbcp$mcmc.means[(burnin+1):lastpos,]
      mcmc.rburn <- gbcp$mcmc.rhos[1:burnin,]
      mcmc.rhos <- gbcp$mcmc.rhos[(burnin+1):lastpos,]
      for (i in 2:numworkers) {
        lastpos <- mcmc.w[i] + burnin
        mcmc.mburn <- rbind(mcmc.mburn, ans[[i]]$mcmc.means[1:burnin,])
        mcmc.means <- rbind(mcmc.means,
                            ans[[i]]$mcmc.means[(burnin+1):lastpos,])
        mcmc.rburn <- rbind(mcmc.rburn, ans[[i]]$mcmc.rhos[1:burnin,])
        mcmc.rhos <- rbind(mcmc.rhos,
                            ans[[i]]$mcmc.rhos[(burnin+1):lastpos,])
      }
      gbcp$mcmc.means <- rbind(mcmc.mburn, mcmc.means)
    }

    return(gbcp)

  } else {

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
	
	# Do the work in C: 
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
        posterior.prob <- out$pchange
        posterior.prob[length(x)] <- NA
	
	# RETURN RESULTS
	y <- (list(data=x,
		return.mcmc=return.mcmc,
		mcmc.means=results,
		mcmc.rhos=rhos,
		blocks=out$blocks,
		posterior.mean=out$pmean,
		posterior.var=out$pvar,
		posterior.prob=posterior.prob,
		burnin=burnin,  
		mcmc=mcmc,     
		p0=p0,		 
                w0=w0)		 
                )
	class(y) <- "bcp"
	return(y)

  } # End the single processing.

}
