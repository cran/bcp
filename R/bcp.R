#
# bcp: an R package for performing a Bayesian analysis
# of change point problems.
#
# Copyright (C) 2011 Chandra Erdman and John W. Emerson
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, a copy is available at
# http://www.r-project.org/Licenses/
#
#-------------------
# FILE: bcp.R

"bcp" <- function(x, w0=0.2, p0=0.2, burnin=50, mcmc=500, return.mcmc=FALSE) {

######################################################
########################### BEGIN THE WORKER FUNCTION:
######################################################

"worker.bcp" <- function(mcmc, x, w0, p0, burnin, return.mcmc) {

  require(bcp)

  # INITIALIZATION
  if (is.data.frame(x)) x <- matrix(as.double(x), nrow=nrow(x), ncol=ncol(x))
  if (is.vector(x)) x <- matrix(as.double(x), ncol=1)
  if (!is.matrix(x)) stop("x must be a vector, matrix or a data frame")
  if (nrow(x)==1) {
    warning("coercing data to a single series")
    x <- matrix(as.vector(x), ncol=1)
  }

  # Do the work in C:
  out <- .Call("rcpp_bcp", 
		PACKAGE="bcp", 
		data = x,
		mcmcreturn = as.integer(return.mcmc),
		burnin = as.integer(burnin), 
		mcmc = as.integer(mcmc),     
		a = as.double(p0),
		c = as.double(w0)
      	      )

  out$pchange[nrow(x)] <- NA     # Fix up the last position, always NA
	
  # RETURN RESULTS
  y <- list(data=x,
		return.mcmc=return.mcmc,
		mcmc.means=out$mcmc.means,
		mcmc.rhos=out$mcmc.rhos,
		blocks=out$blocks,
		posterior.mean=out$pmean,
		posterior.var=out$pvar,
		posterior.prob=out$pchange,
		burnin=burnin,            
		mcmc=mcmc,     
		p0=p0,		 
		w0=w0)		 
  class(y) <- "bcp"

  return(y)

}

###################################################
########################### END THE WORKER FUNCTION
########################### BEGIN THE MAIN SECTION:
###################################################

# Function header and foreach setup, from above:
#
#"bcp" <- function(x, w0=0.2, p0=0.2, burnin=50, mcmc=500, return.mcmc=FALSE) {
#

  require(foreach)
  if (is.null(getDoParName())) {
    registerDoSEQ() # A little hack to avoid the foreach warning 1st time.
  }

  # Divy up the work across available workers:
  numworkers <- getDoParWorkers()
  mcmc.w <- rep(floor(mcmc/numworkers), numworkers)
  if (numworkers>1) { # Give the first worker a little extra if needed.
    mcmc.w[1] <- mcmc - sum(mcmc.w[2:numworkers])
  }

  ######################################################
  # The call to do the work, perhaps in parallel:

  ans <- foreach(mcmc=mcmc.w) %dopar% {
    worker.bcp(mcmc, x=x, w0=w0, p0=p0, burnin=burnin, return.mcmc=return.mcmc)
  }

  ######################################
  # Process the result and return. 

  gmeans <- ans[[1]]$posterior.mean
  gprob <- ans[[1]]$posterior.prob
  gvar <- ans[[1]]$posterior.var * (mcmc.w[1]-1)
  if (numworkers>1) {
    for (i in 2:numworkers) {
      gmeans <- gmeans + ans[[i]]$posterior.mean
      gprob <- gprob + ans[[i]]$posterior.prob
      gvar <- gvar + ans[[i]]$posterior.var * (mcmc.w[i]-1)
    }
    gmeans <- gmeans / numworkers
    gprob <- gprob / numworkers
  }

  # gvar at this point is the within-group SS (WSS)
  bss <- mcmc.w[1] * (ans[[1]]$posterior.mean - gmeans)^2
  if (numworkers>1) {
    for (i in 2:numworkers) {
      bss <- bss + mcmc.w[i] * (ans[[i]]$posterior.mean - gmeans)^2
    }
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
  if (numworkers>1) {
    for (i in 2:numworkers) {
      bblocks <- c(bblocks, ans[[i]]$blocks[1:burnin])
      allblocks <- c(allblocks, ans[[i]]$blocks[(burnin+1):length(ans[[i]]$blocks)])
    }
  }
  gbcp$blocks <- c(bblocks, allblocks)

  if (return.mcmc) {                ### Non-trivial because of multivariate as well
                                    ### as possible multiple workers...
    lastpos <- mcmc.w[1] + burnin

    # Handle the rhos, burnin and real mcmc partition states:
    mcmc.rburn <- gbcp$mcmc.rhos[,1:burnin]
    mcmc.rhos <- gbcp$mcmc.rhos[,(burnin+1):lastpos]

    nsamples <- ncol(gbcp$data)
    n <- nrow(gbcp$data)

    # Ditto for the results, but possibly multivariate:
    temp <- matrix(gbcp$mcmc.means[,1], nrow=n, ncol=lastpos)
    mcmc.mburn <- temp[,1:burnin]
    mcmc.means <- temp[,(burnin+1):lastpos]
    if (nsamples>1) {
      mcmc.mburn <- list(mcmc.mburn)
      mcmc.means <- list(mcmc.means)
      for (i in 2:nsamples) {
        temp <- matrix(gbcp$mcmc.means[,i], nrow=n, ncol=lastpos)
        mcmc.mburn[[i]] <- temp[,1:burnin]
        mcmc.means[[i]] <- temp[,(burnin+1):lastpos]
      }
    }

    if (numworkers>1) {
      for (j in 2:numworkers) {
        lastpos <- mcmc.w[j] + burnin
        temp <- matrix(ans[[j]]$mcmc.means[,1], nrow=n, ncol=lastpos)
        # if mcmc.mburn is a matrix, simply cbind stuff from the workers
        if (nsamples==1) {
          mcmc.mburn <- cbind(mcmc.mburn, temp[,1:burnin])
          mcmc.means <- cbind(mcmc.means, temp[,(burnin+1):lastpos])
        } else {
          # if mcmc.mburn is a list of matrices, do the work across samples,
          # being very careful of the list structure: a list with matrices for
          # each data series, whereas we're unpacking the worker results from
          # a worker list, as well... be very careful!
          mcmc.mburn[[1]] <- cbind(mcmc.mburn[[1]], temp[,1:burnin])
          mcmc.means[[1]] <- cbind(mcmc.means[[1]], temp[,(burnin+1):lastpos])
          for (i in 2:nsamples) {
            temp <- matrix(ans[[j]]$mcmc.means[,i], nrow=n, ncol=lastpos)
            mcmc.mburn[[i]] <- cbind(mcmc.mburn[[i]], temp[,1:burnin])
            mcmc.means[[i]] <- cbind(mcmc.means[[i]], temp[,(burnin+1):lastpos])
          }
        }
      }
    }

    if (nsamples==1) {
      gbcp$mcmc.means <- cbind(mcmc.mburn, mcmc.means)
    } else {
      gbcp$mcmc.means <- list(burnins=mcmc.mburn, means=mcmc.means)
    }
  }

  return(gbcp)

}

