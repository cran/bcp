"bcp" <-
function(x, w0=0.2, p0=0.2, burnin=10, mcmc=100) {

  # DEFINITIONS AND INITIALIZATION
  n = length(x)             # n = sample size.
  M = burnin + mcmc         # M = number of iterations or "passes" through
                            #     the data.
  mu0 = mean(x)             # mu0 = mean of data.
  # DEFINITION              # p = probability of a change point.
  rho = rep(0,n)            # rho = vector of 0/1 that specifies a partition.
  rho[n] = 1
  i = 1                     #       if rho[i] = 1 then a block ends in this
                            #       position, starting a new one in i+1.
  b.num = 1 + c(0,cumsum(rho[1:(n-1)])) # b.num = vector of length n of block numbers. 
  b.size = table(b.num)     # b.size = vector of block sizes. 
  b = max(b.num)            # b = number of blocks we obtain if rho[i] = 0.
  # DEFINITION              # b.mean = vector of block means.
  # DEFINITION              # mu.hat = vector of length b of posterior means.
  theta = rep(0,n)          # theta = vector of length n of muhats.  
  # DEFINITION sqd = rep(0,n)           # sqd = vector of squared deviations from muhat. 
  W=0; B=0; W0=0; B0=0; W1=0; B1=0      # W,B = within and between block sums of squares
                                        #       given a partition, rho.

  integrand1 = function(p,n,b)       {p^b*(1-p)^(n-b-1)}
  integrand2 = function(w,n,b,W1,B1) {w^(b/2)/((W1 + B1*w)^((n-1)/2))}
  integrand3 = function(p,n,b)       {p^(b-1)*(1-p)^(n-b)}
  integrand4 = function(w,n,b,W0,B0) {w^((b-1)/2)/((W0 + B0*w)^((n-1)/2))}
  integrand5 = function(w,n,b,W,B)   {w^((b+1)/2)/(W + B*w)^((n-1)/2)}
  integrand6 = function(w,n,b,W,B)   {w^((b-1)/2)/(W + B*w)^((n-1)/2)}
  integrand7 = function(p,n,b)       {p^((b+2)/2-1)*(1-p)^((n-b-3)/2-1)}
  integrand8 = function(p,n,b)       {p^((b+1)/2-1)*(1-p)^((n-b-2)/2-1)}
  integrand9 = function(p,n,b,W,B)   {p^((b+3)/2-1)*(1-p)^((n-b-4)/2-1)}
  integrand10 = function(p,n,b,W,B)  {p^((b+1)/2-1)*(1-p)^((n-b-2)/2-1)} 

  # Results
  results = array(0,c(M, dim(t(rho))[2]))
  blocks = c(0,M)
  rhos = array(0,c(M, dim(t(rho))[2]))
  dimnames(results) = list(c(1:M), 1:n)

  #############################################################################
  # Start the big loop
  for (m in 1:M) {
    for (i in 1:(n-1)) {

      ###################
      # CONSIDER rho[i]=0
      rho[i] = 0
      b.num = 1 + c(0,cumsum(rho[1:(n-1)]))
      b = max(b.num)
      bsqd = rep(0,b)

      # Introduce an artificial change point far from 'i' if we have only one block.
      if (b==1) {
        if (i<n/2) rho[n-2] <- 1
        else rho[2] <- 1
        b.num = 1 + c(0,cumsum(rho[1:(n-1)]))
        b = max(b.num)
        bsqd = rep(0,b)
        temp.cp <- TRUE
      }


      # CALCULATE MU-HATS AND SUMS OF SQUARES
      if (m>1) wstar = (W/B) * integrate(integrand9, 0, (B*w0/W)/(1+(B*w0/W)), n=n, b=b, W=W, B=B)$value /
                               integrate(integrand10, 0, (B*w0/W)/(1+(B*w0/W)), n=n, b=b, W=W, B=B)$value

      b.mean = unlist(lapply(split(x,b.num),mean))
      b.size = table(b.num)
      if (m==1) muhat = b.mean
      else      muhat = (1-wstar)*b.mean + wstar*mu0
      bsqd = (b.size - 1) * (muhat - mu0)^2
      sqd <- (x - muhat[b.num])^2
      B0 = sum(bsqd)
      W0 = sum(sqd)

      #############################################################
      # NOW CONSIDER rho[i] = 1, and do essentially the same thing.

	# first, take out artificial change point if it exists.
	if (temp.cp) {
        if (i<n/2) rho[n-2] <- 0
        else rho[2] <- 0
        temp.cp <- FALSE
      }

      rho[i] = 1
      b.num = 1 + c(0,cumsum(rho[1:(n-1)]))
      bsqd = rep(0,b+1)
      b.mean = unlist(lapply(split(x,b.num),mean))
      b.size = table(b.num)
      if (m==1) muhat = b.mean
      else      muhat = (1-wstar)*b.mean + wstar*mu0
      bsqd = (b.size - 1) * (muhat - mu0)^2
      sqd = (x - muhat[b.num])^2
      B1 = sum(bsqd)
      W1 = sum(sqd)

      #####################################################################################
      # ok, now start doing the real work

      # THE RATIO -- SECTION 4.2
      ratio = (integrate(integrand1, 0, p0, n=n, b=b)$value / integrate(integrand3, 0, p0, n=n, b=b)$value) *
              (W0/W1)^((n-b-2)/2) * (B0/B1)^((b+1)/2) * sqrt(W1/B1) *
              ( integrate(integrand7, 0, (B1*w0/W1)/(1+(B1*w0/W1)), n=n, b=b)$value / 
                integrate(integrand8, 0, (B0*w0/W0)/(1+(B0*w0/W0)), n=n, b=b)$value )      
      p = ratio/(1 + ratio)                 # ratio = p/(1-p) implies this
      if(is.na(p)) {   # For debugging purposes if we get NA.
        print(b)
        print(rho)
        print(B0)
        print(ratio)
        print(i)
      }

      # COMPARE 'p' TO RANDOM VALUE FROM U[0,1] AND UPDATE EVERYTHING.
      rho[i] = 1*(runif(1) < p)
      b.num = 1 + c(0,cumsum(rho[1:(n-1)]))
      b.mean = unlist(lapply(split(x,b.num),mean))
      b.size = table(b.num)
      if (rho[i]==0)  { W=W0; B=B0; b.size[b+1]=0 }
      else       { W=W1; B=B1 }

      # CALCULATE MU-HATS
      wstar = (W/B) * (integrate(integrand9, 0, (B*w0/W)/(1+(B*w0/W)), n=n, b=b, W=W, B=B)$value /
                       integrate(integrand10, 0, (B*w0/W)/(1+(B*w0/W)), n=n, b=b, W=W, B=B)$value)
      muhat = rep(0,b)
      muhat = (1 - wstar)*b.mean + wstar*mu0

      # FILL IN THETAS (and there are n of them)
      theta = muhat[b.num]

    } # End this pass through all n.

    # DONE this pass, save these results:
    results[m,1:as.numeric(dim(results)[2])] = t(theta)
    blocks[m] = b
    rhos[m,1:length(x)] = t(rho)
  } # Done all the passes

    return(list(data=x,
              results=results,
              rhos=rhos,
              blocks=blocks,
              posterior.means=apply(results[burnin:M,1:n],2,mean)))
  
} # End the changepoint() function

