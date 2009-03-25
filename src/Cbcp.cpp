/*
bcp: an R package for performing a Bayesian analysis
of change point problems.

Copyright (C) 2007, 2008 Chandra Erdman and John W. Emerson

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, a copy is available at
http://www.r-project.org/Licenses/

-------------------
FILE: Cbcp.cpp  */

/*  LIBRARIES  */
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <R_ext/Random.h>
#include <vector>

using namespace std;

extern "C" { 

/* MAIN   */
void Cbcp(double *data, 
  int *mcmcreturn,
  int *n,          
  int *burnin,     
  int *mcmc, 
  int *rho,			
  int *rhos,			
  int *blocks,		       
  double *results,		
  double *a,
  double *c,
  double *pmean,				
  double *pvar,		
  double *pchange	
){

  GetRNGstate();

  /* INITIALIZATION */
  int i, m, j, k; 											
  int flag, nn, MM, b, MCMC, BURNIN; 
  int prevb, thisblock, prevblock, myblock;  
  double W, B, W0, W1, B0, B1, trivB, SStot; 			       
  double mu0; 						
  double wstar, ratio, p, p0, w0;
  double xmax1, xmax0, xmax;
  double myrand;
	
  /* SET UP LOCAL COPIES OF VARIABLE FOR CONVENIENCE. */			
  nn = n[0];
  MM = burnin[0] + mcmc[0];	
  p0 = a[0]; 
  if (a[0]==0) p0 = 0.001;
  w0 = c[0]; 
  MCMC = mcmc[0];	
  BURNIN = burnin[0];	

  /* DYNAMIC CREATION OF VECTORS AND MATRICES */								 		       	
  typedef vector<double> DoubleVec;	 
  typedef vector<int> IntVec; 
  DoubleVec bmean(nn); 
  DoubleVec muhat(nn); 
  DoubleVec betai13(nn); 
  DoubleVec ss(nn); 
  IntVec bsize0(nn); 
  DoubleVec bsqd0(nn);	 
  DoubleVec wsqd0(nn); 
  IntVec bsize1(nn); 
  DoubleVec bsqd1(nn); 
  DoubleVec wsqd1(nn); 
  IntVec prevbend(nn); 
  IntVec bend0(nn);	 
  IntVec bend1(nn); 
  DoubleVec cumsums(nn); 
  DoubleVec SScum(nn); 
  DoubleVec sqd(nn); 	 
		
  /* INITIALIZATION OF VARIABLES.------------------------------------------------------------------ */
		
  for(j=0; j<nn; j++) {	
    betai13[j] = 0; 			
    bmean[j] = 0; 						 
    muhat[j] = 0;						
    ss[j] = 0;
    cumsums[j] = 0;					
    SScum[j] = 0;			
    bsize0[j] = 0; 					
    bsize1[j] = 0; 					
    prevbend[j] = 0; 						
    bend0[j] = 0; 					
    bend1[j] = 0; 													
    bsqd0[j] = 0; 					
    bsqd1[j] = 0; 					
    wsqd0[j] = 0; 					
    wsqd1[j] = 0; 	
    sqd[j] = 0;	  
  }

  /* CUMULATIVE SUMS AND SUMS OF SQUARES*/
  cumsums[0] = data[0];						
  SScum[0] = data[0]*data[0];				
  for(j=1; j<nn; j++) {						
    cumsums[j] +=  data[j]+cumsums[j-1];		
    SScum[j] +=  data[j]*data[j] + SScum[j-1];	
  }									
  mu0 = 0;								
  for(j=0; j<nn; j++) mu0 +=data[j];
  mu0 /= (double) nn;	  
  SStot = mu0*mu0*nn*nn;	

  /* TRIVIAL B0 */
  trivB = (data[0] - mu0)*(data[0] - mu0);	
  for(j=0; j<nn; j++) {
    sqd[j] =  (data[j] - mu0)*(data[j] - mu0);
    if (sqd[j] < trivB) trivB = sqd[j];
  }
  if (trivB < 0.00001) trivB = 0.00001;
  
  W1 = 0; B1 = 1; wstar = 0; xmax = 0; myblock=0; prevblock=0; thisblock = 0; prevb=1;	

  /* EVALUATE BETA INTEGRALS THAT ONLY DEPEND ON b. */
  for(j=1; j<(nn-4); j++) {
    betai13[j] = (exp(lbeta((double) j+1, (double) nn-j))*pbeta((double) p0, (double) j+1, (double) nn-j, 1, 0))/ 
      (exp(lbeta((double) j, (double) nn-j+1))*pbeta((double) p0, (double) j, (double) nn-j+1, 1, 0)); 
  }

  b = 1;  
  flag = 1; 
  bend0[0] = nn-1;
  bsize0[0] = nn;
  bmean[0] = cumsums[bend0[0]] / (double) bsize0[0];
  bsqd0[0] = 0;
  wsqd0[0] = SScum[bend0[0]] - SStot/ (double) nn;
  W0 = wsqd0[0];
  W = W0;
  B0 = trivB;
  B = B0; 

  /* -----------------------------------------------------START THE BIG LOOP-------------------------------------------------- */
  for(m=0; m<MM; m++) {

    for(i=0; i<=(nn-2); i++) { 
      
      if ((m > 0) && (prevbend[prevblock] <= i)) prevblock = prevblock+1;

      /* ----------------------------------- CONSIDER REMOVING A CHANGE POINT ----------------------------------*/
      if (rho[i]==1) { 	

	if (b==2) flag = 1;

	if (flag==1) {
	  bend0[0] = nn-1;
	  bsize0[0] = nn;
	  bmean[0] = cumsums[bend0[0]] / (double) bsize0[0];
	  bsqd0[0] = 0;
	  wsqd0[0] = SScum[bend0[0]] - SStot/nn;
	  W0 = wsqd0[0];
	  W = W0;
	  B0 = trivB;
	  B = B0;
	} else {					
	  /* bend0 */
	  bend0[myblock] = bend0[myblock + 1];
	  if (m > 0) bend0[myblock] = prevbend[prevblock];				
	  
	  /* bsize0 */
	  if (myblock==0) bsize0[myblock] = bend0[myblock + 1] + 1;			
	  else {
	    bsize0[myblock] = bend0[myblock] - bend0[myblock - 1];									
	    bsize0[myblock + 1] = bend0[myblock + 1] - bend0[myblock];		
	  }

	  /* bsqd0 and wsqd0 */
	  if (myblock==0) { 
	    bsqd0[myblock] = (cumsums[bend0[myblock]]*cumsums[bend0[myblock]])/ (double) bsize0[myblock];
	    wsqd0[myblock] = SScum[bend0[myblock]] - bsqd0[myblock];						
	  } else {
	    bsqd0[myblock] = ((cumsums[bend0[myblock]] - cumsums[bend0[myblock - 1]])
			      *(cumsums[bend0[myblock]] - cumsums[bend0[myblock - 1]])
			      / (double) bsize0[myblock]);
	    wsqd0[myblock] = (SScum[bend0[myblock]] - SScum[bend0[myblock - 1]]) - bsqd0[myblock];
	  }	  
	  W0 = W + wsqd0[myblock] - wsqd1[myblock] - wsqd1[myblock + 1];
	  B0 = B + bsqd0[myblock] - bsqd1[myblock] - bsqd1[myblock + 1];
	}
      } /* end if rho[i]==1 */			
      /* ----------------------------------- CONSIDER ADDING A CHANGE POINT ---------------------------------*/	
      if (rho[i]==0) { 

	if (b==1) flag = 1; 		
	
	bend1[myblock] = i;													
	bend1[myblock + 1] = bend0[myblock];				
	
	/* bsize1 */
	if (myblock==0) bsize1[myblock] = i + 1;					
	else bsize1[myblock] = i - bend1[myblock - 1];								
	bsize1[myblock + 1] = bend0[myblock] - bend1[myblock];						

	/* bsqd1 and wsqd1 */
	if (myblock==0) {													
	  bsqd1[myblock] = (cumsums[i]*cumsums[i])/ (double) bsize1[myblock];
	  bsqd1[myblock + 1] = (cumsums[bend1[myblock + 1]] - cumsums[bend1[myblock]])
	    *(cumsums[bend1[myblock + 1]] - cumsums[bend1[myblock]])
	    / (double) bsize1[myblock + 1];
	  wsqd1[myblock] = SScum[i] - bsqd1[myblock];
	  wsqd1[myblock + 1] = (SScum[bend1[myblock + 1]] - SScum[bend1[myblock]]) - 
	    bsqd1[myblock + 1];  
	} else {
	  bsqd1[myblock] = (cumsums[bend1[myblock]] - cumsums[bend1[myblock - 1]])
	    *(cumsums[bend1[myblock]] - cumsums[bend1[myblock - 1]])
	    / (double) bsize1[myblock];
	  bsqd1[myblock + 1] = (cumsums[bend1[myblock + 1]] - cumsums[bend1[myblock]])
	    *(cumsums[bend1[myblock + 1]] - cumsums[bend1[myblock]])
	    / (double) bsize1[myblock + 1];
	  if (bsize1[myblock]==1) {					       								
	    wsqd1[myblock] = 0;								
	    wsqd1[myblock + 1] = (SScum[bend1[myblock + 1]] - SScum[bend1[myblock]]) - 
	      bsqd1[myblock + 1];
	  } else{
	    wsqd1[myblock] = (SScum[bend1[myblock]] - SScum[bend1[myblock - 1]]) - 
	      bsqd1[myblock];
	    wsqd1[myblock + 1] = (SScum[bend1[myblock + 1]] - SScum[bend1[myblock]]) - 
	      bsqd1[myblock + 1];
	  }
	}
	W1 = W - wsqd0[myblock] + wsqd1[myblock] + wsqd1[myblock + 1];					
	B1 = B - bsqd0[myblock] + bsqd1[myblock] + bsqd1[myblock + 1];		
      } /* end if rho[i]==0 */

      /* --------------------------------------------------------------------------------------------------------------*/	

      xmax1 = (B1*w0/W1)/(1+(B1*w0/W1));
      xmax0 = (B0*w0/W0)/(1+(B0*w0/W0));

      if (rho[i]==1) b = b - 1;		/* since b is # of blocks with rho[i]=0 */					 	
	
      /* THE RATIO */
      ratio = betai13[b]*pow(W0/W1, (double) (nn-b-2)/2) * pow(B0/B1, (double) (b+1)/2) * pow(W1/B1, 0.5)*	
	(exp(lbeta((double) (b+2)/2, (double) (nn-b-3)/2))*pbeta((double) xmax1, (double) (b+2)/2, (double) (nn-b-3)/2, 1, 0))/
	(exp(lbeta((double) (b+1)/2, (double) (nn-b-2)/2))*pbeta((double) xmax0, (double) (b+1)/2, (double) (nn-b-2)/2, 1, 0));
      p = ratio/(1 + ratio);
      if(b>=nn-5) p = 0;  /* since we use nn-b-4 as beta in wstar */

							
      /* COMPARE p TO RANDOM VALUE FROM U[0,1] AND UPDATE EVERYTHING */
      myrand = runif(0.0, 1.0);
				
      if (myrand < p) { /* IF A CHANGE POINT IS ADDED */

	/* RESET FLAG IF NECESSARY */
	if (flag==1) { 						
	  flag = 0;

	  bend0[0] = i;
	  bend0[1] = nn-1;
	  bsize0[0] = i + 1;
	  bsize0[1] = nn - bsize0[0];
	
	  bsqd0[0] = cumsums[bend0[0]]*cumsums[bend0[0]]/ (double) bsize0[0];
	  bsqd0[1] = (cumsums[bend0[1]] - cumsums[bend0[0]])*(cumsums[bend0[1]] - cumsums[bend0[0]])
	    / (double) bsize0[1];							 
	  wsqd0[0] = SScum[bend0[0]] - bsqd0[0];
	  wsqd0[1] = (SScum[bend0[1]] - SScum[bend0[0]]) - bsqd0[1];

	  B = bsqd0[0] + bsqd0[1] - (SStot/nn);
	  W = wsqd0[0] + wsqd0[1];					
	  B0 = B;
	  B1 = B;
	  W0 = W;
	} else {
	  /* UPDATE ALL "0" VARIABLES */
	  bend0[myblock] = bend1[myblock];
	  bend0[myblock + 1] = prevbend[prevblock + 1];				
	  bsize0[myblock] = bsize1[myblock];
	  bsize0[myblock + 1] = bend0[myblock + 1] - bend0[myblock];
	  if (myblock==0) { 
	    bsqd0[myblock] = (cumsums[bend0[myblock]]*cumsums[bend0[myblock]])/ (double) bsize0[myblock];
	    wsqd0[myblock] = SScum[bend0[myblock]] - bsqd0[myblock];						
	  }
	  else{
	  bsqd0[myblock] = ((cumsums[bend0[myblock]] - cumsums[bend1[myblock - 1]])
			    *(cumsums[bend0[myblock]] - cumsums[bend1[myblock - 1]])
			    / (double) bsize0[myblock]);
	  bsqd0[myblock + 1] = ((cumsums[bend0[myblock + 1]] - cumsums[bend1[myblock]])
				*(cumsums[bend0[myblock + 1]] - cumsums[bend1[myblock]])
				/ (double) bsize0[myblock + 1]);				 		
	  wsqd0[myblock] = (SScum[bend0[myblock]] - SScum[bend1[myblock - 1]]) - bsqd0[myblock];
	  wsqd0[myblock + 1] = (SScum[bend0[myblock + 1]] - SScum[bend1[myblock]]) - bsqd0[myblock + 1];
	  }	 
	}
															
	b = b + 1; 									 
	myblock = myblock + 1;	
	rho[i] = 1; 							
	W = W1;
	B = B1;		
	
      } /* end if change point added */
      else {
	rho[i] = 0;
	W = W0;
	B = B0;
      }
      
      if (rho[i+1]==1) {
	W1 = W; 	
	B1 = B;
	//bend1
	bend1[myblock] = i+1;									
	bend1[myblock + 1] = bend0[myblock + 1];					
	
	//bsize1
	if (myblock==0) bsize1[myblock] = i+2;						
	else bsize1[myblock] = bend1[myblock] - bend1[myblock - 1];		
	bsize1[myblock + 1] = bend1[myblock + 1] - bend1[myblock];	
					
	if (m > 0) {
	  bend1[myblock + 1] = prevbend[prevblock + 1];			
	  bsize1[myblock + 1] = bend1[myblock + 1] - bend1[myblock];	
	}
					
	bsqd1[myblock] = bsqd0[myblock];							
	bsqd1[myblock + 1] = (cumsums[bend1[myblock + 1]] - cumsums[bend1[myblock]])
	  *(cumsums[bend1[myblock + 1]] - cumsums[bend1[myblock]])
	  / (double) bsize1[myblock + 1];
	
	wsqd1[myblock] = wsqd0[myblock];						
	wsqd1[myblock + 1] = (SScum[bend1[myblock + 1]] - SScum[bend1[myblock]]) - 
	  bsqd1[myblock + 1];
   		
	bend0[myblock] = bend1[myblock];						
	bend0[myblock + 1] = bend1[myblock + 1];						
      }
      
      if (rho[i]==1) {
	W0 = W;								
	B0 = B;										
	bend0[myblock - 1] = bend1[myblock - 1];					
	bend0[myblock] = bend1[myblock];						
	bend0[myblock + 1] = prevbend[prevblock + 1];			
	bsize0[myblock - 1] = bsize1[myblock - 1];					
	bsize0[myblock] = bsize1[myblock];						
	bsize0[myblock + 1] = bend0[myblock + 1] - bend0[myblock];		
	bsqd0[myblock - 1] = bsqd1[myblock - 1];					
	bsqd0[myblock] = bsqd1[myblock];							
	bsqd0[myblock + 1] = bsqd1[myblock + 1];					
	wsqd0[myblock - 1] = wsqd1[myblock - 1];					
	wsqd0[myblock] = wsqd1[myblock];					
	wsqd0[myblock + 1] = wsqd1[myblock + 1];					
	
	if (rho[i+1]==1) {
	  bend0[myblock] = prevbend[prevblock + 1];								
	  bend0[myblock + 1] = prevbend[prevblock + 2];							
	  bsize0[myblock] = bend0[myblock] - bend0[myblock - 1];					
	  bsize0[myblock + 1] = bend0[myblock + 1] - bend0[myblock];				
	  bsqd0[myblock] = ((cumsums[bend0[myblock]] - cumsums[bend1[myblock - 1]])
			    *(cumsums[bend0[myblock]] - cumsums[bend1[myblock - 1]])
			    / (double) bsize0[myblock]);
	  wsqd0[myblock] = (SScum[bend0[myblock]] - SScum[bend1[myblock - 1]]) - 
	    ((cumsums[bend0[myblock]] - cumsums[bend1[myblock - 1]])
	     *(cumsums[bend0[myblock]] - cumsums[bend1[myblock - 1]])
	     / (double) bsize0[myblock]);
	  bsqd1[myblock] = (cumsums[bend1[myblock]] - cumsums[bend1[myblock - 1]])		
	    *(cumsums[bend1[myblock]] - cumsums[bend1[myblock - 1]])
	    / (double) bsize1[myblock];
	  bsqd1[myblock + 1] = (cumsums[bend1[myblock + 1]] - cumsums[bend1[myblock]])
	    *(cumsums[bend1[myblock + 1]] - cumsums[bend1[myblock]])
	    / (double) bsize1[myblock + 1];
	  if (bsize1[myblock]==1) {				    									
	    wsqd1[myblock] = 0;								
	    wsqd1[myblock + 1] = (SScum[bend1[myblock + 1]] - SScum[bend1[myblock]]) - 
	      bsqd1[myblock + 1];
	  }
	  else{
	    wsqd1[myblock] = (SScum[bend1[myblock]] - SScum[bend1[myblock - 1]]) - 
	      bsqd1[myblock];
	    wsqd1[myblock + 1] = (SScum[bend1[myblock + 1]] - SScum[bend1[myblock]]) - 
	      bsqd1[myblock + 1];
	  }					   											
	}
      }													

      xmax1 = (B1*w0/W1)/(1+(B1*w0/W1));
      xmax0 = (B0*w0/W0)/(1+(B0*w0/W0));
      
    } /* end this pass through all n. */

    /* RESET FOR NEXT ITERATION */
    myblock = 0; 													
    prevblock = 0; 													
    if (rho[nn-2]==1) for (k=0; k<=b; k++) {
	prevbend[k] = bend1[k];
	bend0[k] = bend1[k];
	bsize0[k] = bsize1[k];
	bsqd0[k] = bsqd1[k];
	wsqd0[k] = wsqd1[k];
    } else {
      for (k=0; k<b; k++) {
	prevbend[k] = bend0[k];
	bend1[k] = bend0[k];
	bsize1[k] = bsize0[k];
	bsqd1[k] = bsqd0[k];
	wsqd1[k] = wsqd0[k];
      }
      bend1[b + 1] = 0;
      bsize1[b] = 0;
      prevbend[b] = nn-1;
    }
    
    if ((b+1) < prevb) {
      for (k = (b+1); k<(nn-2); k++) {
	prevbend[k] = nn-1;
	bend0[k] = 0;
	bsize0[k] = 0;
	bsqd0[k] = 0;
	wsqd0[k] = 0;
	bend1[k] = 0;
	bsize1[k] = 0;
	bsqd1[k] = 0;
	wsqd1[k] = 0;
      }
    }

    prevb = b;

    /* CALCULATE BLOCK MEANS */       
    bmean[0] = cumsums[bend0[0]] / (double) bsize0[0];
    if (b>1) {				
      for (k=1; k<b; k++) {
        bmean[k] = (cumsums[bend0[k]] - cumsums[bend0[k-1]])/ (double) bsize0[k];
      }
    }   

    /* GET xmax */				 
    xmax = (B*w0/W)/(1+(B*w0/W));
    
    /* CALCULATE WSTAR */
    wstar = (W/B)*(exp(lbeta((double) (b+3)/2, (double) (nn-b-4)/2))*pbeta( (double) xmax, (double) (b+3)/2, (double) (nn-b-4)/2,1, 0))/
      (exp(lbeta((double) (b+1)/2, (double) (nn-b-2)/2))*pbeta((double) xmax, (double) (b+1)/2, (double) (nn-b-2)/2, 1, 0));
    if ((wstar<=0) || (wstar>=1)) printf("%d %d %f %f the sky has fallen!!! %20.15f\n", m, i, W, B, wstar);
    
    /* MUHATS */
    for(j=0; j<nn; j++) {
      muhat[j] = (1 - wstar)*bmean[thisblock] + wstar*mu0; 
      if (bend0[thisblock] <= j) thisblock = thisblock+1;
    }
    thisblock = 0; 

    
    /* STORE RESULTS */
    blocks[m] = b;
    if (mcmcreturn[0]==1) {		
      for (j=0; j<nn; j++) {
	rhos[nn*m + j ] = rho[j];
	results[nn*m + j] = muhat[j];
      }
    }

    if (m >= BURNIN) {					
      for (j=0; j<nn; j++) {				
	pchange[j] += rho[j];			
	pmean[j] += muhat[j];			
	ss[j] += muhat[j]*muhat[j];		
      }
    }

  } /* done all iterations. */

  for (j=0; j<nn; j++) {									
    pchange[j] = pchange[j]/MCMC;						
    pmean[j] = pmean[j]/MCMC;							
    pvar[j] = (ss[j]/MCMC - pmean[j]*pmean[j])*(nn/(nn-1)); 	
    if (pvar[j] < 0) pvar[j] = 0;
  }

  PutRNGstate();
		
}  /* END MAIN  */

}
