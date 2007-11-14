/*---------------------------------------------------------------------- 
Chandra Erdman, Yale University, Dec. 2006                                                        

C Implementation of                  
"A BAYESIAN ANALYSIS FOR CHANGE POINT PROBLEMS"  
Journal of the American Statistical Association     
Daniel Barry and J. A. Hartigan, 1993                                                          
 
----------------------------------------------------------------------*/
/*  LIBRARIES  */
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <malloc.h>

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
	) {
		
		/* INITIALIZATION */
        
		int i; 							
		int m; 							
		int j; 							
		int k; 							
		int flag, nn, MM, b, b1, cursize, curblock, MCMC, BURNIN;  /* NEW: MCMC AND BURNIN */
		double W, B, W0, W1, B0, B1; 			       
		double mu0; 						
		double wstar, ratio, p, p0, w0;
		double xmax1, xmax2, xmax3, xmax4; 		        
	
		/* SET UP LOCAL COPIES OF VARIABLE FOR CONVENIENCE. */	
		
		nn = n[0];		
		MM = burnin[0] + mcmc[0];	
		p0 = a[0]; 		
		w0 = c[0]; 
		MCMC = mcmc[0];	
		BURNIN = burnin[0];	

		/* DYNAMIC CREATION OF VECTORS AND MATRICES: ---------------------------------------- */
        
		int *bsize = (int *) malloc(nn*sizeof(int)); 		
		int *bnum = (int *) malloc(nn*sizeof(int)); 		
		double *bmean = (double *) malloc(nn*sizeof(double)); 	
		double *bsqd = (double *) malloc(nn*sizeof(double)); 	
		double *sqd = (double *) malloc(nn*sizeof(double)); 	
		double *muhat = (double *) malloc(nn*sizeof(double)); 
		double *ss = (double *) malloc(nn*sizeof(double)); 				
		double *betai13 = (double *) malloc(nn*sizeof(double));	
		
		/* INITIALIZATION OF VARIABLES.------------------------------------------------------------------ */
  
		for(j=0; j<nn; j++) { 
			bsize[j] = 0; 
			bnum[j] = 0;
			bmean[j] = 0; 
			bsqd[j] = 0; 
			sqd[j] = 0; 
			muhat[j] = 0;
			ss[j] = 0;			
			betai13[j] = 0;
        	}	

		W = 0; B = 0; W0 = 0; B0 = 0; W1 = 0; B1 = 0; b = 0; mu0 = 0; wstar = 0;
		xmax1 = 0.2; xmax4 = 0;
		flag = 0;
		
		for(j=0; j<nn; j++) mu0 +=  data[j]/nn;
		
		/* EVALUATE BETA INTEGRALS THAT ONLY DEPEND ON b. */
		for(j=1; j<nn; j++) {
                  betai13[j] = (exp(lbeta((double) j+1, (double) nn-j))*pbeta((double) p0, (double) j+1, (double) nn-j, 1, 0))/ 
				(exp(lbeta((double) j, (double) nn-j+1))*pbeta((double) p0, (double) j, (double) nn-j+1, 1, 0)); 
		}
           
                /* START THE BIG LOOP--------------------------------------------------------------------------------- */
		for(m=0; m<MM; m++) { 			
			for(i=0; i<=(nn-2); i++) { 	
				
				/* CONSIDER rho[i] = 0 */
				rho[i] = 0;
				
				/* COMPUTE BLOCK NUMBER VECTOR */
				curblock = 1;
				for(j=0; j<nn; j++) {
					bnum[j] = curblock;
					curblock = curblock + rho[j];
				}
				b = bnum[nn-1];	
						
				/* INTRODUCE AN ARTIFICIAL CHANGE POINT IF WE HAVE ONLY ONE BLOCK */
				if (b==1) {
					flag = 1;
					if (i < (nn)/2) rho[nn - 2] = 1;	
					else rho[1] = 1;
					
					/* REVISE BLOCK NUMBER VECTOR */
					curblock = 1;
					for(j=0; j<nn; j++) {
						bnum[j] = curblock;
						curblock = curblock + rho[j];
					}
					b = 2;
				}
    
				/* FIND BLOCK SIZES */
				cursize = 0;
				for(j=0; j<nn; j++) {
					cursize++;
					if (rho[j]==1) {
						bsize[bnum[j]-1] = cursize;
						cursize = 0;
					}
				}
				
				/* CALCULATE BLOCK MEANS */
				for(j=0; j<nn; j++) bmean[j] = 0;
				for(j=0; j<nn; j++) bmean[bnum[j]-1] +=  data[j] / (double) bsize[bnum[j]-1];

					
				/* CALCULATE SQUARED DEVIATIONS */
				B0 = 0;
				for(k=0; k<b; k++) {	    
					bsqd[k] = ((double) bsize[k])*(bmean[k] - mu0)*(bmean[k] - mu0);
					B0 += bsqd[k];
				}
				W0 = 0;
				for (j=0; j<nn; j++) { 
					sqd[j] = (data[j] - bmean[bnum[j]-1])*(data[j] - bmean[bnum[j]-1]);
					W0 += sqd[j]; 		  
				}
                                xmax3 = (B0*w0/W0)/(1+(B0*w0/W0));

	/*------------------------------------------------------------------------------------------------------------------------*/
				/* NOW CONSIDER rho[i] = 1 */
				
				rho[i] = 1;
	  
				/* COMPUTE BLOCK NUMBER VECTOR */
				curblock = 1;
				for(j=0; j<nn; j++) {
					bnum[j] = curblock;
					curblock = curblock + rho[j];
				}
				b1 = bnum[nn-1];
				
				/* FIND BLOCK SIZES */
				cursize = 0;		          
				for(j=0; j<nn; j++) { 		
					cursize++;			
					if (rho[j]==1) {			
						bsize[bnum[j]-1] = cursize;    
						cursize = 0;			
					}					
				}					
				
				/* CALCULATE BLOCK MEANS */
				for(j=0; j<nn; j++) bmean[j] =0;
				for(j=0; j<nn; j++) bmean[bnum[j]-1] +=  data[j] / (double) bsize[bnum[j]-1]; 
								
				/* CALCULATE SUMS OF SQUARES */ 
				B1 = 0;									
				for(k=0; k<b1; k++) {	    						
					bsqd[k] = ((double) bsize[k])*(bmean[k] - mu0)*(bmean[k] - mu0);	
					B1 += bsqd[k];								
				}				
				W1 = 0;									
				for (j=0; j<nn; j++) { 							
					sqd[j] = (data[j] - bmean[bnum[j]-1])*(data[j] - bmean[bnum[j]-1]);	
					W1 += sqd[j]; 		  
				}
				xmax2 = (B1*w0/W1)/(1+(B1*w0/W1));
				
				/* THE RATIO */
				ratio = betai13[b]*pow(W0/W1, (double) (nn-b-2)/2) * pow(B0/B1, (double) (b+1)/2) * pow(W1/B1, 0.5)*	
					(exp(lbeta((double) (b+2)/2, (double) (nn-b-3)/2))*pbeta((double) xmax2, (double) (b+2)/2, (double) (nn-b-3)/2, 1, 0))/
					(exp(lbeta((double) (b+1)/2, (double) (nn-b-2)/2))*pbeta((double) xmax3, (double) (b+1)/2, (double) (nn-b-2)/2, 1, 0));
				p = ratio/(1 + ratio);
                                if(b>=nn-5) p = 0;		
  
				/* COMPARE p TO RANDOM VALUE FROM U[0,1] AND UPDATE EVERYTHING */
				if (runif(0.0, 1.0) < p) rho[i] = 1; else rho[i] = 0;		
                                
				/* RESET FLAG IF NECESSARY */
				if (flag==1) { 	
					if (i < (nn)/2) rho[nn - 2] = 0;	
					else  rho[1] = 0;
					flag = 0;
				}
				
				} /* end this pass through all n. */
				
				/* COMPUTE BLOCK NUMBER VECTOR */ 				 
				curblock = 1;
				for(j=0; j<nn; j++) {
					bnum[j] = curblock;
					curblock = curblock + rho[j];
				}
				b = bnum[nn-1];
				
				/* FIND BLOCK SIZES */
				cursize = 0;
				for(j=0; j<nn; j++) bsize[j] = 0;
				for(j=0; j<nn; j++) { 		
					cursize++;			
					if (rho[j]==1) {			
						bsize[bnum[j]-1] = cursize;    
						cursize = 0;			
					}					
				}					
				
				/* CALCULATE BLOCK MEANS */ 
				for(j=0; j<nn; j++) bmean[j] =0;
				for(j=0; j<nn; j++) bmean[bnum[j]-1] +=  data[j] / ((double) bsize[bnum[j]-1]); 	
				
				/* GET xmax4 */				 
				if (rho[nn-2]==0) {W = W0; B = B0;}
				else {W = W1; B = B1;}
				xmax4 = (B*w0/W)/(1+(B*w0/W));
								
				/* CALCULATE WSTAR */
				wstar = (W/B)*(exp(lbeta((double) (b+3)/2, (double) (nn-b-4)/2))*pbeta( (double) xmax4, (double) (b+3)/2, (double) (nn-b-4)/2,1, 0))/
                                              (exp(lbeta((double) (b+1)/2, (double) (nn-b-2)/2))*pbeta((double) xmax4, (double) (b+1)/2, (double) (nn-b-2)/2, 1, 0));
				if (wstar<=0 | wstar>=1) printf("%d %d %f %f the sky has fallen!!! %20.15f\n", m, i, W, B, wstar);
									
				for(j=0; j<nn; j++) muhat[j] = (1 - wstar)*bmean[bnum[j]-1] + wstar*mu0;  
				
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
		}

		free(bsize);
		free(bnum); 
		free(bmean); 
		free(bsqd); 
		free(sqd); 
		free(muhat); 
		free(ss);		
		free(betai13);
		
 	}  /* END MAIN  */

	
