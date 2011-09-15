/*
bcp: an R package for performing a Bayesian analysis
of change point problems.

Copyright (C) 2011 Chandra Erdman and John W. Emerson

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
#include <Rcpp.h>                               // Are all these needed???
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <R_ext/Random.h>
#include <R.h>
#include <Rdefines.h>
#include <vector>

using namespace std;
using namespace Rcpp;

RcppExport SEXP rcpp_bcp(SEXP pdata, SEXP pmcmcreturn, SEXP pburnin,
                         SEXP pmcmc, SEXP pa, SEXP pc) {

  NumericMatrix data(pdata);
  int nn = data.nrow();
  int mm = data.ncol();

  int mcmcreturn = INTEGER_DATA(pmcmcreturn)[0];
  int burnin = INTEGER_DATA(pburnin)[0];
  int mcmc = INTEGER_DATA(pmcmc)[0];
  int nn2 = nn;
  int MM2 = burnin + mcmc;
  if (mcmcreturn==0) { MM2 = 1; nn2 = 1; }

  double w0 = NUMERIC_DATA(pc)[0];
  double p0 = NUMERIC_DATA(pa)[0];
  if (p0==0) p0 = 0.001;

  // INITIALIZATION OF LOCAL VARIABLES
  int i, m, j, k, l;
  int flag, MM, b;
  int prevb, thisblock, prevblock, myblock;
  double W, B, W0, W1, B0, B1, trivB, SStot;
  double mu0;
  double wstar, ratio, p;
  double xmax1, xmax0, xmax;
  double myrand;

  MM = burnin + mcmc;

  // Things to be returned to R:
  NumericMatrix pmean(nn, mm);
  NumericMatrix pvar(nn, mm);
  NumericVector pchange(nn);
  NumericVector blocks(burnin+mcmc);
  NumericMatrix rhos(nn2, MM2);
  NumericMatrix results(nn2*MM2,mm);

  GetRNGstate();                      // Consider Dirk's comment on this.

  /* DYNAMIC CREATION OF LOCAL VECTORS AND MATRICES */
  typedef vector<double> DoubleVec;
  typedef vector<int> IntVec;
  typedef vector<DoubleVec> DoubleMatrix;

  //double bmean[mm][nn];         /* why these this way? */
  //double bsqd0[mm][nn];
  //double wsqd0[mm][nn];
  //double bsqd1[mm][nn];
  //double wsqd1[mm][nn];
  //double cumsums[mm][nn];
  //double SScum[mm][nn];
  //double sqd[mm][nn];
  //double muhat[mm][nn];
  //double ss[mm][nn];

  DoubleVec cumsumsc(nn);
  DoubleMatrix cumsums(mm, cumsumsc);
  DoubleVec SScumc(nn);
  DoubleMatrix SScum(mm, SScumc);
  DoubleVec sqdc(nn);
  DoubleMatrix sqd(mm, sqdc);
  DoubleVec muhatc(nn);
  DoubleMatrix muhat(mm, muhatc);
  DoubleVec ssc(nn);
  DoubleMatrix ss(mm, ssc);

  DoubleVec bmeanc(nn);
  DoubleMatrix bmean(mm, bmeanc);
  DoubleVec bsqd0c(nn);
  DoubleMatrix bsqd0(mm, bsqd0c);
  DoubleVec wsqd0c(nn);
  DoubleMatrix wsqd0(mm, wsqd0c);
  DoubleVec bsqd1c(nn);
  DoubleMatrix bsqd1(mm, bsqd1c);
  DoubleVec wsqd1c(nn);
  DoubleMatrix wsqd1(mm, wsqd1c);

  IntVec bsize0(nn);
  IntVec bsize1(nn);
  IntVec prevbend(nn);
  IntVec bend0(nn);
  IntVec bend1(nn);
  IntVec rho(nn);					
  DoubleVec betai13(nn);
  DoubleVec mu0vec(mm);
  DoubleVec SStotvec(mm);
  DoubleVec W0vec(mm);
  DoubleVec B0vec(mm);
  DoubleVec W1vec(mm);
  DoubleVec B1vec(mm);
  DoubleVec Wvec(mm);
  DoubleVec Bvec(mm);
  DoubleVec trivBvec(mm);

  /* INITIALIZATION OF VARIABLES.-------------------- */
  for (j=0; j<nn; j++) {
    betai13[j] = 0;
    bsize0[j] = 0;
    bsize1[j] = 0;
    prevbend[j] = 0;
    bend0[j] = 0;
    bend1[j] = 0;
    rho[j] = 0;		
    for (l=0; l<mm; l++) {
      bmean[l][j] = 0;
      cumsums[l][j]  = 0;
      SScum[l][j]  = 0;
      bsqd0[l][j] = 0;
      bsqd1[l][j] = 0;
      wsqd0[l][j]  = 0;
      wsqd1[l][j]  = 0;
      sqd[l][j]  = 0;
      muhat[l][j] = 0;
      ss[l][j] = 0;
    }
  }
  rho[nn-1] = 1;							

  for (l=0; l<mm; l++) {
    mu0vec[l] = 0;
    SStotvec[l] = 0;
    W0vec[l] = 0;
    B0vec[l] = 0;
    W1vec[l] = 0;
    B1vec[l] = 0;
    Wvec[l] = 0;
    Bvec[l] = 0;
    trivBvec[l] = 0;
  }

  /* CUMULATIVE SUMS AND SUMS OF SQUARES*/
  for (l=0; l<mm; l++) {
    cumsums[l][0] = data(0,l);
    SScum[l][0] = data(0,l) * data(0,l);
  }
  for (j=1; j<nn; j++) {
    for (l=0; l<mm; l++) {
      cumsums[l][j] += data(j,l) + cumsums[l][j-1];
      SScum[l][j] += data(j,l)*data(j,l) + SScum[l][j-1];
    }
  }
  mu0 = 0;
  SStot = 0;
  for (j=0; j<nn; j++) for(l=0; l<mm; l++) mu0vec[l] += data(j,l);
  for (l=0; l<mm; l++) {
    mu0vec[l] /= (double) nn;
    SStotvec[l] = mu0vec[l]*mu0vec[l]*nn*nn;
    mu0 += mu0vec[l];
    SStot += SStotvec[l];
  }
  mu0 = mu0/mm;
  SStot = SStot/mm;

  /* TRIVIAL B0 */
  trivB = 0;
  for (l=0; l<mm; l++) {
    trivBvec[l] = (data(0,l) - mu0vec[l])*(data(0,l) - mu0vec[l]);
    for (j=0; j<nn; j++) {
      sqd[l][j] = (data(j,l) - mu0vec[l])*(data(j,l) - mu0vec[l]);
      if (sqd[l][j] < trivBvec[l]) trivBvec[l] = sqd[l][j];
    }
    if (trivBvec[l] < 0.00001) trivBvec[l] = 0.00001;
    trivB += trivBvec[l];
  }

  W1 = 0; B1 = 0; wstar = 0; xmax = 0; myblock=0; prevblock=0; thisblock = 0; prevb=1;
  W0 = 0;

  /* EVALUATE BETA INTEGRALS THAT ONLY DEPEND ON b. */
  for (j=1; j<(nn-4); j++) {
    betai13[j] = ( exp(Rf_lbeta((double) j+1, (double) nn-j)) *
                   Rf_pbeta((double) p0, (double) j+1, (double) nn-j, 1, 0) ) /
                 ( exp(Rf_lbeta((double) j, (double) nn-j+1)) *
                   Rf_pbeta((double) p0, (double) j, (double) nn-j+1, 1, 0) );
  }

  b = 1;
  flag = 1;
  bend0[0] = nn-1;
  bsize0[0] = nn;
  for (l=0; l<mm; l++) {
    bmean[l][0] = cumsums[l][bend0[0]] / (double) bsize0[0];
    bsqd0[l][0] = 0;
    wsqd0[l][0] = SScum[l][bend0[0]] - SStotvec[l]/ (double) nn;
    W0vec[l] = wsqd0[l][0];
    Wvec[l] = W0vec[l];
    W0 += Wvec[l];
  }
  W = W0;
  B0 = trivB;
  B = B0;

  /* ----------------START THE BIG LOOP--------------------------- */
  for (m=0; m<MM; m++) {

    for (i=0; i<=(nn-2); i++) {

      if ((m > 0) && (prevbend[prevblock] <= i)) prevblock = prevblock + 1;

      /* --------- CONSIDER REMOVING A CHANGE POINT ---------------------*/
      if (rho[i]==1) {
	W0 = 0;
	B0 = 0;
	if (b==2) flag = 1;

	if (flag==1) {
	  bend0[0] = nn-1;
          bsize0[0] = nn;
	  for (l=0; l<mm; l++) {
	    bmean[l][0] = cumsums[l][bend0[0]] / (double) bsize0[0];
	    bsqd0[l][0] = 0;
	    wsqd0[l][0] = SScum[l][bend0[0]] - SStotvec[l]/ (double) nn;
            W0vec[l] = wsqd0[l][0];
            W0 += W0vec[l];
          }
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
	  for (l=0; l<mm; l++) {
	    if (myblock==0) {
	      bsqd0[l][myblock] = cumsums[l][bend0[myblock]] * cumsums[l][bend0[myblock]] /
                                  (double) bsize0[myblock];
	      wsqd0[l][myblock] = SScum[l][bend0[myblock]] - bsqd0[l][myblock];
	    } else {
	      bsqd0[l][myblock] = (cumsums[l][bend0[myblock]] - cumsums[l][bend0[myblock - 1]]) *
			          (cumsums[l][bend0[myblock]] - cumsums[l][bend0[myblock - 1]]) /
			          (double) bsize0[myblock];
	      wsqd0[l][myblock] = SScum[l][bend0[myblock]] - SScum[l][bend0[myblock - 1]] - 
                                  bsqd0[l][myblock];
	    }	  
	    W0vec[l] = Wvec[l] + wsqd0[l][myblock] - wsqd1[l][myblock] - wsqd1[l][myblock + 1];
	    B0vec[l] = Bvec[l] + bsqd0[l][myblock] - bsqd1[l][myblock] - bsqd1[l][myblock + 1];
	    W0 += W0vec[l];
	    B0 += B0vec[l];
	  }
        }
      } /* end if rho[i]==1 */
      /* ------------ CONSIDER ADDING A CHANGE POINT ---------------*/
      if (rho[i]==0) {
	B1 = 0;
	W1 = 0;
	if (b==1) flag = 1;
	
	bend1[myblock] = i;
	bend1[myblock + 1] = bend0[myblock];
	
	/* bsize1 */
	if (myblock==0) bsize1[myblock] = i + 1;
	else bsize1[myblock] = i - bend1[myblock - 1];
	bsize1[myblock + 1] = bend0[myblock] - bend1[myblock];

        for (l=0; l<mm; l++) {
	  /* bsqd1 and wsqd1 */
	  if (myblock==0) {
	    bsqd1[l][myblock] = (cumsums[l][i]*cumsums[l][i])/ (double) bsize1[myblock];
	    bsqd1[l][myblock + 1] = (cumsums[l][bend1[myblock + 1]] - cumsums[l][bend1[myblock]])
	      *(cumsums[l][bend1[myblock + 1]] - cumsums[l][bend1[myblock]])
	      / (double) bsize1[myblock + 1];
	    wsqd1[l][myblock] = SScum[l][i] - bsqd1[l][myblock];
	    wsqd1[l][myblock + 1] = (SScum[l][bend1[myblock + 1]] - SScum[l][bend1[myblock]]) -
	      bsqd1[l][myblock + 1];
	  } else {
	    bsqd1[l][myblock] = (cumsums[l][bend1[myblock]] - cumsums[l][bend1[myblock - 1]])
	      *(cumsums[l][bend1[myblock]] - cumsums[l][bend1[myblock - 1]])
	      / (double) bsize1[myblock];
	    bsqd1[l][myblock + 1] = (cumsums[l][bend1[myblock + 1]] - cumsums[l][bend1[myblock]])
	      *(cumsums[l][bend1[myblock + 1]] - cumsums[l][bend1[myblock]])
	      / (double) bsize1[myblock + 1];
	    if (bsize1[myblock]==1) {
	      wsqd1[l][myblock] = 0;
	      wsqd1[l][myblock + 1] = (SScum[l][bend1[myblock + 1]] - SScum[l][bend1[myblock]]) -
	        bsqd1[l][myblock + 1];
	    } else{
	      wsqd1[l][myblock] = (SScum[l][bend1[myblock]] - SScum[l][bend1[myblock - 1]]) -
	        bsqd1[l][myblock];
	      wsqd1[l][myblock + 1] = (SScum[l][bend1[myblock + 1]] - SScum[l][bend1[myblock]]) -
	        bsqd1[l][myblock + 1];
	    }
	  }
	  W1vec[l] = Wvec[l] - wsqd0[l][myblock] + wsqd1[l][myblock] + wsqd1[l][myblock + 1];
	  B1vec[l] = Bvec[l] - bsqd0[l][myblock] + bsqd1[l][myblock] + bsqd1[l][myblock + 1];
	  W1 += W1vec[l];
	  B1 += B1vec[l];
        }
      } /* end if rho[i]==0 */

      /* --------------------------------------------------------*/

      xmax1 = (B1*w0/W1)/(1+(B1*w0/W1));
      xmax0 = (B0*w0/W0)/(1+(B0*w0/W0));

      if (rho[i]==1) b = b - 1;	/* since b is # of blocks with rho[i]=0 */

      /* THE RATIO */
      ratio = betai13[b] * pow(W0/W1, (double) (nn-b-2)/2) * 
              pow(B0/B1, (double) (b+1)/2) * pow(W1/B1, 0.5) *	
	      ( exp(Rf_lbeta((double) (b+2)/2, (double) (nn-b-3)/2)) *
                Rf_pbeta((double) xmax1, (double) (b+2)/2, (double) (nn-b-3)/2, 1, 0) ) /
	      ( exp(Rf_lbeta((double) (b+1)/2, (double) (nn-b-2)/2)) *
                Rf_pbeta((double) xmax0, (double) (b+1)/2, (double) (nn-b-2)/2, 1, 0) );
      p = ratio/(1 + ratio);
      if (b>=nn-5) p = 0;  /* since we use nn-b-4 as beta in wstar */

      /* COMPARE p TO RANDOM VALUE FROM U[0,1] AND UPDATE EVERYTHING */
      myrand = Rf_runif(0.0, 1.0);

      if (myrand < p) { /* IF A CHANGE POINT IS ADDED */

	/* RESET FLAG IF NECESSARY */
	if (flag==1) {
	  flag = 0;
	  B=0;
	  W=0;
	  bend0[0] = i;
	  bend0[1] = nn - 1;
	  bsize0[0] = i + 1;
	  bsize0[1] = nn - bsize0[0];

          for (l=0; l<mm; l++) {
	    bsqd0[l][0] = cumsums[l][bend0[0]]*cumsums[l][bend0[0]] / (double) bsize0[0];
	    bsqd0[l][1] = (cumsums[l][bend0[1]] - cumsums[l][bend0[0]]) *
                          (cumsums[l][bend0[1]] - cumsums[l][bend0[0]]) / (double) bsize0[1];
	    wsqd0[l][0] = SScum[l][bend0[0]] - bsqd0[l][0];
	    wsqd0[l][1] = (SScum[l][bend0[1]] - SScum[l][bend0[0]]) - bsqd0[l][1];

            Bvec[l] = bsqd0[l][0] + bsqd0[l][1] - (SStotvec[l]/nn);
	    Wvec[l] = wsqd0[l][0] + wsqd0[l][1];
	    B1vec[l] = Bvec[l];
	    W1vec[l] = Wvec[l];
            B += Bvec[l];
            W += Wvec[l];
	  }
	  B0 = B;
	  B1 = B;
	  W0 = W;
	} else {
	  /* UPDATE ALL "0" VARIABLES */
	  bend0[myblock] = bend1[myblock];
	  bend0[myblock + 1] = prevbend[prevblock + 1];
	  bsize0[myblock] = bsize1[myblock];
	  bsize0[myblock + 1] = bend0[myblock + 1] - bend0[myblock];
		
          for (l=0; l<mm; l++) {
	    if (myblock==0) { 
	      bsqd0[l][myblock] = (cumsums[l][bend0[myblock]]*cumsums[l][bend0[myblock]])
                                  / (double) bsize0[myblock];
	      wsqd0[l][myblock] = SScum[l][bend0[myblock]] - bsqd0[l][myblock];
	    }
	    else{
	    bsqd0[l][myblock] = ((cumsums[l][bend0[myblock]] - cumsums[l][bend1[myblock - 1]])
			    *(cumsums[l][bend0[myblock]] - cumsums[l][bend1[myblock - 1]])
			    / (double) bsize0[myblock]);
	    bsqd0[l][myblock + 1] = ((cumsums[l][bend0[myblock + 1]] - cumsums[l][bend1[myblock]])
				* (cumsums[l][bend0[myblock + 1]] - cumsums[l][bend1[myblock]])
				/ (double) bsize0[myblock + 1]);
	    wsqd0[l][myblock] = (SScum[l][bend0[myblock]] - SScum[l][bend1[myblock - 1]])
                                - bsqd0[l][myblock];
	    wsqd0[l][myblock + 1] = (SScum[l][bend0[myblock + 1]] - SScum[l][bend1[myblock]])
                                    - bsqd0[l][myblock + 1];
	    }
	  }
	}
	b = b + 1;
	myblock = myblock + 1;
	rho[i] = 1;
        for (l=0; l<mm; l++) {
	  Wvec[l] = W1vec[l];
	  Bvec[l] = B1vec[l];
	  W0vec[l] = Wvec[l];
	  B0vec[l] = Bvec[l];
        }
        W = W1;
	B = B1;

      } /* end if change point added */
      else {
	rho[i] = 0;
         for (l=0; l<mm; l++) {
	  Wvec[l] = W0vec[l];
	  Bvec[l] = B0vec[l];
	  W1vec[l] = Wvec[l];
	  B1vec[l] = Bvec[l];
        }
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

        for (l=0; l<mm; l++) {
	  bsqd1[l][myblock] = bsqd0[l][myblock];
	  bsqd1[l][myblock + 1] = (cumsums[l][bend1[myblock + 1]] - cumsums[l][bend1[myblock]])
	    *(cumsums[l][bend1[myblock + 1]] - cumsums[l][bend1[myblock]])
	    / (double) bsize1[myblock + 1];
	
	  wsqd1[l][myblock] = wsqd0[l][myblock];
	  wsqd1[l][myblock + 1] = (SScum[l][bend1[myblock + 1]] - SScum[l][bend1[myblock]]) -
	    bsqd1[l][myblock + 1];
        }

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
	      
        for (l=0; l<mm; l++) {      
	  bsqd0[l][myblock - 1] = bsqd1[l][myblock - 1];					
	  bsqd0[l][myblock] = bsqd1[l][myblock];
	  bsqd0[l][myblock + 1] = bsqd1[l][myblock + 1];					
	  wsqd0[l][myblock - 1] = wsqd1[l][myblock - 1];					
	  wsqd0[l][myblock] = wsqd1[l][myblock];					
	  wsqd0[l][myblock + 1] = wsqd1[l][myblock + 1];					
        }

	if (rho[i+1]==1) {
	  bend0[myblock] = prevbend[prevblock + 1];
	  bend0[myblock + 1] = prevbend[prevblock + 2];
	  bsize0[myblock] = bend0[myblock] - bend0[myblock - 1];
	  bsize0[myblock + 1] = bend0[myblock + 1] - bend0[myblock];
		
          for(l=0; l<mm; l++) {		
	    bsqd0[l][myblock] = ((cumsums[l][bend0[myblock]] - cumsums[l][bend1[myblock - 1]])
			    *(cumsums[l][bend0[myblock]] - cumsums[l][bend1[myblock - 1]])
			    / (double) bsize0[myblock]);
	    wsqd0[l][myblock] = (SScum[l][bend0[myblock]] - SScum[l][bend1[myblock - 1]]) - 
	      ((cumsums[l][bend0[myblock]] - cumsums[l][bend1[myblock - 1]])
	       *(cumsums[l][bend0[myblock]] - cumsums[l][bend1[myblock - 1]])
	       / (double) bsize0[myblock]);
	    bsqd1[l][myblock] = (cumsums[l][bend1[myblock]] - cumsums[l][bend1[myblock - 1]])
	      *(cumsums[l][bend1[myblock]] - cumsums[l][bend1[myblock - 1]])
	      / (double) bsize1[myblock];
	    bsqd1[l][myblock + 1] = (cumsums[l][bend1[myblock + 1]] - cumsums[l][bend1[myblock]])
	      *(cumsums[l][bend1[myblock + 1]] - cumsums[l][bend1[myblock]])
	      / (double) bsize1[myblock + 1];
          }  

	  if (bsize1[myblock]==1) {	
            for (l=0; l<mm; l++) {		  
	      wsqd1[l][myblock] = 0;								
	      wsqd1[l][myblock + 1] = (SScum[l][bend1[myblock + 1]] - SScum[l][bend1[myblock]]) - 
	        bsqd1[l][myblock + 1];
            }
	  }
	  else{
            for (l=0; l<mm; l++) {			  
	      wsqd1[l][myblock] = (SScum[l][bend1[myblock]] - SScum[l][bend1[myblock - 1]]) - 
	        bsqd1[l][myblock];
	      wsqd1[l][myblock + 1] = (SScum[l][bend1[myblock + 1]] - SScum[l][bend1[myblock]]) - 
	        bsqd1[l][myblock + 1];
            } 
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
	for(l=0; l<mm; l++) {
	  bsqd0[l][k] = bsqd1[l][k];
	  wsqd0[l][k] = wsqd1[l][k];
	}
    } else {
      for (k=0; k<b; k++) {
	prevbend[k] = bend0[k];
	bend1[k] = bend0[k];
	bsize1[k] = bsize0[k];
	for(l=0; l<mm; l++) {
	  bsqd1[l][k] = bsqd0[l][k];
	  wsqd1[l][k] = wsqd0[l][k];
	}
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
        bend1[k] = 0;
	bsize1[k] = 0;
	for(l=0; l<mm; l++) {
	  bsqd0[l][k] = 0;
	  wsqd0[l][k] = 0;
	  bsqd1[l][k] = 0;
	  wsqd1[l][k] = 0;
	}
      }
    }

    prevb = b;

    /* CALCULATE BLOCK MEANS */  
    for(l=0; l<mm; l++) {   
      bmean[l][0] = cumsums[l][bend0[0]] / (double) bsize0[0];
      if (b>1) {				
        for (k=1; k<b; k++) {
          bmean[l][k] = (cumsums[l][bend0[k]] - cumsums[l][bend0[k-1]]) / 
                        (double) bsize0[k];
        }
      }
    }   

    /* GET xmax */				 
    xmax = (B*w0/W)/(1+(B*w0/W));
    
    /* CALCULATE WSTAR */
    wstar = (W/B) * ( exp(Rf_lbeta((double) (b+3)/2, (double) (nn-b-4)/2)) *
                      Rf_pbeta((double) xmax, (double) (b+3)/2, (double) (nn-b-4)/2,1, 0) ) /
                    ( exp(Rf_lbeta((double) (b+1)/2, (double) (nn-b-2)/2)) *
                      Rf_pbeta((double) xmax, (double) (b+1)/2, (double) (nn-b-2)/2, 1, 0) );
    if ((wstar<=0) || (wstar>=1)) printf("%d %d %f %f the sky has fallen!!! %20.15f\n", m, i, W, B, wstar);

    /* MUHATS */
    for (j=0; j<nn; j++) {
      for (l=0; l<mm; l++) {
        muhat[l][j] = (1 - wstar)*bmean[l][thisblock] + wstar*mu0vec[l];
      }
      if (bend0[thisblock] <= j) thisblock++;
    }
    thisblock = 0; 

    /* STORE RESULTS */
    blocks[m] = b;
    if (mcmcreturn==1) {
      for (j=0; j<nn; j++) {
        rhos(j, m) = rho[j];
        for (l=0; l<mm; l++) {
          results(nn*m+j,l) = muhat[l][j];
        }
      }
    }

    if (m >= burnin) {					
      for (j=0; j<nn; j++) {				
	pchange[j] += (double) rho[j];

        for (l=0; l<mm; l++) {
          ss[l][j] += muhat[l][j] * muhat[l][j];
          pmean(j, l) += muhat[l][j];
        }
      }
    }

  } /* done all iterations. */

  /* New initialization: ss (local C), pvar (in R, same as pmean) */

  for (j=0; j<nn; j++) {									
    pchange[j] = pchange[j] / (double) mcmc;						
    for (l=0; l<mm; l++) {
      pmean(j, l) = pmean(j, l) / (double) mcmc;
      pvar(j, l) = ( ss[l][j]/mcmc - pmean(j, l) * pmean(j, l) ) *
                   (mcmc/(mcmc-1)); 
      if (pvar(j, l) < 0) pvar(j, l) = 0;
    }
  }

  PutRNGstate();

  List z;
  z["pmean"] = pmean;
  z["pvar"] = pvar;
  z["pchange"] = pchange;
  z["blocks"] = blocks;
  z["mcmc.rhos"] = rhos;
  z["mcmc.means"] = results;
    
  return z;

}  /* END MAIN  */

