#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/**
* <description>
*
* <details>
* @title Sample 1 variate from an exponentially-tilted Stable distribution
* \tilde{S}(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),V0*Indicator(alpha==1),
* Indicator(alpha!=1);1) with corresponding Laplace-Stieltjes transform
* exp(-V0*((1+t)^alpha-1))
* Note: should be *fast*
* @param n  sample size
* @param alpha parameter
* @return numeric(n) vector
*/
SEXP retstable(SEXP V0, SEXP alpha, SEXP rstable1)
{    
    /*get arguments*/
    double V=asReal(V0);
    double alp=asReal(alpha);
    
    /*declare result*/
    SEXP res = PROTECT(allocVector(REALSXP,1));
    double* V01 = REAL(res);
  
    /*find the optimal number of summands*/
    double m;
  	if(V<=1) 
  	  m=1.0;
  	else{
      double fV=floor(V);
      double cV=ceil(V);
  	  if(pow(exp(-V),1.0/cV-1.0/fV)<=cV/fV){
   		  m=fV;
   		}
  	  else{
  	    m=cV;
  	  }
    }
    
    /*apply standard rejection m times*/
    int k;
    double V01k;
    double U;
    double mygamma=pow(cos(M_PI_2*alp)*V/m,1.0/alp);
    GetRNGstate();
    for(k=0;k<m;k++){
  		/*standard rejection*/
  		do{
  		  /*sample from the distribution corresponding to the Laplace-Stieltjes
  		    transform exp(-(V_0/m)*t^alpha)*/
        V01k=REAL(eval(rstable1,1,alpha=alpha,beta=1,gamma=myamma));
  			U=unif_rand();
  		}
  		while(U>exp(-V01k));
  		V01[0]+=V01k;/*update sum*/
  	}
  	PutRNGstate();
  	UNPROTECT(1);
    
    return(res);
}  	
    
