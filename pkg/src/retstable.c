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
SEXP retstable(SEXP V0, SEXP alpha)
{
    double V01, V=asReal(V0), alpha=asReal(alpha);
    
    /*find the optimal number of summands*/
    int m;
  	if(V<=1) 
  	  m=1;
  	else if(floor(V)*pow(exp(-V),-1.0/floor(V))<=ceil(V)*pow(exp(-V),-1.0/ceil(V)))
  		m=int(floor(V));
  	else 
  	  m=int(ceil(V));

    /*apply standard rejection m times*/
    int k,V01k;
    double U, mygamma=pow(cos((M_PI_2)*alpha)*V0/double(m),1.0/alpha);
    GetRNGstate();
    for(k=0;k<m;k++){
  		//standard rejection
  		do{
  		  /*sample from the distribution corresponding to the Laplace-Stieltjes
  		  /*transform exp(-(V_0/m)*t^alpha)*/
  			V01k=rtstable1(1,alpha,beta=1,mygamma);
  			U=unif_rand();
  		}
  		while(U>exp(-V01k));
  		V01+=V01k;/*update sum*/
  	}
  	PutRNGstate();
    return(V01);
}  	
    
