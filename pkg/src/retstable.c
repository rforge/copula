#include <Rmath.h>

#include "nacopula.h"

/* lang5 will probably be defined in future versions of R :*/
#ifndef lang5
static SEXP lang5(SEXP s, SEXP t, SEXP u, SEXP v, SEXP w)
{
    PROTECT(s);
    s = LCONS(s, list4(t, u, v, w));
    UNPROTECT(1);
    return s;
}
#endif

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
    double V = asReal(V0), alp = asReal(alpha),
	V01 = 0.; /* result */

    /* find m := optimal number of summands */
    double m;
    if(V <= 1)
	m=1.;
    else {
	double fV=floor(V);
	double cV=ceil(V);
	if(pow(exp(-V), 1./cV - 1./fV) <= cV/fV)
	    m= fV;
	else
  	    m= cV;
    }

    /* apply standard rejection m times */
    int k;
    SEXP gamma_ = ScalarReal(pow(cos(M_PI_2*alp)*V/m, 1./alp));
    GetRNGstate();
    for(k= 0; k < m; k++) {
	double V01k, U;
	SEXP one = ScalarReal(1.);
	/* standard rejection */
	do {
	    /* sample from the distribution corresponding to the Laplace-Stieltjes
	      transform exp(-(V_0/m)*t^alpha) */
	    V01k = asReal(eval(lang5(install("rstable1"), one, /* alpha=*/ alpha,
				     /* beta=*/ one, /* gamma= */ gamma_),
			       R_GlobalEnv));
	    U = unif_rand();

	} while (U > exp(- V01k));

	V01 += V01k; /* V01 := sum_{k'=1}^k  X  */
    }
    PutRNGstate();

    return(ScalarReal(V01));
}

