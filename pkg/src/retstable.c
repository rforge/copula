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

/* MM: this is really not (yet) useful:- CPU time saving is close to 0
 * -- if you had used  Rprof() as I said,  you'd have seen that really only
 *  rstable1() is the costly part ! */

/**
* <description>
* @title Sample  length(V0)  variates from exponentially-tilted Stable distribution
*
* <details> Sample from
* S~(alpha,1, (cos(alpha*pi/2)V0/m)^(1/alpha), V0*I(alpha==1), I(alpha!=1);1)
* with corresponding Laplace-Stieltjes transform
* exp(-V0*((1+t)^alpha-1))

* @param V0_  numeric(n) {R vector}
* @param alpha parameter {R numeric(1)}
* @return R numeric(n) vector
*/
SEXP retstable_c(SEXP V0_, SEXP alpha)
{
    int n = LENGTH(PROTECT(V0_ = coerceVector(V0_, REALSXP))), i;
    double *V0 = REAL(V0_), alp = asReal(alpha), c0 = cos(M_PI_2 * alp);
    SEXP one, V__;
    double *V = REAL(PROTECT(V__ = allocVector(REALSXP, n))); /* the result */

    PROTECT(one = ScalarReal(1.));
    GetRNGstate();
    for(i = 0; i < n; i++) {
	/* find m := optimal number of summands */
	double m;
	if(V0[i] <= 1)
	    m= 1.;
	else {
	    double fV= floor(V0[i]), cV= ceil(V0[i]);
	    if(pow(exp(-V0[i]), 1./cV - 1./fV) <= cV/fV)
		m= fV;
	    else
		m= cV;
	}

	/* apply standard rejection m times */
	SEXP gamma_, rstableCall;
	int k;
	PROTECT(gamma_ = ScalarReal(pow(c0*V0[i]/m, 1./alp)));
	PROTECT(rstableCall = lang5(install("rstable1"), one, /* alpha=*/ alpha,
				    /* beta=*/ one, /* gamma= */ gamma_));
	for(k= 0; k < m; k++) {
	    double V01_k, U;
	    /* standard rejection */
	    do {
		/* sample from the distribution corresponding to the Laplace-Stieltjes
		   transform exp(-(V0_0/m)*t^alpha) */
		V01_k = asReal(eval(rstableCall, R_GlobalEnv));
		U = unif_rand();

	    } while (U > exp(- V01_k));

	    V[i] += V01_k; /* V__ := sum_{k'=1}^k  X  */
	}
	UNPROTECT(2);
    } /* end for(i .. n) */
    PutRNGstate();
    UNPROTECT(3);
    return(V__);
}

