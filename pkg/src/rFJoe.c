#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/**
* <description>
*
* <details>
* @title Sample F with probability mass function p_k=\binom{alpha}{k}(-1)^{k-1},
*   k in IN, for Joe's family
* 	Note: should be *fast* as it is used as building block in many places
* @param n sample size
* @param alpha parameter
* @return numeric(n) vector
*/
SEXP rFJoe(SEXP n, SEXP alpha)
{
    int nn = asInteger(n), i;
    double alp = asReal(alpha);
    /* double* U = R_alloc(nn, sizeof(double)); */
    SEXP res = PROTECT(allocVector(REALSXP, nn));
    double* V = REAL(res);

    if(nn >= 1) {
	if(alp == 1) {
	    for(i=0; i < nn; i++)
		V[i] = 1.;
	} else {
	    double U, G1_a = gammafn(1-alp), I_al = 1./alp;
	    GetRNGstate();

	    for(i=0; i < nn; i++)  {
		/* FIXME(MM): (for alp not too close to 1): re-express
		 * using 1-U !*/
		U = unif_rand();
		if(U <= alp)
		    V[i] = 1.;
		else {
		    double Ginv = pow((1-U)*G1_a, -I_al);
		    double fGinv = floor(Ginv);
		    if(1-U < 1./(fGinv*beta(fGinv,1.-alp)))
			V[i] = ceil(Ginv);
		    else V[i] = fGinv;
		}
	    }
	    PutRNGstate();
	}
    }
    UNPROTECT(1);
    return(res);
}
