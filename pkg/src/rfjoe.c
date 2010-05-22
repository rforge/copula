#include <Rmath.h>

#include "nacopula.h"

/**
 * Sample V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with 
 * Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the algorithm of 
 * Hofert (2010).
 * Note: The caller of this function must use GetRNGstate() and PutRNGstate().
 * @param alpha parameter theta0/theta1 in (0,1]
 * @param iAlpha 1-alpha 
 * @param gamma_1_a Gamma(1-alpha)
 * @return a random variate from F
 * @author Marius Hofert, Martin Maechler
*/
double rFJoe(double alpha, double iAlpha /**< := 1 - alpha */, 
	     double gamma_1_a /**< == Gamma(1 - alpha) == Gamma(iALpha) */){ 
    double U, I_al = 1./alpha;

    /**< FIXME(MM): (for alpha not too close to 1): re-express using 1-U */
    U = unif_rand();
    if(U <= alpha)
	return 1.;
    else { /**< alpha < U < 1 */
	double Ginv = pow((1-U)*gamma_1_a, -I_al);
	double fGinv = floor(Ginv);
	if(1-U < 1./(fGinv*beta(fGinv, iAlpha)))
	    return ceil(Ginv);
	else return fGinv;
    }
}

/** 
 * Vectorize rFJoe. Generate a vector of variates 
 * V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with Laplace-Stieltjes 
 * transform 1-(1-exp(-t))^alpha.
 * @param V vector of random variates from F (result)
 * @param n length of the vector V
 * @param alpha parameter theta0/theta1 in (0,1]
 * @param iAlpha 1-alpha 
 * @return none
 * @author Marius Hofert, Martin Maechler
*/
void rFJoe_vec(double V[], const int n,
	       const double alpha, const double iAlpha /**< := 1 - alpha */){
    if(n >= 1) {
	double G1_a = gammafn(iAlpha);
	GetRNGstate();

	for(int i=0; i < n; i++)
	    V[i] = rFJoe(alpha, iAlpha, G1_a);

	PutRNGstate();
    }
    return;
}

/**
 * Generate a vector of variates 
 * V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with Laplace-Stieltjes 
 * transform 1-(1-exp(-t))^alpha.
 * Note: Should be fast as it is used as a building block in different places.
 * @param n sample size
 * @param alpha parameter theta0/theta1 in (0,1]
 * @return vector of random variates V
 * @author Martin Maechler
*/
SEXP rFJoe_c(SEXP n, SEXP alpha)
{
    int nn = asInteger(n);
    double alp = asReal(alpha);
    SEXP res = PROTECT(allocVector(REALSXP, nn));

    rFJoe_vec(REAL(res), nn, alp, 1. - alp);

    UNPROTECT(1);
    return(res);
}
