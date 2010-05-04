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

/* sample S0 ~ S(alpha, 1, (cos(alpha*pi/2))^{1/alpha}, I_{alpha == 1}; 1)
* with Laplace-Stieltjes transform exp(-t^alpha) via the algorithm as described
* in Chambers, Mallows, Stuck (1976); S0 = S(alpha,1) in this reference
*/
double rstable0(double alpha){
	/* S(1, 1, 0, 1; 1) with Laplace-Stieltjes transform exp(-t) */
	if(alpha == 1.) return 1.;

	/* alpha in (0,1) */
	double U = unif_rand();
	double W;
	do { W = exp_rand(); } while (W == 0.);
	return pow(A_(M_PI*U,alpha)/pow(W,1.-alpha),1./alpha);
}

/* sample S ~ S(alpha, 1, gamma, 0; 1)
* with characteristic function exp(-gamma^alpha*|t|^alpha*(1-i*sgn(t)
*                                  *tan(alpha*pi/2))), see Nolan's book for the
* parameterization
*/
double rstable(double alpha, double gamma){
	if(alpha == 1.) return gamma/0.;
	else return (gamma/pow(cos(M_PI_2*alpha),1./alpha))*rstable0(alpha);
}

/* for vectorizing rstable */
void rstable_vec(double S[], const int n, const double alpha, const double gamma){
	if(n >= 1){
	GetRNGstate();
		for(int i=0; i < n; i++) S[i] = rstable(alpha, gamma);
	PutRNGstate();
	}
	return;
}

/**
* <description>
*
* <details>
* @title Sample Stable S(alpha, beta = 1, gamma)
* @param n sample size
* @param alpha parameter
* @param gamma parameter
* @return numeric(n) vector
*/
SEXP rstable_c(SEXP n, SEXP alpha, SEXP gamma)
{
    int nn = asInteger(n);
    double alp = asReal(alpha), gamm = asReal(gamma);
    SEXP res = PROTECT(allocVector(REALSXP, nn));

    rstable_vec(REAL(res), nn, alp, gamm);
    UNPROTECT(1);
    return(res);
}

/* sample St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0)^{1/alpha},
*			  V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
* with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)) via double
* rejection, see Hofert (2010)
*/
void retstable_MH(double *St, const double V0[], double h, double alpha, int n)
{
    double gamma0 = cos(M_PI_2 * alpha);

    GetRNGstate();

    for(int i = 0; i < n; i++) { /* for each of the n required variates */

	/* alpha ==1 => St corresponds to a point mass at V0 with
	* Laplace-Stieltjes transform exp(-V0*t) */
	if(alpha == 1.){
		St[i] = V0[i];
		continue;
	}

	/* find m := optimal number of summands using the asymptotic formula */
	int m;
	if(h == 1.) m = imax2(1, (int)round(V0[i]));
	else m = imax2(1, (int)round(V0[i] * pow(h,alpha)));

	/* apply standard rejection m times and sum of the variates St_k*/
	double gamma = pow(gamma0*V0[i]/m, 1./alpha); /* define gamma */
	St[i] = 0.; /* will be the result after the summation */
	for(int k= 0; k < m; k++) {
	    double St_k, U;
	    /* standard rejection */
	    do {
		/* sample St_k~S(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha},
		*		 (V_0/m)*I_{\alpha = 1}; 1) with
		* Laplace-Stieltjes transform exp(-(V_0/m)*t^alpha) */
		St_k = rstable(alpha, gamma);
		U = unif_rand();

	    } while (U > exp(- St_k));
	    /* on acceptance, St_k ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
	    *				       *V_0/m)^{1/alpha},
	    *				       (V_0/m)*I_{alpha = 1},
	    *				       h*I_{alpha != 1}; 1) with
	    * Laplace-Stieltjes transform exp(-(V_0/m)((h+t)^alpha-h^alpha)) */
	    St[i] += St_k;
	} /* complexity m*exp((V_0/m)*h^alpha) => roughly e*V_0*h^alpha */

    } /* end for(i=0 .. n-1) */

    PutRNGstate();
    return;
}

/**
* <description>
* @title Sample length(V0) variates from exponentially-tilted Stable distribution
*
* <details> call retstable_MH for sampling St ~ \tilde{S}(alpha, 1,
*	(cos(alpha*pi/2)*V_0)^{1/alpha}, V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
* with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha))
*
* @param V0_  numeric(n) {R vector}
* @param h parameter {R numeric(1)}
* @param alpha parameter {R numeric(1)}
* @return R numeric(n) vector
*/
SEXP retstable_MH_c(SEXP V0_, SEXP h, SEXP alpha)
{
    int n = LENGTH(PROTECT(V0_ = coerceVector(V0_, REALSXP)));
    double h_ = asReal(h), alp = asReal(alpha);
    SEXP St = allocVector(REALSXP, n); /* the result */
    PROTECT(St);

    retstable_MH(REAL(St), REAL(V0_), h_, alp, n);

    UNPROTECT(2);
    return(St);
}

/**  Sinc(x) :  sin(x)/x , including the limit at x = 0, fast and accurately
 *
 * @param x any (double precision) number
 * @return  sinc(x)
 * @author Martin Maechler; 28 Apr 2010
 */
double sinc_MM(double x) {
    double ax = fabs(x);
    if(ax < 0.006) {
	if(x == 0.) return 1;
	double x2 = x*x;
	if(ax < 2e-4)
	     return 1. - x2/6.;
	else return 1. - x2/6.*(1 - x2/20.);
    }
    /* else */
    return sin(x)/x;
}

/* to be called from R */
SEXP sinc_c(SEXP x_) {
    int n = LENGTH(PROTECT(x_ = coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n); /* the result */
    double *x = REAL(x_), *r = REAL(r_);

    for(int i=0; i < n; i++)
	r[i] = sinc_MM(x[i]);

    UNPROTECT(1);
    return r_;
}

/* Zolotarev's function to the power 1-alpha
* A_(x,a) :=  A(x, a) ^ (1-a( */

/* the 3-arg. version allows more precision for  alpha ~=~ 1 : */
#define _A_3(_x, _alpha_, _I_alpha)				\
    pow(_I_alpha* sinc_MM(_I_alpha*_x), _I_alpha) *		\
    pow(_alpha_ * sinc_MM(_alpha_ *_x), _alpha_ ) / sinc_MM(_x)

double A_(double x, double alpha) {
  double Ialpha = 1.-alpha;
  return _A_3(x, alpha, Ialpha);
}

/* to be called from R --- see experiments in ../tests/retstable-ex.R */
SEXP A__c(SEXP x_, SEXP alpha, SEXP I_alpha) {
    int n = LENGTH(PROTECT(x_ = coerceVector(x_, REALSXP)));
    double alp = asReal(alpha), I_alp = asReal(I_alpha);
    if(fabs(alp + I_alp - 1.) > 1e-12)
	error("'I_alpha' must be == 1 - alpha more accurately");
    SEXP r_ = allocVector(REALSXP, n); /* the result */
    double *x = REAL(x_), *r = REAL(r_);

    for(int i=0; i < n; i++)
	r[i] = _A_3(x[i], alp, I_alp);

    UNPROTECT(1);
    return r_;
}

/* function B(x)/B(0), see Devroye (2009) */
double BdB0(double x,double alpha) {
    double Ialpha = 1.-alpha;
    double den = pow(sinc_MM(alpha*x),alpha) * pow(sinc_MM(Ialpha*x),Ialpha);
    return sinc_MM(x) / den;
}

/* sample St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0)^{1/alpha},
*			  V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
* with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)) via double
* rejection, see Devroye (2009)
*/
void retstable_LD(double *St, const double V0[], double h, double alpha, int n)
{
	/* variables not depending on V0 : */
    	const double c1 = sqrt(M_PI_2);
    	const double c2 = 2.+c1;
    	double b = (1.-alpha)/alpha;

	for(int i = 0; i < n; i++) { /* for each of the n required variates */

		/* set lambda for our parameterization */
		double V0alpha = pow(V0[i],1./alpha);
		double lambda = h*V0alpha;

		/* apply the algorithm of Devroye (2009) to draw from
		* \tilde{S}(alpha, 1, (cos(alpha*pi/2))^{1/alpha}, I_{alpha = 1},
		* lambda*I_{alpha != 1};1) with Laplace-Stieltjes transform
		* exp(-((lambda+t)^alpha-lambda^alpha)) */
		double gamma = pow(lambda,alpha)*alpha*(1.-alpha);
		double sgamma = sqrt(gamma);
		double c3 = c2* sgamma;
		double xi = 1. + M_SQRT2/M_PI * c3;
		double psi = c3*exp(-gamma*M_PI*M_PI/8.)/M_SQRT_PI;
		double w1 = c1*xi/sgamma;
		double w2 = 2.*M_SQRT_PI * psi;
		double w3 = xi*M_PI;
		double X, c, E;
	 	do {
			double U, z, Z; /* <--- will be the "product" of the
					* inner rejection sample */
			do {
		    	    double W_ = unif_rand(), V = unif_rand();
			    if(gamma >= 1) {
				if(V < w1/(w1+w2)) U = fabs(norm_rand())/sgamma;
				else U = M_PI*(1.-W_*W_);
			    }
			    else{
				if(V < w3/(w2+w3)) U = M_PI*W_;
				else U = M_PI*(1.-W_*W_);
			    }
			    double W = unif_rand();
			    double zeta = sqrt(BdB0(U,alpha));
			    double phi = pow(sgamma+alpha*zeta,1./alpha);
			    z = phi/(phi-pow(sgamma,1./alpha));
			    /* compute rho */
			    double rho = M_PI*exp(-pow(lambda,alpha)*(1.-1. \
				         /(zeta*zeta))) / ((1.+c1)*sgamma/zeta \
				 	 + z);
			    double d = 0.;
			    if(U >= 0 && gamma >= 1) d += xi*exp(-gamma*U*U/2.);
			    if(U > 0 && U < M_PI) d += psi/sqrt(M_PI-U);
			    if(U >= 0 && U <= M_PI && gamma < 1) d += xi;
			    rho *= d;
			    Z = W*rho;
			    /* check rejection condition */
			} while( !(U < M_PI && Z <= 1.));

			double
			    a = pow(A_(U,alpha), 1./(1.-alpha)),
			    m = pow(b*lambda/a,alpha),
			    delta = sqrt(m*alpha/a),
			    a1 = delta*c1,
			    a3 = z/a,
			    s = a1+delta+a3;
			double V_ = unif_rand(), N_ = 0., E_ = 0. /* -Wall */;
			if(V_ < a1/s) {
			    N_ = norm_rand();
			    X = m-delta*fabs(N_);
			} else {
			    if(V_ < delta/s)
				X = m+delta*unif_rand();
			    else {
				E_ = exp_rand();
				X = m+delta+E_*a3;
			    }
			}
			E = -log(Z);
			/* check rejection condition */
			c = a*(X-m)+lambda*(pow(X,-b)-pow(m,-b));
			if(X < m) c -= N_*N_/2.;
			else if(X > m+delta) c -= E_;

	   	} while (!(X >= 0 && c <= E));

		/* transform variates from exp(-((lambda+t)^alpha-lambda^alpha))
		* to those of exp(-V_0((h+t)^alpha-h^alpha)) */
		St[i] = V0alpha / pow(X,b);

	} /* end for(i=0 .. n-1) */
	return;
}

/**
* <description>
* @title Sample length(V0) variates from exponentially-tilted Stable distribution
*
* <details> call retstable_LD for sampling St ~ \tilde{S}(alpha, 1,
*	(cos(alpha*pi/2)*V_0)^{1/alpha}, V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
* with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha))
*
* @param V0_  numeric(n) {R vector}
* @param h parameter {R numeric(1)}
* @param alpha parameter {R numeric(1)}
* @return R numeric(n) vector
*/
SEXP retstable_LD_c(SEXP V0_, SEXP h, SEXP alpha)
{
    int n = LENGTH(PROTECT(V0_ = coerceVector(V0_, REALSXP)));
    double h_ = asReal(h);
    double alp = asReal(alpha);
    SEXP St = allocVector(REALSXP, n); /* the result */
    PROTECT(St);

    retstable_LD(REAL(St), REAL(V0_), h_, alp, n);

    UNPROTECT(2);
    return(St);
}
