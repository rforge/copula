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
    int n = LENGTH(PROTECT(V0_ = coerceVector(V0_, REALSXP)));
    double alp = asReal(alpha);
    SEXP V_ = allocVector(REALSXP, n); /* the result */
    PROTECT(V_);

    C_retstable(REAL(V_), REAL(V0_), alp, n);
    /* or
       C_retstable_R(REAL(V_), REAL(V0_), alp, n);
    */

    UNPROTECT(2);
    return(V_);
}


/* MM: this is really not useful, as the time is spent in the call to R's
 *     rstable1() !
 * keep it as example / for didactical reasons */
void C_retstable_R(double *V, const double V0[], double alpha, int n)
{
    double c0 = cos(M_PI_2 * alpha);
    SEXP one, alpha_R;
    one     = PROTECT(ScalarReal(1.));
    alpha_R = PROTECT(ScalarReal(alpha));

    GetRNGstate();
    for(int i = 0; i < n; i++) {
	/* find m := optimal number of summands -- now using asymptotic formula :*/
	double m = fmax(1., floor(0.5 + V0[i]));
	/* Previously:
	 * double m;
	 * if(V0[i] <= 1)
	 *     m= 1.;
	 * else {
	 *     double fV= floor(V0[i]), cV= ceil(V0[i]);
	 *     if(pow(exp(-V0[i]), 1./cV - 1./fV) <= cV/fV)
	 * 	m= fV;
	 *     else
	 * 	m= cV;
	 * } */

	/* apply standard rejection m times */
	SEXP gamma_, rstableCall;
	PROTECT(gamma_ = ScalarReal(pow(c0*V0[i]/m, 1./alpha)));
	PROTECT(rstableCall = lang5(install("rstable1"), one, /* alpha=*/ alpha_R,
				    /* beta=*/ one, /* gamma= */ gamma_));
	for(int k= 0; k < m; k++) {
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
    UNPROTECT(1);
    return;
}

/* computes a stable  S(alpha, beta=1, gamma, \I_{\{\alpha=1\}}; 1) distributed random variate
 * with exponent alpha, skewness parameter beta=1, scale parameter gamma, and
 * shift parameter delta = I_{\alpha = 1} --- adapted from gsl_ran_levy_skew */
double S(const double alpha, const double gamma)
{
    double delta=0.;
    double Theta=M_PI*(unif_rand() - 0.5), W;
    do {
	/* computes an variate from an Exp(lambda = 1) - distribution: */
	W = exp_rand();
    } while (W == 0.); /* U = unif_rand() == 1  should never happen */

    if(alpha == 1.) { /* NOTA BENE: Marius had this wrong .. */
#define p2 M_PI_2
	double p2bT = p2 + /* beta * */ Theta;
	return (p2bT*tan(Theta) - /* beta * */ log((p2 * W * cos(Theta))/p2bT))/p2;
#undef p2
    }
    else{
	double t = /* beta * */ tan(M_PI_2*alpha);
	double th0 = atan(t)/alpha;
	double s = pow(1.+t*t, 1./(2.*alpha));
	double x = s * sin(alpha*(Theta+th0)) / pow(cos(Theta), 1./alpha) *
	    pow(cos(Theta-alpha*(Theta+th0)) / W, (1.-alpha)/alpha);
	return delta + gamma*x;
    }
}


/*-- [TODO]  this is really "the same" as Stilde_MH() below ------------------- */
void C_retstable(double *V, const double V0[], double alpha, int n)
{
    double c0 = cos(M_PI_2 * alpha);

    GetRNGstate();
    for(int i = 0; i < n; i++) {
	/* find m := optimal number of summands -- now using asymptotic formula :*/
	int m = imax2(1, (int)round(V0[i]));

	/* apply standard rejection m times */
	double gamma = pow(c0*V0[i]/m, 1./alpha);
	V[i] = 0.;
	for(int k= 0; k < m; k++) {
	    double V01_k, U;
	    /* standard rejection */
	    do {
		/* sample from the distribution corresponding to the Laplace-Stieltjes
		   transform exp(-(V0_0/m)*t^alpha) */
		V01_k = S(alpha, gamma); /* rstable(1, alpha, beta=1, gamma) : */
		U = unif_rand();

	    } while (U > exp(- V01_k));

	    V[i] += V01_k; /* V__ := sum_{k'=1}^k  X  */
	}
    } /* end for(i .. n) */
    PutRNGstate();
    return;
}


/**  Sinc(x) :  sin(x)/x , including the limit at x = 0, fast and accurately
 *
 * @param x any (double precision) number
 * @return  sinc(x)
 * @author Martin Maechler; 28 Apr 2010
 */
double R_sinc(double x) {
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

SEXP sinc(SEXP x_) {
    int n = LENGTH(PROTECT(x_ = coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n); /* the result */
    double *x = REAL(x_), *r = REAL(r_);

    for(int i=0; i < n; i++)
	r[i] = R_sinc(x[i]);

    UNPROTECT(1);
    return r_;
}

// Zolotarev's function to the power 1-alpha -------------------------
//  A_(x,a) == A~(x,a)  :=  A(x, a) ^ (1-a(

// the 3-arg. version allows more precision for  alpha ~=~ 1 :
#define _A_3(_x, _alpha_, _I_alpha)				\
    pow(_I_alpha* R_sinc(_I_alpha*_x), _I_alpha) *		\
    pow(_alpha_ * R_sinc(_alpha_ *_x), _alpha_ ) / R_sinc(_x)

double A_(double x, double alpha) {
  double Ialpha = 1.-alpha;
  return _A_3(x, alpha, Ialpha);
}

// to be called from R --- see experiments in ../tests/retstable-ex.R
SEXP A_Zolotarev(SEXP x_, SEXP alpha, SEXP I_alpha) {
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

// sample S(alpha,1,gamma,0;1), see Nolan's book for the parameterization
// reference: Chambers, Mallows, Stuck (1976)
double S_(double alpha,double gamma) {
    if(alpha == 1.) return DBL_MAX;
    double W;
    do { W = exp_rand(); } while (W == 0.);
    double U = unif_rand();// real number in [0,1)
    return gamma * pow(A_(M_PI*U, alpha) / (cos(M_PI_2*alpha)*pow(W, 1.-alpha)),
		       1./alpha);
}

void rstable_b1_vec(double V[], const int n,
		    const double alpha, const double gamma)
{
    if(n >= 1) {
	GetRNGstate();
	for(int i=0; i < n; i++)
/* S_() is *WRONG*
 *           V[i] = S_(alpha, gamma);
 */
	    V[i] = S(alpha, gamma);
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
SEXP rstable_b1(SEXP n, SEXP alpha, SEXP gamma)
{
    int nn = asInteger(n);
    double alp = asReal(alpha), gamm = asReal(gamma);
    SEXP res = PROTECT(allocVector(REALSXP, nn));

    rstable_b1_vec(REAL(res), nn, alp, gamm);
    UNPROTECT(1);
    return(res);
}


// sample  S~  ~ Stilde(alpha, beta = 1, gamma = (cos(alpha*pi/2)*V_0)^{1/alpha},
//                                    V_0\I_{\{\alpha = 1\}}, h\I_{\{\alpha\neq1\}}; 1)
// with Laplace-Stieltjes transform exp(-V_0((h+t)^\alpha-h^\alpha))
double Stilde_MH(double V0,double h,double alpha)
{
    if(alpha == 1.)
	return V0;// point mass at V0 with Laplace-Stieltjes transform exp(-V0*t)
    int m = imax2(1, (int)round(V0 * pow(h,alpha)));
    double gamma = pow(cos(M_PI_2*alpha)*V0/m,1./alpha),
	St = 0.;

    for(int i = 0; i < m; i++) {
	double St_k, U;
	do {
	    // sample S~S(alpha,1,(cos(alpha*pi/2)*V_0/m)^{1/alpha},(V_0/m)\I_{\{\alpha = 1\}};1)
	    // with Laplace-Stieltjes transform exp(-(V_0/m)*t^alpha)
	    St_k = S_(alpha,gamma);// beta = 1, delta = 0 (since alpha in (0,1))
	    U = unif_rand();
	}
	while(U > exp(-St_k));// expected number of iterations = exp((V_0/m)*h^alpha)
	St += St_k;
    }// complexity = m*exp((V_0/m)*h^alpha) which is roughly e*V_0*h^alpha
    return St;
}


// function B(x)/B(0), see Devroye (2009)
double BdB0(double x,double alpha) {
    double Ialpha = 1.-alpha;
    double den = pow(R_sinc(alpha*x),alpha) * pow(R_sinc(Ialpha*x),Ialpha);
    return R_sinc(x) / den;
}

// computes a Stilde(alpha,1,(cos(alpha*pi/2)*V_0)^{1/alpha},V_0\I_{\{\alpha = 1\}},h\I_{\{\alpha\neq1\}};1)
// distributed random variate with Laplace-Stieltjes transform exp(-V_0((h+t)^\alpha-h^\alpha))
double Stilde_LD(double V0,double h,double alpha)
{
    double c1 = sqrt(M_PI_2);
    double c2 = 2.+c1;
    double b = (1.-alpha)/alpha;

    // all these depend on V0 :
    double V0alpha = pow(V0,1./alpha);
    double lambda = h*V0alpha;
    // apply algorithm of Devroye (2009) to generate Stilde_{alpha,lambda}
    // ~Stilde(alpha,1,(cos(alpha*pi/2))^{1/alpha},\I_{\{\alpha = 1\}},lambda\I_{\{\alpha\neq1\}};1)
    // with Laplace-Stieltjes transform exp(-((lambda+t)^\alpha-lambda^\alpha))
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
	double U, z, Z; // <--- will be the "product" of the inner rejection sample
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
	    // compute rho
	    double rho = M_PI*exp(-pow(lambda,alpha)*(1.-1./(zeta*zeta))) /
		((1.+c1)*sgamma/zeta + z);
	    double d = 0.;
	    if(U >= 0 && gamma >= 1) d += xi*exp(-gamma*U*U/2.);
	    if(U > 0 && U < M_PI) d += psi/sqrt(M_PI-U);
	    if(U >= 0 && U <= M_PI && gamma < 1) d += xi;
	    rho *= d;
	    Z = W*rho;
	    // check rejection condition
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
	// check rejection condition
	c = a*(X-m)+lambda*(pow(X,-b)-pow(m,-b));
	if(X < m) c -= N_*N_/2.;
	else if(X > m+delta) c -= E_;

    } while (!(X >= 0 && c <= E));

    // return V0a * S~;  S~ = Stilde_{alpha,lambda} = 1 / X^b
    return V0alpha / pow(X,b);
}
