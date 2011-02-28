/*
 Copyright (C) 2010 Marius Hofert and Martin Maechler

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 3 of the License, or (at your option) any later 
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <Rmath.h>

#include "nacopula.h"

/**
 * Sample V ~ F with Laplace-Stieltjes transform 
 * (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0))
 * via the algorithm of Hofert (2010). C version.
 * @param p parameter 1-e^(-theta1)
 * @param alpha parameter theta0/theta1 in (0,1]
 * @param iAlpha 1-alpha 
 * @param theta0_le_1 in {0,1} with 1 if and only if theta0 <= 1
 * @return a random variate from F
 * @author Marius Hofert, Martin Maechler
*/
double rFFrank(double p, double alpha, double iAlpha /**< == 1 - alpha */, 
	       int theta0_le_1){
    double U, V;
    if(theta0_le_1) {
	do {
	    U = unif_rand();
	    V = rLog(p);
	} while (U*(V-alpha) > 1./beta(V, iAlpha));

    } else {
	double gamma_1_a = gammafn(iAlpha);
	do {
	    U = unif_rand();
	    V = rFJoe(alpha, iAlpha, gamma_1_a);
	} while(U > pow(p, V-1.));
    }
    return V;
}

/** 
 * Vectorize rFFrank. Generate a vector of variates 
 * V ~ F with Laplace-Stieltjes transform 
 * (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)).
 * @param V vector of random variates from F (result)
 * @param n length of the vector V
 * @param theta_0 parameter theta0 in (0,infinity)
 * @param theta_1 parameter theta1 in [theta0, infinity)
 * @return none
 * @author Marius Hofert, Martin Maechler
*/
void rFFrank_vec(double *V, const int n, const double theta_0, 
		 const double theta_1){
    double p  =  - expm1(-theta_1),
	alpha = theta_0 / theta_1, iAlpha = (theta_1 - theta_0) / theta_1;
    int th_0_le_1 = (theta_0 <= theta_1);

    if(n >= 1) {
	GetRNGstate();

	for(int i=0; i < n; i++)
	    V[i] = rFFrank(p, alpha, iAlpha, th_0_le_1);

	PutRNGstate();
    }
    return;
}

/**
 * Generate a vector of variates V ~ F with Laplace-Stieltjes transform 
 * (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)).
 * @param n sample size
 * @param theta_0 parameter theta0 in (0,infinity)
 * @param theta_1 parameter theta1 in [theta0, infinity)
 * @return vector of random variates V
 * @author Martin Maechler
*/
SEXP rFFrank_c(SEXP n_, SEXP theta_0_, SEXP theta_1_){
    int n = asInteger(n_);
    double theta_0 = asReal(theta_0_),
	   theta_1 = asReal(theta_1_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    if(n >= 1)
	rFFrank_vec(REAL(res), n, theta_0, theta_1);
    UNPROTECT(1);
    return res;
}
