/*#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2008, 2009
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################*/


/*****************************************************************************

  Multivariate serial independence test based on the empirical
  copula process

  Ivan Kojadinovic, December 2007

*****************************************************************************/


// empcop.stat.c  --- "Utilities" ----------------------------------------------
//void K_array(int n, int p, double *J, double *K);
void K_array(int n, int p, const double J[], double *K);
//void L_array(int n, int p, double *K, double *L);
void L_array(int n, int p, const double K[], double *L);

double I_n(int n, int p, double *J, double *K, double *L);
double M_A_n(int n, int p, double *J, double *K, double *L, int A);

// different versions of "Array J" -----------------------
void J_m (int n, int p, const int b[], const double U[], const int R[], double *J);
void J_s (int n, int p,                const double U[],                double *J);

void J_sm(int n, int p, int q, const double U[], const int B[], double *J);
void J_u (int n, int p,        const double R[],                double *J);


// empcopm.test.c  --- Independence test among random vectors ------------------
void bootstrap(int *n, int *N, int *p, int *b, double *U, int *m,
	       double *MA0, double *I0, int *subset, char **subset_char,
	       int *pe);
void empirical_copula_test_rv(double *U, int *n, int *p, int *b, int *m, double *MA0,
			      double *I0, int *N, int *subset, double *MA, double *I,
			      double *pval, double *fisher, double *tippett, double *Ipval);

// empcops.test.c  --- Serial Indep. test Genest-Remillard (2004) --------------
void simulate_empirical_copula_serial(int *n, int *N, int *p, int *m,
				      double *TA0, double *G0, int *subset,
				      char **subset_char, double *fisher0,
				      double *tippett0, int *pe);
void empirical_copula_test_serial(double *U, int *n, int *p, int *m, double *TA0, double *G0,
				  int *N, int *subset, double *TA, double *G, double *pval,
				  double *fisher, double *tippett, double *globpval,
				  double *fisher0, double *tippett0);


// empcopsm.test.c  --- Multivar. Serial Indep. test ---------------------------
void bootstrap_serial(int *n, int *N, int *p, int *q, double *U, int *m,
		      double *MA0, double *I0, int *subset, char **subset_char,
		      int *pe);
void empirical_copula_test_rv_serial(double *U, int *n, int *p, int *q, int *m, double *MA0,
				     double *I0, int *N, int *subset, double *MA, double *I,
				     double *pval, double *fisher, double *tippett, double *Ipval);

// empcopu.test.c  --- Multivar.Indep. test Genest-Remillard (2004) ------------
void simulate_empirical_copula(int *n, int *N, int *p, int *m, double *TA0,
			       double *G0, int *subset, char **subset_char,
			       double *fisher0, double *tippett0, int *pe);
void empirical_copula_test(double *R, int *n, int *p, int *m, double *TA0, double *G0,
			   int *N, int *subset, double *TA, double *G, double *pval,
			   double *fisher, double *tippett, double *globpval,
			   double *fisher0, double *tippett0);




