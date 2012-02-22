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
#################################################################################
*/

#ifndef COPULA_GOF_H
#define COPULA_GOF_H

#include <R.h>

void cramer_vonMises(int *n, int *p, double *U, double *Ctheta,
		     double *stat);
void cramer_vonMises_2(int *p, double *U, int *n, double *V, int *m,
		       double *Ctheta, double *stat);

void multiplier(int *p, double *u0, int *m, double *u, int *n,
		double *influ, int *N, double *s0);

void cramer_vonMises_Pickands(int *n, int *m, double *S,
			      double *T, double *Atheta,
			      double *stat);

void cramer_vonMises_CFG(int *n, int *m, double *S,
			 double *T, double *Atheta,
			 double *stat);

void cramer_vonMises_Afun(int *n, int *m, double *S,
			  double *T, double *Atheta,
			  double *stat, int *CFG);

#endif
