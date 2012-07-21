/*
  Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan

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

/**
 * @file   empcop.h
 * @author Ivan Kojadinovic 
 * @date   Sat Jul 21 17:25:20 2012
 * 
 * @brief  Bivariate and multivariate versions of the empirical copula 
 *         and related estimators of the partial derivatives 
 * 
 */

#ifndef EMPCOP_H
#define EMPCOP_H

/// Bivariate versions; used by exchTest and evTestA
double bivCn(double *U, double *V, int n, double u, double v);
double der1bivCn(double *U, double *V, int n, double u, double v);
double der2bivCn(double *U, double *V, int n, double u, double v);

/// Multivariate versions; used by evTestC and by the multiplier gof tests
double multCn(double *U, int n, int p, double *V, int m, int k, double o);
double der_multCn(double *U, int n, int p, double *u, double *v, double denom);

#endif
