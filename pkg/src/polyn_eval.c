#include "nacopula.h"

SEXP polyn_eval(SEXP coef, SEXP x)
{
 SEXP result;
 int n = LENGTH(x), i, j,
     m = LENGTH(coef);
 // deal with integer or numeric -- NULL cannot (yet?) be coerced
 if(isNull(x)) { result = allocVector(REALSXP, 0); return result; }
 if(!isNull(coef)) coef = coerceVector(coef, REALSXP);
 PROTECT(coef);
 PROTECT(x = coerceVector(x, REALSXP));
 PROTECT(result = Rf_duplicate(x));
 double *cf = REAL(coef), *xx = REAL(x), *res = REAL(result);
 for(i = 0; i < n; i++) {
   double r, xi = xx[i];
   if(m == 0) {
      r = 0.;
   } else {
      j = m-1;
      r = cf[j];
      for (j--; j >= 0; j--)
          r = cf[j] + r * xi;
   }
   res[i] = r;
 }
 UNPROTECT(3);
 return result;
}
