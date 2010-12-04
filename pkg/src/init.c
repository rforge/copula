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

#include <R_ext/Rdynload.h>

#include "nacopula.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(sinc_c, 1),
    CALLDEF(A__c, 3),
    CALLDEF(polyn_eval, 2),

    CALLDEF(rstable_c, 2),
    CALLDEF(retstable_c, 4),

    CALLDEF(rLog_vec_c, 3),
    CALLDEF(rSibuya_vec_c, 2),

    CALLDEF(rF01Frank_vec_c, 5),
    CALLDEF(rF01Joe_vec_c, 3),

    {NULL, NULL, 0}
};

/**
 * register routines
 * @param dll pointer
 * @return none
 * @author Martin Maechler
 */
void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_nacopula(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
