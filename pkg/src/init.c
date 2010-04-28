#include <R_ext/Rdynload.h>

#include "nacopula.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(sinc, 1),
    CALLDEF(A_Zolotarev, 3),

    CALLDEF(rLog_c, 2),
    CALLDEF(rFJoe_c, 2),
    CALLDEF(rFFrank_c, 3),

    CALLDEF(rstable_b1, 3),
    CALLDEF(retstable_c, 2),

    {NULL, NULL, 0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_nacopula(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
