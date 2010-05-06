#include <R_ext/Rdynload.h>

#include "nacopula.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(sinc_c, 1),
    CALLDEF(A__c, 3),

    CALLDEF(rstable_c, 2),
    CALLDEF(retstable_c, 4),

    CALLDEF(rLog_c, 2),
    CALLDEF(rFJoe_c, 2),
    CALLDEF(rFFrank_c, 3),

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
