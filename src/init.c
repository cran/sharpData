
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "monomat.h"

static R_NativePrimitiveArgType kern_types[6] =
    {REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType afun_types[9] =
    {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, 
      INTSXP};

static R_FortranMethodDef fortranMethods[] = {
    {"kern", (DL_FUNC) &F77_SUB(kern), 6, kern_types},
    {"afun", (DL_FUNC) &F77_SUB(afun), 9, afun_types},
    {NULL, NULL, 0, NULL}
};
 
void attribute_visible R_init_sharpData(DllInfo *info)
{
    R_registerRoutines(info, NULL, NULL, fortranMethods, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
