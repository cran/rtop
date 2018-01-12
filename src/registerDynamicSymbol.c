#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .Fortran calls */
extern void F77_NAME(varioex)(void *, void *, void *, void *, void *);
extern void F77_NAME(vredhyp)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(vredind)(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"varioex", (DL_FUNC) &F77_NAME(varioex), 5},
  {"vredhyp", (DL_FUNC) &F77_NAME(vredhyp), 8},
  {"vredind", (DL_FUNC) &F77_NAME(vredind), 8},
  {NULL, NULL, 0}
};

void R_init_rtop(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
