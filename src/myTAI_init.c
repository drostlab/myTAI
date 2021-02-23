#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _myTAI_cpp_bootMatrix(SEXP, SEXP, SEXP);
extern SEXP _myTAI_cpp_geom_mean(SEXP);
extern SEXP _myTAI_cpp_harmonic_mean(SEXP);
extern SEXP _myTAI_cpp_omitMatrix(SEXP, SEXP);
extern SEXP _myTAI_cpp_pMatrix(SEXP, SEXP);
extern SEXP _myTAI_cpp_std_error(SEXP);
extern SEXP _myTAI_cpp_TAI(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
        {"_myTAI_cpp_bootMatrix",    (DL_FUNC) &_myTAI_cpp_bootMatrix,    3},
        {"_myTAI_cpp_geom_mean",     (DL_FUNC) &_myTAI_cpp_geom_mean,     1},
        {"_myTAI_cpp_harmonic_mean", (DL_FUNC) &_myTAI_cpp_harmonic_mean, 1},
        {"_myTAI_cpp_omitMatrix",    (DL_FUNC) &_myTAI_cpp_omitMatrix,    2},
        {"_myTAI_cpp_pMatrix",       (DL_FUNC) &_myTAI_cpp_pMatrix,       2},
        {"_myTAI_cpp_std_error",     (DL_FUNC) &_myTAI_cpp_std_error,     1},
        {"_myTAI_cpp_TAI",           (DL_FUNC) &_myTAI_cpp_TAI,           2},
        {NULL, NULL, 0}
};

void R_init_myTAI(DllInfo *dll)
{
        R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
        R_useDynamicSymbols(dll, FALSE);
}