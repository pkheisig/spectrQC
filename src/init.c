#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _spectrQC_calc_wls_mat_cpp(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_spectrQC_calc_wls_mat_cpp", (DL_FUNC) &_spectrQC_calc_wls_mat_cpp, 4},
    {NULL, NULL, 0}
};

void R_init_spectrQC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
