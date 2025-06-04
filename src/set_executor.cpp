#include "Rtatami.h"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
SEXP set_executor(SEXP ptr) {
    Rtatami::set_executor(ptr);
    return R_NilValue;
}
