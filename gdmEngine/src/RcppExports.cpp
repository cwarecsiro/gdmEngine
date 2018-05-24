// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// BigGridTransform
int BigGridTransform(int nrows, int ncols, int nlayers, Rcpp::StringVector InFilePaths, Rcpp::StringVector OutFilePaths);
RcppExport SEXP _gdmEngine_BigGridTransform(SEXP nrowsSEXP, SEXP ncolsSEXP, SEXP nlayersSEXP, SEXP InFilePathsSEXP, SEXP OutFilePathsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nlayers(nlayersSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type InFilePaths(InFilePathsSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type OutFilePaths(OutFilePathsSEXP);
    rcpp_result_gen = Rcpp::wrap(BigGridTransform(nrows, ncols, nlayers, InFilePaths, OutFilePaths));
    return rcpp_result_gen;
END_RCPP
}
// PairsDissim
NumericVector PairsDissim(IntegerMatrix site_spp, IntegerMatrix pair_rows, IntegerVector site_rich, int max_richness);
RcppExport SEXP _gdmEngine_PairsDissim(SEXP site_sppSEXP, SEXP pair_rowsSEXP, SEXP site_richSEXP, SEXP max_richnessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type site_spp(site_sppSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type pair_rows(pair_rowsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type site_rich(site_richSEXP);
    Rcpp::traits::input_parameter< int >::type max_richness(max_richnessSEXP);
    rcpp_result_gen = Rcpp::wrap(PairsDissim(site_spp, pair_rows, site_rich, max_richness));
    return rcpp_result_gen;
END_RCPP
}
// PairsDist
NumericVector PairsDist(NumericMatrix env_vars, NumericMatrix pair_rows);
RcppExport SEXP _gdmEngine_PairsDist(SEXP env_varsSEXP, SEXP pair_rowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type env_vars(env_varsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pair_rows(pair_rowsSEXP);
    rcpp_result_gen = Rcpp::wrap(PairsDist(env_vars, pair_rows));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gdmEngine_BigGridTransform", (DL_FUNC) &_gdmEngine_BigGridTransform, 5},
    {"_gdmEngine_PairsDissim", (DL_FUNC) &_gdmEngine_PairsDissim, 4},
    {"_gdmEngine_PairsDist", (DL_FUNC) &_gdmEngine_PairsDist, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_gdmEngine(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
