// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// ReadParentVector
String ReadParentVector(String file_path, size_t mutation_count);
RcppExport SEXP _PhylExR_ReadParentVector(SEXP file_pathSEXP, SEXP mutation_countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type file_path(file_pathSEXP);
    Rcpp::traits::input_parameter< size_t >::type mutation_count(mutation_countSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadParentVector(file_path, mutation_count));
    return rcpp_result_gen;
END_RCPP
}
// GetChains
NumericVector GetChains(NumericVector parent_vector);
RcppExport SEXP _PhylExR_GetChains(SEXP parent_vectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parent_vector(parent_vectorSEXP);
    rcpp_result_gen = Rcpp::wrap(GetChains(parent_vector));
    return rcpp_result_gen;
END_RCPP
}
// get_parent_name
std::string get_parent_name(std::string node);
RcppExport SEXP _PhylExR_get_parent_name(SEXP nodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type node(nodeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_parent_name(node));
    return rcpp_result_gen;
END_RCPP
}
// GetConfigMatrix
NumericMatrix GetConfigMatrix(DataFrame datum2node, std::vector<std::string> ordered_nodes);
RcppExport SEXP _PhylExR_GetConfigMatrix(SEXP datum2nodeSEXP, SEXP ordered_nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type datum2node(datum2nodeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type ordered_nodes(ordered_nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetConfigMatrix(datum2node, ordered_nodes));
    return rcpp_result_gen;
END_RCPP
}
// LogAdd
double LogAdd(double x, double y);
RcppExport SEXP _PhylExR_LogAdd(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(LogAdd(x, y));
    return rcpp_result_gen;
END_RCPP
}
// LogSumExp
double LogSumExp(NumericVector x);
RcppExport SEXP _PhylExR_LogSumExp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(LogSumExp(x));
    return rcpp_result_gen;
END_RCPP
}
// ScLikelihoodWithDropout
double ScLikelihoodWithDropout(bool has_snv, size_t var_reads, size_t total_reads, double dropout_alpha, double dropout_beta, double bursty_alpha, double bursty_beta, double biallelic_alpha, double biallelic_beta, double seq_err, double dropout_mixing_proportion, double bursty_mixing_proportion, double biallelic_mixing_proportion);
RcppExport SEXP _PhylExR_ScLikelihoodWithDropout(SEXP has_snvSEXP, SEXP var_readsSEXP, SEXP total_readsSEXP, SEXP dropout_alphaSEXP, SEXP dropout_betaSEXP, SEXP bursty_alphaSEXP, SEXP bursty_betaSEXP, SEXP biallelic_alphaSEXP, SEXP biallelic_betaSEXP, SEXP seq_errSEXP, SEXP dropout_mixing_proportionSEXP, SEXP bursty_mixing_proportionSEXP, SEXP biallelic_mixing_proportionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type has_snv(has_snvSEXP);
    Rcpp::traits::input_parameter< size_t >::type var_reads(var_readsSEXP);
    Rcpp::traits::input_parameter< size_t >::type total_reads(total_readsSEXP);
    Rcpp::traits::input_parameter< double >::type dropout_alpha(dropout_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type dropout_beta(dropout_betaSEXP);
    Rcpp::traits::input_parameter< double >::type bursty_alpha(bursty_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type bursty_beta(bursty_betaSEXP);
    Rcpp::traits::input_parameter< double >::type biallelic_alpha(biallelic_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type biallelic_beta(biallelic_betaSEXP);
    Rcpp::traits::input_parameter< double >::type seq_err(seq_errSEXP);
    Rcpp::traits::input_parameter< double >::type dropout_mixing_proportion(dropout_mixing_proportionSEXP);
    Rcpp::traits::input_parameter< double >::type bursty_mixing_proportion(bursty_mixing_proportionSEXP);
    Rcpp::traits::input_parameter< double >::type biallelic_mixing_proportion(biallelic_mixing_proportionSEXP);
    rcpp_result_gen = Rcpp::wrap(ScLikelihoodWithDropout(has_snv, var_reads, total_reads, dropout_alpha, dropout_beta, bursty_alpha, bursty_beta, biallelic_alpha, biallelic_beta, seq_err, dropout_mixing_proportion, bursty_mixing_proportion, biallelic_mixing_proportion));
    return rcpp_result_gen;
END_RCPP
}
// ScLikelihood
double ScLikelihood(bool has_snv, size_t var_reads, size_t total_reads, double biallelic_alpha, double biallelic_beta, double bursty_alpha, double bursty_beta, double seq_err, double bursty_mixture_prob);
RcppExport SEXP _PhylExR_ScLikelihood(SEXP has_snvSEXP, SEXP var_readsSEXP, SEXP total_readsSEXP, SEXP biallelic_alphaSEXP, SEXP biallelic_betaSEXP, SEXP bursty_alphaSEXP, SEXP bursty_betaSEXP, SEXP seq_errSEXP, SEXP bursty_mixture_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type has_snv(has_snvSEXP);
    Rcpp::traits::input_parameter< size_t >::type var_reads(var_readsSEXP);
    Rcpp::traits::input_parameter< size_t >::type total_reads(total_readsSEXP);
    Rcpp::traits::input_parameter< double >::type biallelic_alpha(biallelic_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type biallelic_beta(biallelic_betaSEXP);
    Rcpp::traits::input_parameter< double >::type bursty_alpha(bursty_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type bursty_beta(bursty_betaSEXP);
    Rcpp::traits::input_parameter< double >::type seq_err(seq_errSEXP);
    Rcpp::traits::input_parameter< double >::type bursty_mixture_prob(bursty_mixture_probSEXP);
    rcpp_result_gen = Rcpp::wrap(ScLikelihood(has_snv, var_reads, total_reads, biallelic_alpha, biallelic_beta, bursty_alpha, bursty_beta, seq_err, bursty_mixture_prob));
    return rcpp_result_gen;
END_RCPP
}
// IdentifyCellMutationStatus
NumericMatrix IdentifyCellMutationStatus(DataFrame datum2node, std::vector<std::string> ordered_nodes, std::vector<std::string> ordered_mutations, NumericMatrix var_counts, NumericMatrix total_counts, NumericMatrix dropout_hp, NumericMatrix bursty_hp, NumericMatrix biallelic_hp);
RcppExport SEXP _PhylExR_IdentifyCellMutationStatus(SEXP datum2nodeSEXP, SEXP ordered_nodesSEXP, SEXP ordered_mutationsSEXP, SEXP var_countsSEXP, SEXP total_countsSEXP, SEXP dropout_hpSEXP, SEXP bursty_hpSEXP, SEXP biallelic_hpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type datum2node(datum2nodeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type ordered_nodes(ordered_nodesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type ordered_mutations(ordered_mutationsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type var_counts(var_countsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type total_counts(total_countsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dropout_hp(dropout_hpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bursty_hp(bursty_hpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type biallelic_hp(biallelic_hpSEXP);
    rcpp_result_gen = Rcpp::wrap(IdentifyCellMutationStatus(datum2node, ordered_nodes, ordered_mutations, var_counts, total_counts, dropout_hp, bursty_hp, biallelic_hp));
    return rcpp_result_gen;
END_RCPP
}
// IdentifyCellMutationStatusBursty
NumericMatrix IdentifyCellMutationStatusBursty(DataFrame datum2node, std::vector<std::string> ordered_nodes, std::vector<std::string> ordered_mutations, NumericMatrix var_counts, NumericMatrix total_counts, NumericMatrix bursty_hp, NumericMatrix biallelic_hp);
RcppExport SEXP _PhylExR_IdentifyCellMutationStatusBursty(SEXP datum2nodeSEXP, SEXP ordered_nodesSEXP, SEXP ordered_mutationsSEXP, SEXP var_countsSEXP, SEXP total_countsSEXP, SEXP bursty_hpSEXP, SEXP biallelic_hpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type datum2node(datum2nodeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type ordered_nodes(ordered_nodesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type ordered_mutations(ordered_mutationsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type var_counts(var_countsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type total_counts(total_countsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bursty_hp(bursty_hpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type biallelic_hp(biallelic_hpSEXP);
    rcpp_result_gen = Rcpp::wrap(IdentifyCellMutationStatusBursty(datum2node, ordered_nodes, ordered_mutations, var_counts, total_counts, bursty_hp, biallelic_hp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PhylExR_ReadParentVector", (DL_FUNC) &_PhylExR_ReadParentVector, 2},
    {"_PhylExR_GetChains", (DL_FUNC) &_PhylExR_GetChains, 1},
    {"_PhylExR_get_parent_name", (DL_FUNC) &_PhylExR_get_parent_name, 1},
    {"_PhylExR_GetConfigMatrix", (DL_FUNC) &_PhylExR_GetConfigMatrix, 2},
    {"_PhylExR_LogAdd", (DL_FUNC) &_PhylExR_LogAdd, 2},
    {"_PhylExR_LogSumExp", (DL_FUNC) &_PhylExR_LogSumExp, 1},
    {"_PhylExR_ScLikelihoodWithDropout", (DL_FUNC) &_PhylExR_ScLikelihoodWithDropout, 13},
    {"_PhylExR_ScLikelihood", (DL_FUNC) &_PhylExR_ScLikelihood, 9},
    {"_PhylExR_IdentifyCellMutationStatus", (DL_FUNC) &_PhylExR_IdentifyCellMutationStatus, 8},
    {"_PhylExR_IdentifyCellMutationStatusBursty", (DL_FUNC) &_PhylExR_IdentifyCellMutationStatusBursty, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_PhylExR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
