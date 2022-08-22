#include <Rcpp.h>

#include <gsl/gsl_sf.h>

using namespace Rcpp;

#include "numerical_utils.hpp"

double log_add(double x, double y)
{
    // make x the max
    if (y > x) {
        double temp = x;
        x = y;
        y = temp;
    }
    // now x is bigger
    if (x == R_NegInf) {
        return x;
    }
    double negDiff = y - x;
    if (negDiff < -20) {
        return x;
    }
    return x + log(1.0 + exp(negDiff));
}

double log_beta_binomial_pdf(size_t x, size_t n, double alpha, double beta) {
    if (x > n)
    {
        return 0;
    }

    double ret = gsl_sf_lnchoose(n, x);
    ret += gsl_sf_lnbeta(x+alpha, n-x+beta);
    ret -= gsl_sf_lnbeta(alpha, beta);
    return ret;
}
