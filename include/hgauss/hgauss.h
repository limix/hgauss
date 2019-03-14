#ifndef HGAUSS_H
#define HGAUSS_H

#include "hcephes.h"

double hgauss_logcdf(double x);

static inline double hgauss_logpdf(double x) {
  static const double logC = 0.91893853320467266954096885456237941980361938476;
  return -(x * x) / 2.0 - logC;
}

static inline double hgauss_pdf(double x) { return exp(hgauss_logpdf(x)); }

static inline double hgauss_cdf(double x) { return hcephes_ndtr(x); }

static inline double hgauss_icdf(double x) { return hcephes_ndtri(x); }

static inline double hgauss_logsf(double x) { return hgauss_logcdf(-x); }

static inline double hgauss_sf(double x) { return hgauss_cdf(-x); }

static inline double hgauss_isf(double x) { return -hgauss_icdf(x); }

#endif /* end of include guard: HGAUSS_H */
