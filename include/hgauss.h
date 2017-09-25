#ifndef HGAUSS_H
#define HGAUSS_H

#include "hcephes/hcephes.h"

double hgauss_logcdf(double x);

static inline double hgauss_pdf(double x) { return exp(hgauss_logpdf(x)); }

static inline double hgauss_cdf(double x) { return hcephes_ndtr(x); }

static inline double hgauss_icdf(double x) { return hcephes_ndtri(x); }

static inline double hgauss_logsf(double x) { return hgauss_logcdf(-x); }

static inline double hgauss_sf(double x) { return hgauss_cdf(-x); }

static inline double hgauss_isf(double x) { return -hgauss_icdf(x); }

static inline double hgauss_logpdf(double x) {
  static const double _norm_pdf_logC =
      0.9189385332046726695409688545623794198036193847656250;
  return -(x * x) / 2.0 - _norm_pdf_logC;
}

#endif /* end of include guard: HGAUSS_H */
