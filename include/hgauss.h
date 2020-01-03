#ifndef HGAUSS_H
#define HGAUSS_H

#ifdef _WIN32
#ifdef HGAUSS_EXPORTS
#define HGAUSS_API __declspec(dllexport)
#else
#define HGAUSS_API __declspec(dllimport)
#endif
#else
#define HGAUSS_API
#endif

/** Major liknorm version. */
#define HGAUSS_VERSION_MAJOR 0
/** Minor liknorm version. */
#define HGAUSS_VERSION_MINOR 1
/** Minor liknorm version. */
#define HGAUSS_VERSION_PATCH 0
/** Liknorm version. */
#define HGAUSS_VERSION "0.1.0"

#ifdef __cplusplus
extern "C"
{
#endif

    HGAUSS_API double hgauss_logcdf(double x);
    HGAUSS_API double hgauss_logpdf(double x);
    HGAUSS_API double hgauss_pdf(double x);
    HGAUSS_API double hgauss_cdf(double x);
    HGAUSS_API double hgauss_icdf(double x);
    HGAUSS_API double hgauss_logsf(double x);
    HGAUSS_API double hgauss_sf(double x);
    HGAUSS_API double hgauss_isf(double x);

#ifdef __cplusplus
}
#endif

#endif
