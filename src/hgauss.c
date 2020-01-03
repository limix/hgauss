#include "hgauss.h"

#include <float.h>

#include "hcephes.h"
#define _USE_MATH_DEFINES
#include <math.h>

/*
 * double logcdf(double a)
 *
 * For a > -20, use the existing ndtr technique and take a log.
 * for a <= -20, we use the Taylor series approximation of erf to compute
 * the log CDF directly. The Taylor series consists of two parts which we will
 * name "left" and "right" accordingly.  The right part involves a summation
 * which we compute until the difference in terms falls below the
 * machine-specific EPSILON.
 *
 * \Phi(z) &=&
 *   \frac{e^{-z^2/2}}{-z\sqrt{2\pi}}  * [1 +  \sum_{n=1}^{N-1}  (-1)^n
 * \frac{(2n-1)!!}{(z^2)^n}]
 *   + O(z^{-2N+2})
 *   = [\mbox{LHS}] * [\mbox{RHS}] + \mbox{error}.
 *
 */
HGAUSS_API double hgauss_logcdf(double x)
{
    // we compute the left hand side of the approx (LHS) in one shot
    double log_LHS;
    // variable used to check for convergence
    double last_total = 0;
    // includes first term from the RHS summation
    double right_hand_side = 1;
    // numerator for RHS summand
    double numerator = 1;
    // use reciprocal for denominator to avoid division
    double denom_factor = 1;
    // the precomputed division we use to adjust the denominator
    double denom_cons = 1.0 / (x * x);
    long   sign = 1, i = 0;

    if (x > 6) {
        return -hcephes_ndtr(-x); /* log(1+x) \approx x */
    }
    if (x > -20) {
        return log(hcephes_ndtr(x));
    }
    log_LHS = -0.5 * x * x - log(-x) - 0.5 * log(2 * M_PI);

    while (fabs(last_total - right_hand_side) > DBL_EPSILON) {
        i += 1;
        last_total = right_hand_side;
        sign = -sign;
        denom_factor *= denom_cons;
        numerator *= 2 * i - 1;
        right_hand_side += sign * numerator * denom_factor;
    }
    return log_LHS + log(right_hand_side);
}

HGAUSS_API double hgauss_logpdf(double x)
{
    static const double logC = 0.91893853320467266954096885456237941980361938476;
    return -(x * x) / 2.0 - logC;
}

HGAUSS_API double hgauss_pdf(double x) { return exp(hgauss_logpdf(x)); }

HGAUSS_API double hgauss_cdf(double x) { return hcephes_ndtr(x); }

HGAUSS_API double hgauss_icdf(double x) { return hcephes_ndtri(x); }

HGAUSS_API double hgauss_logsf(double x) { return hgauss_logcdf(-x); }

HGAUSS_API double hgauss_sf(double x) { return hgauss_cdf(-x); }

HGAUSS_API double hgauss_isf(double x) { return -hgauss_icdf(x); }
