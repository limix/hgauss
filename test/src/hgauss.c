#include "hgauss.h"
#include <math.h>

#define assert_close(x, y)                                                                    \
    if (fabs((x) - (y)) > 1e-7)                                                               \
        return 1;

int main()
{
    assert_close(hgauss_logpdf(1.3), -1.7639385332);
    assert_close(hgauss_logpdf(-0.23558029373548589), -0.94668757);
    assert_close(hgauss_logcdf(1.3), -0.101811802668);
    assert_close(hgauss_logcdf(-0.23558029373548589), -0.89923898);

    return 0;
}
