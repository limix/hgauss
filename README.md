# hgauss

[![Travis](https://travis-ci.com/limix/hgauss.svg?branch=master)](https://travis-ci.com/limix/hgauss)

Auxiliary functions for Normal distribution.

## Install

You can install it via [conda](https://conda.io)

```bash
conda install -c conda-forge hgauss
```

Alternatively, one can compile and install it.
From Linux, MacOS, or Windows (bash terminal) systems, enter

```bash
# DO_CMD=sudo
curl -fsSL https://git.io/JerYI | GITHUB_USER=limix GITHUB_PROJECT=hgauss bash
```

## Documentation

```c
double hgauss_logcdf(double x);
double hgauss_pdf(double x);
double hgauss_cdf(double x);
double hgauss_icdf(double x);
double hgauss_logsf(double x);
double hgauss_sf(double x);
double hgauss_isf(double x);
double hgauss_logpdf(double x);
```

## Authors

- [Danilo Horta](https://github.com/horta)

## License

This project is licensed under the MIT License - see the [license file](https://raw.githubusercontent.com/limix/hgauss/master/LICENSE.md) for details.
