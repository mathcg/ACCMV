# ACCMV

[![Python and R CI](https://github.com/mathcg/ACCMV/actions/workflows/ci.yml/badge.svg)](https://github.com/mathcg/ACCMV/actions/workflows/ci.yml)
[![GitHub release](https://img.shields.io/github/v/release/mathcg/ACCMV)](https://github.com/mathcg/ACCMV/releases/latest)
[![PyPI](https://img.shields.io/pypi/v/accmv)](https://pypi.org/project/accmv/)

`ACCMV` implements inference under the **available complete-case missing value
assumption** for nonmonotone missing data. The repository contains matched R
and Python packages with IPW, regression-adjustment, and multiply-robust
estimation; exponential-tilt sensitivity analysis; marginal-regression
weights; and nonparametric bootstrap confidence intervals.

The estimators support a single primary variable and two-primary-variable
average, product-moment, and joint-distribution targets. Missing values are
represented by `NA` in R and `numpy.nan` in Python.

## Installation

Version 0.1.1 has been submitted to CRAN and is awaiting CRAN review. Until it
appears in the CRAN package index, install the validated release from GitHub:

```r
# install.packages("remotes")
remotes::install_github("mathcg/ACCMV", subdir = "r")
```

```bash
python -m pip install accmv
```

## Quick start

```r
library(accmv)
fit <- estimate_accmv_single(x, y, method = "mr",
                             n_boot = 499, seed = 1)
fit
```

```python
from accmv import estimate_single

fit = estimate_single(x, y, method="mr", n_boot=499, random_state=1)
print(fit.estimate, fit.conf_int)
```

The paper's marginal linear model is available in both languages through
`fit_accmv_regression()` (R) and `fit_ipw_regression()` (Python). Both use
one-based primary-variable column indices so equivalent calls have the same
arguments.

See [`docs/METHODS.md`](docs/METHODS.md) for the API-to-paper mapping and
[`validation/README.md`](validation/README.md) for reproducibility checks. The
original research scripts and diabetes demonstration remain at the repository
root for provenance.

## Reference

If you use this software or the ACCMV method, please cite:

> Cheng, G., Chen, Y.-C., Smith, M. A., and Zhao, Y.-Q. (2022).
> “Handling Nonmonotone Missing Data with Available Complete-Case Missing
> Value Assumption.” [arXiv:2207.02289](https://arxiv.org/abs/2207.02289).

## License

MIT © Gang Cheng.
