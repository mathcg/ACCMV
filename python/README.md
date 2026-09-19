# accmv

Python implementation of inverse-probability-weighted (IPW), regression-adjustment
(RA), and multiply-robust (MR) inference under the available complete-case missing
value assumption.

```bash
python -m pip install accmv
```

```python
from accmv import estimate_single

fit = estimate_single(x, y, method="mr", n_boot=499, random_state=1)
print(fit.estimate, fit.conf_int)
```

For the paper's marginal linear model:

```python
from accmv import fit_ipw_regression

# Shared with R: response/predictors use one-based y-column indices.
fit = fit_ipw_regression(x, y, response=2, predictors=(1,))
print(fit.coefficients)
```

If you use this software, please cite Cheng, Chen, Smith, and Zhao,
“Handling Nonmonotone Missing Data with Available Complete-Case Missing Value
Assumption” ([paper](https://arxiv.org/abs/2207.02289)).
