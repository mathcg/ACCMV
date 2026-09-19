# accmv

Python implementation of inverse-probability-weighted (IPW), regression-adjustment
(RA), and multiply-robust (MR) inference under the available complete-case missing
value assumption.

```python
from accmv import estimate_single

fit = estimate_single(x, y, method="mr", n_boot=499, random_state=1)
print(fit.estimate, fit.conf_int)
```

If you use this software, please cite Cheng, Chen, Smith, and Zhao,
“Handling Nonmonotone Missing Data with Available Complete-Case Missing Value
Assumption” ([paper](https://arxiv.org/abs/2207.02289)).
