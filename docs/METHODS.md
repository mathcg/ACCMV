# Methods and implementation map

The package follows Cheng, Chen, Smith, and Zhao (2022). For a secondary
variable pattern `r`, available complete cases observe every component in `r`;
this is the partial order `R >= r` used in the paper.

| Paper method | Python | R |
|---|---|---|
| Single-primary RA, equation (5) | `single_regression_adjustment()` | `single_regression_adjustment()` |
| Single-primary IPW, equation (7) | `single_ipw()` | `single_ipw()` |
| Single-primary MR, equation (10) | `single_multiply_robust()` | `single_multiply_robust()` |
| Multiple-primary RA, equation (14) | `multiple_regression_adjustment()` | `multiple_ra_*()` |
| Multiple-primary IPW, equation (16) | `multiple_ipw()` | `multiple_ipw()` |
| Multiple-primary MR, equation (18) | `multiple_multiply_robust()` | `multiple_mr_*()` |
| Marginal-model weights, Section 5 | `ipw_regression_weights()` | `accmv_ipw_weights()` |
| Sensitivity analysis, Section 6 | IPW `delta=` | IPW `delta=` |

Log-odds models include an intercept and all observed primary and secondary
variables permitted by the pattern. Continuous regression models use least
squares and indicator targets use logistic regression.

The high-level entry points are `estimate_single()` and `estimate_multiple()`
in Python and `estimate_accmv_single()` and `estimate_accmv_multiple()` in R.
They return a point estimate and optional percentile-bootstrap confidence
interval.

The RA and MR helpers currently implement the paper's one-primary and
two-primary examples. Positivity remains a scientific assumption: sparse
patterns or separated odds models should be diagnosed and may require a
substantively justified model change.
