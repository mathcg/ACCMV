# Validation

`validate_paper.py` reproduces the correctly specified estimators from Sections
7.1, 7.2, and 7.3 of the paper. It reports Monte Carlo bias and standard deviation
against the analytic truths:

- single-primary mean: `89 / 96`;
- two-primary product moment: `175 / 128`.
- marginal-regression intercept and slope: `-1` and `0.5`.

Run the full paper-sized experiment with:

```bash
cd python
uv run python ../validation/validate_paper.py --replications 1000 --n 2000
```

The committed report uses the paper's 1,000 replications; CI runs unit and
smoke simulations on every supported Python and R version.

## Committed result

With 1,000 replications of size 2,000, all correctly specified estimators had
small Monte Carlo bias:

| Design | IPW | RA | MR |
|---|---:|---:|---:|
| Single-primary bias | -0.0144 | 0.0006 | 0.0014 |
| Multiple-primary bias | 0.0016 | 0.0005 | 0.0004 |

For Section 7.3, the ACCMV-weighted regression biases were `0.0008` for the
intercept and `-0.0004` for the slope. The complete-case biases were `-0.0572`
and `-0.0066`, closely reproducing Table 3's `-0.061` and `-0.008`.

The full machine-readable output, including empirical standard deviations and
the corresponding paper values, is in [`results.json`](results.json).

`parity_single.csv` is shipped with the R package and exercised by both test
suites. All three estimators must agree across R and Python to `1e-5` on that
shared fixture.
