# Validation

`validate_paper.py` reproduces the correctly specified estimators from Sections
7.1 and 7.2 of the paper. It reports Monte Carlo bias and standard deviation
against the analytic truths:

- single-primary mean: `89 / 96`;
- two-primary product moment: `175 / 128`.

Run the full paper-sized experiment with:

```bash
cd python
uv run python ../validation/validate_paper.py --replications 1000 --n 2000
```

The committed report uses 200 replications as a fast independent verification;
CI runs unit and smoke simulations on every supported Python and R version.

## Committed result

With 200 replications of size 2,000, all correctly specified estimators had
small Monte Carlo bias:

| Design | IPW | RA | MR |
|---|---:|---:|---:|
| Single-primary bias | -0.0036 | 0.0009 | 0.0088 |
| Multiple-primary bias | 0.0016 | -0.0058 | 0.0021 |

The full machine-readable output, including empirical standard deviations and
the corresponding paper values, is in [`results.json`](results.json).

`parity_single.csv` is shipped with the R package and exercised by both test
suites. All three estimators must agree across R and Python to `1e-5` on that
shared fixture.
