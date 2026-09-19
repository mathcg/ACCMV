"""Reproduce the correctly specified estimators in paper Sections 7.1 and 7.2."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from accmv import (
    estimate_multiple,
    estimate_single,
    fit_ipw_regression,
    simulate_multiple_paper,
    simulate_regression_paper,
    simulate_single_paper,
)


def summarize(values: list[float], truth: float) -> dict[str, float]:
    array = np.asarray(values)
    return {
        "mean": float(array.mean()),
        "bias": float(array.mean() - truth),
        "sample_sd": float(array.std(ddof=1)),
    }


def run(replications: int, n: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    methods = ("ipw", "ra", "mr")
    single = {method: [] for method in methods}
    multiple = {method: [] for method in methods}
    regression_ipw = [[], []]
    regression_complete = [[], []]
    for _ in range(replications):
        x, y = simulate_single_paper(n, int(rng.integers(2**32 - 1)))
        for method in methods:
            single[method].append(estimate_single(x, y, method=method).estimate)
        x, y = simulate_multiple_paper(n, int(rng.integers(2**32 - 1)))
        for method in methods:
            multiple[method].append(
                estimate_multiple(x, y, method=method, target="product").estimate
            )
        x, y = simulate_regression_paper(n, int(rng.integers(2**32 - 1)))
        coefficients = fit_ipw_regression(x, y).coefficients
        complete = ~np.isnan(y).any(axis=1)
        complete_coefficients = np.linalg.lstsq(
            np.column_stack([np.ones(complete.sum()), y[complete, 0]]),
            y[complete, 1],
            rcond=None,
        )[0]
        for index in range(2):
            regression_ipw[index].append(coefficients[index])
            regression_complete[index].append(complete_coefficients[index])
    return {
        "configuration": {"replications": replications, "sample_size": n, "seed": seed},
        "single_primary": {
            "truth": 89 / 96,
            "paper_sample_sd": {"ipw": 0.216, "ra": 0.044, "mr": 0.105},
            "results": {key: summarize(value, 89 / 96) for key, value in single.items()},
        },
        "multiple_primary": {
            "truth": 175 / 128,
            "paper_sample_sd": {"ipw": 0.075, "ra": 0.065, "mr": 0.066},
            "results": {key: summarize(value, 175 / 128) for key, value in multiple.items()},
        },
        "marginal_regression": {
            "truth": {"intercept": -1.0, "slope": 0.5},
            "paper_sample_sd": {"intercept": 0.039, "slope": 0.046},
            "ipw": {
                "intercept": summarize(regression_ipw[0], -1.0),
                "slope": summarize(regression_ipw[1], 0.5),
            },
            "complete_case": {
                "intercept": summarize(regression_complete[0], -1.0),
                "slope": summarize(regression_complete[1], 0.5),
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=200)
    parser.add_argument("--n", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=220702289)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run(args.replications, args.n, args.seed)
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
