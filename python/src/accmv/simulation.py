"""Data-generating mechanisms from Section 7 of the ACCMV paper."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def _sigma(d: int) -> NDArray[np.float64]:
    return 0.5 * np.eye(d) + 0.5 * np.ones((d, d))


def simulate_single_paper(n: int = 2000, random_state: int | None = None):
    """Simulate Section 7.1; the true mean is ``89 / 96``."""

    if n < 1:
        raise ValueError("n must be positive")
    rng = np.random.default_rng(random_state)
    a = rng.integers(0, 2, n)
    r = rng.integers(0, 4, n)
    patterns = np.array([[False, False], [False, True], [True, False], [True, True]])
    x, y = np.full((n, 2), np.nan), np.full(n, np.nan)
    means = {1: np.array([1.0]), 2: np.array([1.0, -1.0]), 3: np.array([0.0, -1.0, -1.0])}
    for i in range(n):
        cols = np.flatnonzero(patterns[r[i]])
        if a[i]:
            draw = rng.multivariate_normal(means[len(cols) + 1], _sigma(len(cols) + 1))
            y[i], x[i, cols] = draw[0], draw[1:]
        elif len(cols):
            draw = rng.multivariate_normal(means[len(cols)], _sigma(len(cols)))
            x[i, cols] = draw
    return x, y


def simulate_multiple_paper(n: int = 2000, random_state: int | None = None):
    """Simulate Section 7.2; the true product moment is ``175 / 128``."""

    if n < 1:
        raise ValueError("n must be positive")
    rng = np.random.default_rng(random_state)
    a = rng.integers(0, 4, n)
    r = rng.integers(0, 4, n)
    patterns = np.array([[False, False], [False, True], [True, False], [True, True]])
    x, y = np.full((n, 2), np.nan), np.full((n, 2), np.nan)
    for i in range(n):
        xcols, ycols = np.flatnonzero(patterns[r[i]]), np.flatnonzero(patterns[a[i]])
        if a[i] == 3:
            d = 2 + len(xcols)
            draw = rng.multivariate_normal(np.ones(d), _sigma(d))
            y[i], x[i, xcols] = draw[:2], draw[2:]
        elif len(ycols) == 1:
            d = 1 + len(xcols)
            mean = np.full(d, 0.5 if d == 1 else 1.0)
            draw = rng.multivariate_normal(mean, _sigma(d))
            y[i, ycols], x[i, xcols] = draw[:1], draw[1:]
        elif len(xcols):
            d = len(xcols)
            mean = np.full(d, 0.5 if d == 1 else 1.0)
            x[i, xcols] = rng.multivariate_normal(mean, _sigma(d))
    return x, y
