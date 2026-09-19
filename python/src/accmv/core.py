"""Core ACCMV estimators.

The implementation follows equations (5), (7), (10), (14), (16), and (18) in
Cheng et al. (2022). Missingness is represented by ``numpy.nan``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

Method = Literal["ipw", "ra", "mr"]
Target = Literal["average", "product", "indicator"]


@dataclass(frozen=True)
class ACCMVData:
    """Validated data and stable missing-pattern encodings."""

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    x_patterns: NDArray[np.bool_]
    y_patterns: NDArray[np.bool_]
    r: NDArray[np.int64]
    a: NDArray[np.int64]

    @property
    def n(self) -> int:
        return self.y.shape[0]


@dataclass(frozen=True)
class ACCMVResult:
    """Point estimate with optional nonparametric-bootstrap inference."""

    estimate: float
    method: str
    target: str
    nobs: int
    bootstrap_values: NDArray[np.float64] | None = None
    std_error: float | None = None
    conf_int: tuple[float, float] | None = None


@dataclass(frozen=True)
class ACCMVRegressionResult:
    """Weighted linear-regression estimate for a marginal primary-data model."""

    coefficients: NDArray[np.float64]
    response: int
    predictors: tuple[int, ...]
    nobs: int
    weights: NDArray[np.float64]
    bootstrap_values: NDArray[np.float64] | None = None
    std_error: NDArray[np.float64] | None = None
    conf_int: NDArray[np.float64] | None = None


def _matrix(value: ArrayLike, name: str) -> NDArray[np.float64]:
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 1:
        arr = arr[:, None]
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a one- or two-dimensional numeric array")
    return arr


def _stable_patterns(mask: NDArray[np.bool_]) -> tuple[NDArray[np.bool_], NDArray[np.int64]]:
    lookup: dict[tuple[bool, ...], int] = {}
    patterns: list[tuple[bool, ...]] = []
    indices = np.empty(mask.shape[0], dtype=np.int64)
    for i, row in enumerate(mask):
        key = tuple(bool(v) for v in row)
        if key not in lookup:
            lookup[key] = len(patterns)
            patterns.append(key)
        indices[i] = lookup[key]
    return np.asarray(patterns, dtype=bool), indices


def prepare_data(x: ArrayLike, y: ArrayLike) -> ACCMVData:
    """Validate arrays and encode their observed-value patterns."""

    x_arr, y_arr = _matrix(x, "x"), _matrix(y, "y")
    if x_arr.shape[0] != y_arr.shape[0]:
        raise ValueError("x and y must contain the same number of rows")
    if x_arr.shape[0] < 2:
        raise ValueError("at least two observations are required")
    if np.isinf(x_arr).any() or np.isinf(y_arr).any():
        raise ValueError("x and y may contain NaN for missingness, but not infinity")
    xp, r = _stable_patterns(~np.isnan(x_arr))
    yp, a = _stable_patterns(~np.isnan(y_arr))
    return ACCMVData(x_arr, y_arr, xp, yp, r, a)


def _design(*parts: NDArray[np.float64]) -> NDArray[np.float64]:
    usable = [np.asarray(p, dtype=float).reshape(len(p), -1) for p in parts if p.size]
    n = len(parts[0]) if parts else 0
    return np.column_stack([np.ones(n), *usable])


def _quadratic_design(values: NDArray[np.float64]) -> NDArray[np.float64]:
    """Intercept, linear terms, and unique second-order products."""

    values = np.asarray(values, dtype=float).reshape(len(values), -1)
    products = [
        values[:, left] * values[:, right]
        for left in range(values.shape[1])
        for right in range(left, values.shape[1])
    ]
    return np.column_stack([np.ones(len(values)), values, *products])


def _linear_fit(x: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.linalg.lstsq(x, y, rcond=None)[0]


def _logistic_fit(x: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
    if np.unique(y).size != 2:
        raise ValueError("each fitted odds model needs observations in both comparison groups")
    beta = np.zeros(x.shape[1])
    ridge = 1e-10 * np.eye(x.shape[1])
    for _ in range(100):
        eta = np.clip(x @ beta, -30.0, 30.0)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = np.maximum(p * (1.0 - p), 1e-10)
        hessian = x.T @ (w[:, None] * x) + ridge
        step = np.linalg.solve(hessian, x.T @ (y - p))
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return beta


def _outcome_fit(
    x: NDArray[np.float64], y: NDArray[np.float64], binary: bool
) -> NDArray[np.float64]:
    return _logistic_fit(x, y) if binary else _linear_fit(x, y)


def _predict(x: NDArray[np.float64], beta: NDArray[np.float64], binary: bool) -> NDArray[np.float64]:
    eta = x @ beta
    return 1.0 / (1.0 + np.exp(-np.clip(eta, -30.0, 30.0))) if binary else eta


def _available(patterns: NDArray[np.bool_], observed: NDArray[np.int64]) -> NDArray[np.int64]:
    if observed.size == 0:
        return np.arange(patterns.shape[0])
    return np.flatnonzero(patterns[:, observed].all(axis=1))


def _single_arrays(data: ACCMVData) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
    if data.y.shape[1] != 1:
        raise ValueError("single-primary estimators require exactly one y column")
    y = data.y[:, 0]
    return y, ~np.isnan(y)


def single_regression_adjustment(
    data: ACCMVData, transform: Callable[[NDArray[np.float64]], NDArray[np.float64]] | None = None,
    *, binary: bool = False,
) -> float:
    """Regression-adjustment estimate for one primary variable."""

    y, observed_y = _single_arrays(data)
    fun = transform or (lambda z: z)
    outcome = np.asarray(fun(y), dtype=float)
    total = float(np.nansum(outcome[observed_y]))
    for pattern_id, pattern in enumerate(data.x_patterns):
        target = (~observed_y) & (data.r == pattern_id)
        if not target.any():
            continue
        cols = np.flatnonzero(pattern)
        train = observed_y & np.isin(data.r, _available(data.x_patterns, cols))
        design_train = _design(data.x[train][:, cols])
        design_target = _design(data.x[target][:, cols])
        beta = _outcome_fit(design_train, outcome[train], binary)
        total += float(_predict(design_target, beta, binary).sum())
    return total / data.n


def single_ipw(
    data: ACCMVData, transform: Callable[[NDArray[np.float64]], NDArray[np.float64]] | None = None,
    *, delta: float | None = None,
) -> float:
    """IPW estimate; ``delta`` enables exponential-tilt sensitivity analysis."""

    y, observed_y = _single_arrays(data)
    fun = transform or (lambda z: z)
    outcome = np.asarray(fun(y), dtype=float)
    contributions = np.zeros(data.n)
    weights = np.zeros(data.n)
    contributions[observed_y] = outcome[observed_y]
    weights[observed_y] = 1.0
    for pattern_id, pattern in enumerate(data.x_patterns):
        target = (~observed_y) & (data.r == pattern_id)
        if not target.any():
            continue
        cols = np.flatnonzero(pattern)
        donors = observed_y & np.isin(data.r, _available(data.x_patterns, cols))
        selected = donors | target
        design = _design(data.x[selected][:, cols])
        response = target[selected].astype(float)
        beta = _logistic_fit(design, response)
        odds = np.exp(np.clip(_design(data.x[donors][:, cols]) @ beta, -30.0, 30.0))
        if delta is not None:
            odds *= np.exp(np.clip(delta * y[donors], -30.0, 30.0))
        contributions[donors] += odds * outcome[donors]
        weights[donors] += odds
    denominator = weights.sum() if delta is not None else data.n
    return float(contributions.sum() / denominator)


def single_multiply_robust(
    data: ACCMVData, transform: Callable[[NDArray[np.float64]], NDArray[np.float64]] | None = None,
    *, binary: bool = False,
) -> float:
    """Multiply-robust estimate for one primary variable."""

    y, observed_y = _single_arrays(data)
    fun = transform or (lambda z: z)
    outcome = np.asarray(fun(y), dtype=float)
    total = float(outcome[observed_y].sum())
    for pattern_id, pattern in enumerate(data.x_patterns):
        target = (~observed_y) & (data.r == pattern_id)
        if not target.any():
            continue
        cols = np.flatnonzero(pattern)
        donors = observed_y & np.isin(data.r, _available(data.x_patterns, cols))
        selected = donors | target
        beta_o = _logistic_fit(_design(data.x[selected][:, cols]), target[selected].astype(float))
        odds = np.exp(np.clip(_design(data.x[donors][:, cols]) @ beta_o, -30.0, 30.0))
        beta_m = _outcome_fit(_design(data.x[donors][:, cols]), outcome[donors], binary)
        m0 = _predict(_design(data.x[donors][:, cols]), beta_m, binary)
        m1 = _predict(_design(data.x[target][:, cols]), beta_m, binary)
        total += float(((outcome[donors] - m0) * odds).sum() + m1.sum())
    return total / data.n


def _complete_y_pattern(data: ACCMVData) -> int:
    found = np.flatnonzero(data.y_patterns.all(axis=1))
    if found.size != 1:
        raise ValueError("the data must contain at least one complete primary-variable case")
    return int(found[0])


def _target_values(y: NDArray[np.float64], target: Target, threshold: float | None) -> NDArray[np.float64]:
    if y.shape[1] != 2:
        raise ValueError("multiple-primary targets currently require exactly two y columns")
    if target == "average":
        return 0.5 * (y[:, 0] + y[:, 1])
    if target == "product":
        return y[:, 0] * y[:, 1]
    if target == "indicator":
        if threshold is None:
            raise ValueError("threshold is required for the indicator target")
        return ((y[:, 0] <= threshold) & (y[:, 1] <= threshold)).astype(float)
    raise ValueError("target must be 'average', 'product', or 'indicator'")


def _conditional_target(
    data: ACCMVData,
    donors: NDArray[np.bool_],
    recipients: NDArray[np.bool_],
    xcols: NDArray[np.int64],
    ycols: NDArray[np.int64],
    target: Target,
    threshold: float | None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    missing = np.setdiff1d(np.arange(2), ycols)
    cov_train = _design(data.x[donors][:, xcols], data.y[donors][:, ycols])
    cov_new = _design(data.x[recipients][:, xcols], data.y[recipients][:, ycols])
    if ycols.size == 0:
        response = _target_values(data.y[donors], target, threshold)
        binary = target == "indicator"
        if target == "product":
            cov_train = _quadratic_design(data.x[donors][:, xcols])
            cov_new = _quadratic_design(data.x[recipients][:, xcols])
        beta = _outcome_fit(cov_train, response, binary)
        return _predict(cov_train, beta, binary), _predict(cov_new, beta, binary)
    if ycols.size != 1 or missing.size != 1:
        raise ValueError("only two primary variables are currently supported")
    if target == "indicator":
        response = (data.y[donors, missing[0]] <= float(threshold)).astype(float)
        beta = _outcome_fit(cov_train, response, True)
        p_train, p_new = _predict(cov_train, beta, True), _predict(cov_new, beta, True)
        obs_train = (data.y[donors, ycols[0]] <= float(threshold)).astype(float)
        obs_new = (data.y[recipients, ycols[0]] <= float(threshold)).astype(float)
        return obs_train * p_train, obs_new * p_new
    response = data.y[donors, missing[0]]
    beta = _linear_fit(cov_train, response)
    pred_train, pred_new = cov_train @ beta, cov_new @ beta
    obs_train = data.y[donors, ycols[0]]
    obs_new = data.y[recipients, ycols[0]]
    if target == "average":
        return 0.5 * (obs_train + pred_train), 0.5 * (obs_new + pred_new)
    return obs_train * pred_train, obs_new * pred_new


def _multiple_components(data: ACCMVData):
    complete_id = _complete_y_pattern(data)
    complete = data.a == complete_id
    for aid, ypattern in enumerate(data.y_patterns):
        if aid == complete_id:
            continue
        ycols = np.flatnonzero(ypattern)
        for rid, xpattern in enumerate(data.x_patterns):
            recipients = (data.a == aid) & (data.r == rid)
            if not recipients.any():
                continue
            xcols = np.flatnonzero(xpattern)
            donors = complete & np.isin(data.r, _available(data.x_patterns, xcols))
            if not donors.any():
                raise ValueError("a missingness pattern has no available complete-case donors")
            yield donors, recipients, xcols, ycols


def multiple_ipw(
    data: ACCMVData, *, target: Target = "average", threshold: float | None = None,
    delta: float | None = None,
) -> float:
    """IPW estimate for two primary variables."""

    complete = data.a == _complete_y_pattern(data)
    outcome = _target_values(data.y, target, threshold)
    contributions, weights = np.zeros(data.n), np.zeros(data.n)
    contributions[complete], weights[complete] = outcome[complete], 1.0
    for donors, recipients, xcols, ycols in _multiple_components(data):
        selected = donors | recipients
        beta = _logistic_fit(
            _design(data.x[selected][:, xcols], data.y[selected][:, ycols]),
            recipients[selected].astype(float),
        )
        odds = np.exp(np.clip(_design(data.x[donors][:, xcols], data.y[donors][:, ycols]) @ beta, -30, 30))
        if delta is not None:
            missing = np.setdiff1d(np.arange(2), ycols)
            tilt = data.y[donors][:, missing].mean(axis=1)
            odds *= np.exp(np.clip(delta * tilt, -30, 30))
        contributions[donors] += odds * outcome[donors]
        weights[donors] += odds
    denominator = weights.sum() if delta is not None else data.n
    return float(contributions.sum() / denominator)


def multiple_regression_adjustment(
    data: ACCMVData, *, target: Target = "average", threshold: float | None = None,
) -> float:
    """Regression-adjustment estimate for two primary variables."""

    complete = data.a == _complete_y_pattern(data)
    outcome = _target_values(data.y, target, threshold)
    total = float(outcome[complete].sum())
    for donors, recipients, xcols, ycols in _multiple_components(data):
        _, m1 = _conditional_target(data, donors, recipients, xcols, ycols, target, threshold)
        total += float(m1.sum())
    return total / data.n


def multiple_multiply_robust(
    data: ACCMVData, *, target: Target = "average", threshold: float | None = None,
) -> float:
    """Multiply-robust estimate for two primary variables."""

    complete = data.a == _complete_y_pattern(data)
    outcome = _target_values(data.y, target, threshold)
    total = float(outcome[complete].sum())
    for donors, recipients, xcols, ycols in _multiple_components(data):
        selected = donors | recipients
        beta = _logistic_fit(
            _design(data.x[selected][:, xcols], data.y[selected][:, ycols]),
            recipients[selected].astype(float),
        )
        odds = np.exp(np.clip(_design(data.x[donors][:, xcols], data.y[donors][:, ycols]) @ beta, -30, 30))
        m0, m1 = _conditional_target(data, donors, recipients, xcols, ycols, target, threshold)
        total += float(((outcome[donors] - m0) * odds).sum() + m1.sum())
    return total / data.n


def ipw_regression_weights(data: ACCMVData) -> NDArray[np.float64]:
    """Return ACCMV IPW weights for complete primary-variable cases."""

    complete = data.a == _complete_y_pattern(data)
    weights = np.zeros(data.n)
    weights[complete] = 1.0
    for donors, recipients, xcols, ycols in _multiple_components(data):
        selected = donors | recipients
        beta = _logistic_fit(
            _design(data.x[selected][:, xcols], data.y[selected][:, ycols]),
            recipients[selected].astype(float),
        )
        weights[donors] += np.exp(
            np.clip(_design(data.x[donors][:, xcols], data.y[donors][:, ycols]) @ beta, -30, 30)
        )
    return weights


def _regression_columns(
    data: ACCMVData, response: int, predictors: tuple[int, ...]
) -> tuple[int, NDArray[np.int64]]:
    """Validate shared one-based primary-variable column indices."""

    response_index = int(response) - 1
    predictor_indices = np.asarray(predictors, dtype=int) - 1
    if response_index < 0 or response_index >= data.y.shape[1]:
        raise ValueError("response must be a valid one-based y column index")
    if predictor_indices.size == 0:
        raise ValueError("predictors must contain at least one one-based y column index")
    if np.any(predictor_indices < 0) or np.any(predictor_indices >= data.y.shape[1]):
        raise ValueError("predictors contains an invalid one-based y column index")
    if response_index in predictor_indices:
        raise ValueError("response may not also be a predictor")
    if np.unique(predictor_indices).size != predictor_indices.size:
        raise ValueError("predictor indices must be unique")
    return response_index, predictor_indices


def _weighted_regression(
    data: ACCMVData, response: int, predictors: tuple[int, ...]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    response_index, predictor_indices = _regression_columns(data, response, predictors)
    weights = ipw_regression_weights(data)
    complete = ~np.isnan(data.y).any(axis=1)
    design = _design(data.y[complete][:, predictor_indices])
    root_weight = np.sqrt(weights[complete])
    coefficients = np.linalg.lstsq(
        design * root_weight[:, None],
        data.y[complete, response_index] * root_weight,
        rcond=None,
    )[0]
    return coefficients, weights


def fit_ipw_regression(
    x: ArrayLike,
    y: ArrayLike,
    *,
    response: int = 2,
    predictors: tuple[int, ...] = (1,),
    n_boot: int = 0,
    level: float = 0.95,
    random_state: int | None = None,
) -> ACCMVRegressionResult:
    """Fit an ACCMV-weighted linear model among complete primary cases.

    ``response`` and ``predictors`` are one-based column indices in both the R
    and Python APIs. The defaults fit the Section 7.3 model ``y2 ~ y1``.
    """

    data = prepare_data(x, y)
    predictors = tuple(int(value) for value in predictors)
    coefficients, weights = _weighted_regression(data, response, predictors)
    result = ACCMVRegressionResult(coefficients, response, predictors, data.n, weights)
    if n_boot == 0:
        return result
    if n_boot < 1:
        raise ValueError("n_boot must be nonnegative")
    if not 0 < level < 1:
        raise ValueError("level must lie strictly between zero and one")
    rng = np.random.default_rng(random_state)
    values: list[NDArray[np.float64]] = []
    attempts = 0
    while len(values) < n_boot and attempts < max(10 * n_boot, 100):
        attempts += 1
        try:
            sample = _resample(data, rng.integers(0, data.n, data.n))
            values.append(_weighted_regression(sample, response, predictors)[0])
        except (ValueError, np.linalg.LinAlgError):
            continue
    if len(values) != n_boot:
        raise RuntimeError("too many bootstrap samples lacked estimable pattern comparisons")
    bootstrap = np.asarray(values)
    alpha = 1.0 - level
    return replace(
        result,
        bootstrap_values=bootstrap,
        std_error=bootstrap.std(axis=0, ddof=1),
        conf_int=np.quantile(bootstrap, [alpha / 2, 1 - alpha / 2], axis=0).T,
    )


def _resample(data: ACCMVData, indices: NDArray[np.int64]) -> ACCMVData:
    return prepare_data(data.x[indices], data.y[indices])


def bootstrap_single(
    data: ACCMVData, *, method: Method = "mr", n_boot: int = 999,
    random_state: int | None = None, transform=None, binary: bool = False,
    delta: float | None = None,
) -> NDArray[np.float64]:
    """Nonparametric bootstrap estimates for a single-primary estimator."""

    if n_boot < 1:
        raise ValueError("n_boot must be positive")
    estimator = {"ipw": single_ipw, "ra": single_regression_adjustment, "mr": single_multiply_robust}.get(method)
    if estimator is None:
        raise ValueError("method must be 'ipw', 'ra', or 'mr'")
    rng = np.random.default_rng(random_state)
    values = []
    attempts = 0
    while len(values) < n_boot and attempts < max(10 * n_boot, 100):
        attempts += 1
        sample = _resample(data, rng.integers(0, data.n, data.n))
        try:
            kwargs = {"transform": transform}
            if method == "ipw":
                kwargs["delta"] = delta
            else:
                kwargs["binary"] = binary
            values.append(estimator(sample, **kwargs))
        except (ValueError, np.linalg.LinAlgError):
            continue
    if len(values) != n_boot:
        raise RuntimeError("too many bootstrap samples lacked estimable pattern comparisons")
    return np.asarray(values)


def bootstrap_multiple(
    data: ACCMVData, *, method: Method = "mr", target: Target = "average",
    threshold: float | None = None, n_boot: int = 999, random_state: int | None = None,
    delta: float | None = None,
) -> NDArray[np.float64]:
    """Nonparametric bootstrap estimates for a multiple-primary estimator."""

    if n_boot < 1:
        raise ValueError("n_boot must be positive")
    estimator = {"ipw": multiple_ipw, "ra": multiple_regression_adjustment, "mr": multiple_multiply_robust}.get(method)
    if estimator is None:
        raise ValueError("method must be 'ipw', 'ra', or 'mr'")
    rng = np.random.default_rng(random_state)
    values = []
    attempts = 0
    while len(values) < n_boot and attempts < max(10 * n_boot, 100):
        attempts += 1
        try:
            sample = _resample(data, rng.integers(0, data.n, data.n))
            kwargs = {"target": target, "threshold": threshold}
            if method == "ipw":
                kwargs["delta"] = delta
            values.append(estimator(sample, **kwargs))
        except (ValueError, np.linalg.LinAlgError):
            continue
    if len(values) != n_boot:
        raise RuntimeError("too many bootstrap samples lacked estimable pattern comparisons")
    return np.asarray(values)


def _with_bootstrap(result: ACCMVResult, values: NDArray[np.float64], level: float) -> ACCMVResult:
    if not 0 < level < 1:
        raise ValueError("level must lie strictly between zero and one")
    alpha = 1.0 - level
    ci = tuple(float(v) for v in np.quantile(values, [alpha / 2, 1 - alpha / 2]))
    return replace(result, bootstrap_values=values, std_error=float(values.std(ddof=1)), conf_int=ci)


def estimate_single(
    x: ArrayLike, y: ArrayLike, *, method: Method = "mr", transform=None,
    binary: bool = False, delta: float | None = None, n_boot: int = 0,
    level: float = 0.95, random_state: int | None = None,
) -> ACCMVResult:
    """Fit a single-primary ACCMV estimator with optional bootstrap inference."""

    data = prepare_data(x, y)
    if method == "ipw":
        value = single_ipw(data, transform, delta=delta)
    elif method == "ra":
        if delta is not None:
            raise ValueError("delta is only available for IPW sensitivity analysis")
        value = single_regression_adjustment(data, transform, binary=binary)
    elif method == "mr":
        if delta is not None:
            raise ValueError("delta is only available for IPW sensitivity analysis")
        value = single_multiply_robust(data, transform, binary=binary)
    else:
        raise ValueError("method must be 'ipw', 'ra', or 'mr'")
    result = ACCMVResult(value, method, "expectation", data.n)
    if n_boot:
        values = bootstrap_single(
            data, method=method, n_boot=n_boot, random_state=random_state,
            transform=transform, binary=binary, delta=delta,
        )
        result = _with_bootstrap(result, values, level)
    return result


def estimate_multiple(
    x: ArrayLike, y: ArrayLike, *, method: Method = "mr", target: Target = "average",
    threshold: float | None = None, delta: float | None = None, n_boot: int = 0,
    level: float = 0.95, random_state: int | None = None,
) -> ACCMVResult:
    """Fit a two-primary ACCMV estimator with optional bootstrap inference."""

    data = prepare_data(x, y)
    if method == "ipw":
        value = multiple_ipw(data, target=target, threshold=threshold, delta=delta)
    elif method == "ra":
        if delta is not None:
            raise ValueError("delta is only available for IPW sensitivity analysis")
        value = multiple_regression_adjustment(data, target=target, threshold=threshold)
    elif method == "mr":
        if delta is not None:
            raise ValueError("delta is only available for IPW sensitivity analysis")
        value = multiple_multiply_robust(data, target=target, threshold=threshold)
    else:
        raise ValueError("method must be 'ipw', 'ra', or 'mr'")
    result = ACCMVResult(value, method, target, data.n)
    if n_boot:
        values = bootstrap_multiple(
            data, method=method, target=target, threshold=threshold, n_boot=n_boot,
            random_state=random_state, delta=delta,
        )
        result = _with_bootstrap(result, values, level)
    return result
