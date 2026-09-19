from pathlib import Path

import numpy as np
import pytest

from accmv import (
    estimate_multiple,
    estimate_single,
    ipw_regression_weights,
    prepare_data,
    simulate_multiple_paper,
    simulate_single_paper,
)


def test_prepare_data_uses_stable_patterns():
    x = np.array([[1.0, np.nan], [np.nan, 2.0], [3.0, np.nan]])
    y = np.array([1.0, np.nan, 2.0])
    data = prepare_data(x, y)
    assert data.x_patterns.tolist() == [[True, False], [False, True]]
    assert data.r.tolist() == [0, 1, 0]


@pytest.mark.parametrize("method", ["ra", "mr"])
def test_single_paper_estimate_is_near_truth(method):
    x, y = simulate_single_paper(20_000, random_state=47)
    fit = estimate_single(x, y, method=method)
    assert fit.estimate == pytest.approx(89 / 96, abs=0.08)


def test_single_ipw_is_finite():
    x, y = simulate_single_paper(10_000, random_state=2)
    assert np.isfinite(estimate_single(x, y, method="ipw").estimate)


@pytest.mark.parametrize("method", ["ipw", "ra", "mr"])
def test_multiple_paper_product_is_near_truth(method):
    x, y = simulate_multiple_paper(20_000, random_state=91)
    fit = estimate_multiple(x, y, method=method, target="product")
    assert fit.estimate == pytest.approx(175 / 128, abs=0.12)


def test_bootstrap_is_reproducible():
    x, y = simulate_single_paper(1000, random_state=9)
    one = estimate_single(x, y, method="ra", n_boot=19, random_state=3)
    two = estimate_single(x, y, method="ra", n_boot=19, random_state=3)
    np.testing.assert_array_equal(one.bootstrap_values, two.bootstrap_values)
    assert one.std_error > 0
    assert one.conf_int[0] < one.conf_int[1]


def test_invalid_method_is_rejected():
    x, y = simulate_single_paper(100, random_state=1)
    with pytest.raises(ValueError, match="method"):
        estimate_single(x, y, method="wrong")


def test_cross_language_parity_fixture():
    fixture = Path(__file__).parents[2] / "r" / "inst" / "extdata" / "parity_single.csv"
    values = np.genfromtxt(fixture, delimiter=",", names=True, missing_values="NA", filling_values=np.nan)
    x = np.column_stack([values["x1"], values["x2"]])
    expected = {
        "ipw": 0.8258232120833385,
        "ra": 0.9233429113247089,
        "mr": 0.7096728452585138,
    }
    for method, truth in expected.items():
        assert estimate_single(x, values["y"], method=method).estimate == pytest.approx(truth, abs=1e-8)


def test_multiple_targets_sensitivity_and_regression_weights():
    x, y = simulate_multiple_paper(3000, random_state=17)
    indicator = estimate_multiple(x, y, method="mr", target="indicator", threshold=1.0)
    sensitivity = estimate_multiple(x, y, method="ipw", target="average", delta=0.1)
    data = prepare_data(x, y)
    weights = ipw_regression_weights(data)
    complete = ~np.isnan(y).any(axis=1)
    assert 0 <= indicator.estimate <= 1
    assert np.isfinite(sensitivity.estimate)
    assert np.all(weights >= 0)
    assert np.all(weights[~complete] == 0)
    assert np.all(weights[complete] >= 1)
