import pandas as pd

from plaque_assay import consts, stats, utils

THRESHOLD = 50
WEAK_THRESHOLD = 60

# acceptable difference between expected and observed IC50 values in percentage difference
EPSILON = 15
# maximum mean_squared_error value expected for fitted model
MSE_PASS = 100


def perc_difference(x: float, y: float) -> float:
    """percentage difference between two numbers"""
    return (abs(x - y) / abs((x + y) / 2)) * 100


dilutions = [
    0.000025,
    0.000025,
    0.000250,
    0.000250,
    0.002500,
    0.002500,
    0.025000,
    0.025000,
]

perc_inf_weak = [
    116.263735,
    93.992355,
    113.685992,
    97.030688,
    122.177319,
    103.793316,
    52.342949,
    61.026772,
]

df_weak_inhibition = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_inf_weak}
)

perc_inf_no = [
    98.729667,
    100.000000,
    100.000000,
    97.147718,
    94.675382,
    100.000000,
    94.496768,
    100.000000,
]

df_no_inhibition = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_inf_no}
)


# percentage infected for something that should be able to fit a model
perc_good = [
    97.5,
    97.8,
    77.9,
    59.8,
    12.6,
    5.4,
    0.8,
    -0.06,
]


# https://github.com/FrancisCrickInstitute/plaque_assay/issues/60
perc_should_be_complete_or_fail = [
    99.90,
    26.05,
    4.099,
    37.047,
    51.86,
    68.83,
    49.118,
    39.655,
]


df_good_inhibition = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_good}
)


def test_calc_heuristics_dilutions():
    weak_inhibition_out = stats.calc_heuristics_dilutions(
        df_weak_inhibition, THRESHOLD, WEAK_THRESHOLD
    )
    assert utils.INT_TO_RESULT[weak_inhibition_out] == "weak inhibition"
    no_inhibition_out = stats.calc_heuristics_dilutions(
        df_no_inhibition, THRESHOLD, WEAK_THRESHOLD
    )
    assert utils.INT_TO_RESULT[no_inhibition_out] == "no inhibition"
    # check that good inhibition values don't use the heuristics
    # but instead get passed to the curve fitting algorithm
    good_inhib_out = stats.calc_heuristics_dilutions(
        df_good_inhibition, THRESHOLD, WEAK_THRESHOLD
    )
    assert good_inhib_out is None


def good_ic50_fit(perc_infected, expected_ic50) -> bool:
    df = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_infected})
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test",
        df=df,
        threshold=THRESHOLD,
        weak_threshold=WEAK_THRESHOLD,
    )
    good_fit = fit_method == "model fit"
    good_ic50 = perc_difference(result, expected_ic50) < EPSILON
    assert mean_squared_error is not None
    good_mse = mean_squared_error < MSE_PASS
    good_instance = isinstance(model_params, stats.ModelParams)
    return all([good_fit, good_ic50, good_mse, good_instance])


def test_calc_model_results():
    assert good_ic50_fit(perc_good, 2214)


def test_calc_64():
    perc_64 = [
        107.4,
        102.9,
        102.0,
        108.7,
        88.7,
        88.3,
        38.1,
        40.1,
    ]
    assert good_ic50_fit(perc_64, 64)


def test_calc_4795():
    perc_4795 = [
        96.1,
        85.7,
        49.6,
        39.7,
        3.2,
        6.5,
        0.1,
        0.0,
    ]
    assert good_ic50_fit(perc_4795, 4795)


def test_calc_408():
    perc_408 = [
        118.4,
        88.6,
        103.9,
        97.6,
        48.8,
        49.4,
        -0.06,
        -0.01,
    ]
    assert good_ic50_fit(perc_408, 408)


def test_calc_401():
    perc = [
        100.4,
        93.1,
        96.2,
        99.4,
        66.3,
        36.0,
        0.04,
        -0.21,
    ]
    assert good_ic50_fit(perc, 401)


def test_calc_1859():
    perc = [
        111.5,
        94.8,
        78.4,
        72.3,
        14.3,
        7.5,
        0.03,
        -0.006,
    ]
    assert good_ic50_fit(perc, 1859)


def test_calc_65():
    perc = [
        102.5,
        101.0,
        91.2,
        102.9,
        82.9,
        91.5,
        42.9,
        32.2,
    ]
    assert good_ic50_fit(perc, 65)


def test_should_be_complete_or_fail():
    """
    https://github.com/FrancisCrickInstitute/plaque_assay/issues/60
    """
    df = pd.DataFrame(
        {"Dilution": dilutions, "Percentage Infected": perc_should_be_complete_or_fail}
    )
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test",
        df=df,
        threshold=THRESHOLD,
        weak_threshold=WEAK_THRESHOLD,
    )
    result_str = utils.int_to_result(int(result))
    assert result_str in ("complete inhibition", "failed to fit model")
    # high MSE, should be flagged
    assert mean_squared_error is not None
    assert mean_squared_error > MSE_PASS
