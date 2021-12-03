from plaque_assay import stats
from plaque_assay import utils

import pandas as pd


THRESHOLD = 50
WEAK_THRESHOLD = 60

# acceptable difference between expected and observed IC50 values
EPSILON = 10
# maximum mean_squared_error value expected for fitted model
MSE_PASS = 100


dilutions = [
    0.000391,
    0.000391,
    0.001563,
    0.001563,
    0.006250,
    0.006250,
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
    100.556437,
    102.200186,
    80.246412,
    82.365569,
    60.787072,
    54.955933,
    12.517334,
    13.988952,
]

# percentage infected values for an IC50 of ~391
perc_391 = [
    104.163,
    91.075,
    78.954,
    77.688,
    8.487,
    8.092,
    3.657,
    -0.475,
]

# percentage infected values for an IC50 of ~184
perc_184 = [
    120.524,
    118.954,
    123.209,
    119.373,
    20.256,
    14.863,
    0.540,
    -0.412,
]

# percentage infected values for an IC50 of ~381
perc_381 = [
    94.035,
    84.759,
    76.207,
    76.039,
    9.250,
    7.387,
    -0.079,
    -0.354,
]

# percentage infected values for an IC50 of ~47
perc_47 = [
    119.075,
    114.617,
    138.538,
    111.016,
    116.646,
    95.939,
    46.947,
    34.030,
]


perc_276 = [
    102.49,
    81.806,
    55.429,
    98.461,
    40.68,
    18.898,
    3.05,
    6.08,
]


perc_156 = [
    109.148,
    105.54,
    134.41,
    99.04,
    51.856,
    45.926,
    18.601,
    8.999,
]


perc_1675 = [
    59.712,
    58.014,
    13.587,
    6.521,
    0.272,
    0.0245,
    0.541,
    -0.10698,
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


def test_calc_model_results():
    expected_ic50 = 150
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test",
        df=df_good_inhibition,
        threshold=THRESHOLD,
        weak_threshold=WEAK_THRESHOLD,
    )
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < EPSILON
    assert mean_squared_error < MSE_PASS
    assert isinstance(model_params, stats.ModelParams)


def test_calc_391():
    expected_ic50 = 391
    df = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_391})
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test", df=df, threshold=THRESHOLD, weak_threshold=WEAK_THRESHOLD,
    )
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < EPSILON
    assert mean_squared_error < MSE_PASS
    assert isinstance(model_params, stats.ModelParams)


def test_calc_184():
    expected_ic50 = 184
    df = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_184})
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test", df=df, threshold=THRESHOLD, weak_threshold=WEAK_THRESHOLD,
    )
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < EPSILON
    assert mean_squared_error < MSE_PASS
    assert isinstance(model_params, stats.ModelParams)


def test_calc_381():
    expected_ic50 = 381
    df = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_381})
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test", df=df, threshold=THRESHOLD, weak_threshold=WEAK_THRESHOLD,
    )
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < EPSILON
    assert mean_squared_error < MSE_PASS
    assert isinstance(model_params, stats.ModelParams)


def test_calc_47():
    expected_ic50 = 47
    df = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_47})
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test", df=df, threshold=THRESHOLD, weak_threshold=WEAK_THRESHOLD,
    )
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < EPSILON
    assert mean_squared_error < MSE_PASS
    assert isinstance(model_params, stats.ModelParams)


def test_calc_276():
    expected_ic50 = 276.79
    df = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_276})
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test", df=df, threshold=THRESHOLD, weak_threshold=WEAK_THRESHOLD,
    )
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < EPSILON
    assert isinstance(model_params, stats.ModelParams)


def test_calc_156():
    expected_ic50 = 156.71
    df = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_156})
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test", df=df, threshold=THRESHOLD, weak_threshold=WEAK_THRESHOLD,
    )
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < EPSILON
    assert isinstance(model_params, stats.ModelParams)


def test_calc_1675():
    expected_ic50 = 1675
    df = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_1675})
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test", df=df, threshold=THRESHOLD, weak_threshold=WEAK_THRESHOLD,
    )
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < EPSILON
    assert mean_squared_error < MSE_PASS
    assert isinstance(model_params, stats.ModelParams)
