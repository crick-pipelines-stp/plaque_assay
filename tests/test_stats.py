from plaque_assay import stats
from plaque_assay import utils

import pandas as pd


THRESHOLD = 50
WEAK_THRESHOLD = 60


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
    fit_method, result, model_params, mean_squared_error = stats.calc_model_results(
        name="test",
        df=df_good_inhibition,
        threshold=THRESHOLD,
        weak_threshold=WEAK_THRESHOLD,
    )
    expected_ic50 = 150
    assert fit_method == "model fit"
    assert abs(result - expected_ic50) < 50
    assert mean_squared_error < 100
    assert isinstance(model_params, stats.ModelParams)
