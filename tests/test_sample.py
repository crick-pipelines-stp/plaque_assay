from plaque_assay.sample import Sample

import pandas as pd


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

# percentage infected for something that should be able to fit a model
# and have an IC50 value acceptable for a positive control well
# and have close replicates
perc_good = [
    100.556437,
    102.200186,
    50.246412,
    43.365569,
    23.787072,
    21.955933,
    12.517334,
    13.988952,
]


# percentage infected for something that should be able to fit a model
# but with 2 or more dodgy replicates
perc_bad_replicates = [
    100.556437,
    140.200186,
    100.246412,
    43.365569,
    23.787072,
    21.955933,
    12.517334,
    13.988952,
]

# some bad data which still triggers a model fit, but fails
perc_model_fail = [
    90,
    79,
    67,
    57,
    55,
    41,
    48,
    134
]

good_test_data = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_good}
)


bad_replicate_test_data = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_bad_replicates}
)

model_failure_data = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_model_fail}
)


def test_check_positive_control():
    sample_good = Sample(sample_name="A06", data=good_test_data)
    # shouldn't have any positive control failures
    assert len(sample_good.failures) == 0


def test_check_duplicate_differences():
    sample_bad_rep = Sample(sample_name="A01", data=bad_replicate_test_data)
    assert len(sample_bad_rep.failures) > 0
    # check that the failure is actually due to duplicate difference
    assert any(i.reason.startswith("2 or more duplicates differ") for i in sample_bad_rep.failures)


def test_check_for_model_fit_failure():
    sample = Sample(sample_name="A01", data=model_failure_data)
    assert len(sample.failures) > 0
    assert any(i.reason == "failed to fit model to data points" for i in sample.failures)
