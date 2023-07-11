from plaque_assay.sample import Sample
from plaque_assay.failure import WellFailure
from plaque_assay import consts

import pandas as pd


VARIANT = "England2"

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


perc_high_mse = [
    180.556437,
    80.200186,
    60.246412,
    30.365569,
    11.787072,
    30.955933,
    3.517334,
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
perc_model_fail = [90, 79, 67, 57, 55, 41, 48, 134]


perc_bad_no_inhibition = [
    100.556437,
    140.200186,
    100.246412,
    80.365569,
    100.787072,
    160.955933,
    100.517334,
    90.988952,
]

# good data, but IC50 value outside expected range
# for the B117 variant's positive control
perc_wrong_ic50 = [
    64.556437,
    60.200186,
    30.246412,
    23.365569,
    13.787072,
    11.955933,
    12.517334,
    13.988952,
]


good_test_data = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_good})


bad_replicate_test_data = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_bad_replicates}
)

model_failure_data = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_model_fail}
)


bad_replicate_no_inhib_data = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_bad_no_inhibition}
)


good_but_wrong_ic50_for_pos_cntrl = pd.DataFrame(
    {"Dilution": dilutions, "Percentage Infected": perc_wrong_ic50}
)

high_mse = pd.DataFrame({"Dilution": dilutions, "Percentage Infected": perc_high_mse})


def test_check_positive_control(monkeypatch):
    monkeypatch.setattr(consts, "DILUTION_1", 40)
    monkeypatch.setattr(consts, "DILUTION_2", 160)
    monkeypatch.setattr(consts, "DILUTION_3", 480)
    monkeypatch.setattr(consts, "DILUTION_4", 2560)
    sample_good = Sample(sample_name="A06", data=good_test_data, variant=VARIANT)
    # shouldn't have any positive control failures
    assert len(sample_good.failures) == 0


def test_check_positive_control_failure():
    sample_wrong_ic50 = Sample(
        sample_name="A06", data=good_but_wrong_ic50_for_pos_cntrl, variant="B117"
    )
    # check there's some failures
    assert len(sample_wrong_ic50.failures) > 0
    # check they're actually failures for unexpected positive control values
    assert any(
        i.failure_reason.startswith("positive control failure")
        for i in sample_wrong_ic50.failures
    )
    # FIXME: pos control ic50 criteria has changed since these samples
    #        need to update the values used in the test
    # use the india (delta) variant, should pass
    #sample_right_ic50 = Sample(
    #    sample_name="A06",
    #    data=good_but_wrong_ic50_for_pos_cntrl,
    #    variant="B.1.617.2 (India)",
    #)
    #assert len(sample_right_ic50.failures) == 0


def test_check_duplicate_differences(monkeypatch):
    monkeypatch.setattr(consts, "DILUTION_1", 40)
    monkeypatch.setattr(consts, "DILUTION_2", 160)
    monkeypatch.setattr(consts, "DILUTION_3", 480)
    monkeypatch.setattr(consts, "DILUTION_4", 2560)
    sample_bad_rep = Sample(
        sample_name="A01", data=bad_replicate_test_data, variant=VARIANT
    )
    assert len(sample_bad_rep.failures) > 0
    # check that the failure is actually due to duplicate difference
    assert any(
        i.failure_reason.startswith("2 or more duplicates differ")
        for i in sample_bad_rep.failures
    )


def test_check_duplicate_differnces_no_inhibition():
    """Don't flag duplicate difference if the result is "no inhibition"""
    sample = Sample(
        sample_name="A01", data=bad_replicate_no_inhib_data, variant=VARIANT
    )
    assert len(sample.failures) == 0
    assert sample.ic50_pretty == "no inhibition"


def test_check_for_model_fit_failure():
    sample = Sample(sample_name="A01", data=model_failure_data, variant=VARIANT)
    assert len(sample.failures) > 0
    assert any(
        i.failure_reason == "failed to fit model to data points"
        for i in sample.failures
    )


def test_check_for_high_mse(monkeypatch):
    monkeypatch.setattr(consts, "DILUTION_1", 40)
    monkeypatch.setattr(consts, "DILUTION_2", 160)
    monkeypatch.setattr(consts, "DILUTION_3", 480)
    monkeypatch.setattr(consts, "DILUTION_4", 2560)
    sample = Sample(sample_name="A01", data=high_mse, variant=VARIANT)
    assert len(sample.failures) > 0
    failure = list(sample.failures)[0]
    assert isinstance(failure, WellFailure)
    assert failure.failure_reason.startswith("model MSE")
