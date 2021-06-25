import os

import pandas as pd

from plaque_assay import data
from plaque_assay.experiment import Experiment


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "test_data", "NA_raw_data"))
PLATE_LIST = [os.path.join(TEST_DATA_DIR, i) for i in os.listdir(TEST_DATA_DIR)]


def test_integration():
    plate_list = PLATE_LIST
    dataset = data.read_data_from_list(plate_list)
    indexfiles = data.read_indexfiles_from_list(plate_list)
    # add variant information to dataset and indexfiles dataframes
    # normally this would involve a db query
    variant = "England2"
    dataset["variant"] = variant
    indexfiles["variant"] = variant
    experiment = Experiment(dataset)
    normalised_data = experiment.get_normalised_data()
    final_results = experiment.get_results_as_dataframe()
    failures = experiment.get_failures_as_dataframe()
    model_parameters = experiment.get_model_parameters()
    assert isinstance(normalised_data, pd.DataFrame)
    assert isinstance(final_results, pd.DataFrame)
    assert isinstance(failures, pd.DataFrame)
    assert isinstance(model_parameters, pd.DataFrame)
