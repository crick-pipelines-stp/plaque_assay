import os
import shutil

from plaque_assay import main
from plaque_assay import utils


THIS_FILE = os.path.abspath(__file__)
THIS_DIR = os.path.dirname(THIS_FILE)
DATA_DIR = os.path.join(THIS_DIR, "test_data_files")
DATA_DIR_19 = os.path.join(THIS_DIR, "test_data_files/RuthExpt20201019")
DATA_DIR_21 = os.path.join(THIS_DIR, "test_data_files/RuthExpt20201021")


def setup_module():
    """create expected data_dir if it doesn't already exist"""
    rel_path = "../data/RuthExpt20201019"
    full_path = os.path.join(THIS_DIR, rel_path)
    if not os.path.isdir(full_path):
        shutil.copytree(DATA_DIR_19, full_path)
    rel_path = "../data/RuthExpt20201021"
    full_path = os.path.join(THIS_DIR, rel_path)
    if not os.path.isdir(full_path):
        shutil.copytree(DATA_DIR_21, full_path)


def test_main_dataset_19():
    """test main function simply runs without error"""
    results_19 = main.main(dataset=19)
    assert isinstance(results_19, dict)


def test_main_dataset_21():
    """test main function simply runs without error"""
    results_21 = main.main(dataset=21)
    assert isinstance(results_21, dict)


def test_result_sanity_dataset_19():
    """
    Want to determine that our results are sensible for
    a subset of wells.
    This will give a large margin-of-error and is **not** a
    test of reproducibility.
    """
    epsilon = 100
    results = main.main(dataset=19)
    expected_results = {
        "A01": utils.result_to_int("no inhibition"),
        "A02": utils.result_to_int("complete inhibition"),
        "C03": 220,
        "C11": 264,
        "E07": 55,
    }
    for well, expected_value in expected_results.items():
        obtained_value = results["results"][well]
        if isinstance(expected_value, str):
            assert obtained_value == expected_value
        elif isinstance(expected_value, (int, float)):
            if expected_value > 0 and obtained_value < 0:
                assert False, "returned a negative value (code for text)"
            assert abs(expected_value - obtained_value) < epsilon
        else:
            raise ValueError(
                f"unexpected type in expected_values {type(expected_value)}"
            )


def test_result_sanity_dataset_21():
    """
    Want to determine that our results are sensible for
    a subset of wells.
    This will give a large margin-of-error and is **not** a
    test of reproducibility.
    """
    epsilon = 100
    results = main.main(dataset=21)
    expected_results = {
        "A02": utils.result_to_int("complete inhibition"),
        "A05": 2100,
        "C07": utils.result_to_int("failed to fit model"),
        "E02": 150,
        "E03": 150,
        "E04": 110,
        "E05": 170,
        "E07": 100,
        "F12": utils.result_to_int("complete inhibition"),
    }
    for well, expected_value in expected_results.items():
        obtained_value = results["results"][well]
        if isinstance(expected_value, str):
            assert obtained_value == expected_value
        elif isinstance(expected_value, (int, float)):
            if expected_value > 0 and obtained_value < 0:
                assert False, "returned a negative value (code for text)"
            assert abs(expected_value - obtained_value) < epsilon
        else:
            raise ValueError(
                f"unexpected type in expected_values {type(expected_value)}"
            )


def teardown_module():
    """remove data_dir and any output files"""
    rel_path = "../data"
    full_path = os.path.join(THIS_DIR, rel_path)
    if os.path.isdir(full_path):
        shutil.rmtree(full_path)
    results_file_path = os.path.join(THIS_DIR, "../results.json")
    if os.path.exists(results_file_path):
        os.remove(results_file_path)
