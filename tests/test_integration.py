import os
import shutil

from plaque_assay import main
from plaque_assay import utils


THIS_FILE = os.path.abspath(__file__)
THIS_DIR = os.path.dirname(THIS_FILE)
DATA_DIR = os.path.join(THIS_DIR, "test_data_files")


def setup_module():
    """create expected data_dir if it doesn't already exist"""
    rel_path = "../data/RuthExpt20201019"
    full_path = os.path.join(THIS_DIR, rel_path)
    if not os.path.isdir(full_path):
        shutil.copytree(DATA_DIR, full_path)


def test_main():
    """test main function simply runs without error"""
    results = main.main()
    assert isinstance(results, dict)


def test_result_sanity():
    """
    Want to determine that our results are sensible for
    a subset of wells.
    This will give a large margin-of-error and is **not** a
    test of reproducibility.
    """
    epsilon = 100
    results = main.main()
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
