import os
import shutil

from plaque_assay import main


THIS_FILE = os.path.abspath(__file__)
THIS_DIR = os.path.dirname(THIS_FILE)
DATA_DIR = os.path.join(THIS_DIR, "test_data_files")


def setup_module():
    """create expected data_dir if it doesn't already exist"""
    rel_path = "../data/RuthExpt20201019"
    full_path = os.path.join(THIS_DIR, rel_path)
    if not os.path.isdir(full_path):
        shutil.copytree(DATA_DIR, full_path)


def teardown_module():
    """remove data_dir and any output files"""
    rel_path = "../data"
    full_path = os.path.join(THIS_DIR, rel_path)
    if os.path.isdir(full_path):
        shutil.rmtree(full_path)


def test_main():
    results = main.main()
    assert isinstance(results, dict)
