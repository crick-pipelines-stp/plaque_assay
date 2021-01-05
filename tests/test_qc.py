import os
from collections import defaultdict

import pandas as pd

from plaque_assay import qc
from plaque_assay import utils


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(THIS_DIR, "test_data_files/RuthExpt20201019")


def load_data():
    """
    data loader, separate from the included package
    dataloader incase that changes.
    """
    experiments = [os.path.join(DATA_DIR, i) for i in os.listdir(DATA_DIR)]
    plate_name_dict = {os.path.abspath(i): utils.get_plate_num(i) for i in experiments}
    dataframes = []
    for path, plate_num in plate_name_dict.items():
        df = pd.read_csv(
            os.path.join(path, "Evaluation1/PlateResults.txt"), skiprows=8, sep="\t"
        )
        well_labels = []
        for row, col in df[["Row", "Column"]].itertuples(index=False):
            well = utils.row_col_to_well(row, col)
            well_labels.append(well)
        df["Well"] = well_labels
        df["PlateNum"] = plate_num
        dataframes.append(df)
    return pd.concat(dataframes)


def test_detect_low_cells_image_region_area():
    """
    test if we detect wells with low image_region_area in the DAPI
    channel, and fail plates if these failures occur in control wells
    """
    dataframe = load_data()
    failures = qc.detect_low_cells_image_region_area(dataframe, 0.7)
    # assert expected wells failed
    expected_failures_wells = {
        25: ["C05", "C06", "C07", "C08", "C10", "C11"],
        26: ["C06", "C07", "C08", "C10", "C11"],
    }
    assert failures.failed_wells == expected_failures_wells
    # assert no entire plate failures
    assert failures.failed_plates == []


def test_detect_high_background():
    """test if we detect wells with high-background fluorescence"""
    real_high_background_wells = set(["C06", "C07", "C08", "C10", "C11"])
    dataframe = load_data()
    failures = qc.detect_high_background(dataframe)
    # check we've found the high-background wells in both
    # the #25 and #26 plates
    high_background = defaultdict(set)
    for failure in failures:
        high_background[failure["plate"]].add(failure["well"])
    # check only plates #25 #26 have high-background wells
    assert sorted(list(high_background.keys())) == [25, 26]
    # check both plates have the expected failed wells
    for plate_num in high_background:
        assert high_background[plate_num] == real_high_background_wells