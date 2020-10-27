import os
import json
from glob import glob

import pandas as pd

import utils

### TEMP STUFF ###
# this will change as the barcodes are standardised
# but for now we are hard-coding this in the code

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

DATA_DIR = os.path.abspath("../data/RuthExpt20201014")

# NOTE: this is temporary until the barcodes are standardised
plate_mapping = {17: 1 / 40, 19: 1 / 160, 21: 1 / 640, 23: 1 / 2560}

VIRUS_ONLY_WELLS = ("A12", "B12", "C12")
NO_VIRUS_WELLS = ("F12", "G12", "H12")


def read_data(data_dir=DATA_DIR):
    """docstring"""
    # don't yet know how plates will be barcoded to identify the dilution series
    # for now just hard-code this information
    experiments = [os.path.join(DATA_DIR, i) for i in os.listdir(data_dir)]
    plate_name_dict = {os.path.abspath(i): utils.get_plate_num(i) for i in experiments}
    dataframes = []
    for path, plate_num in plate_name_dict.items():
        df = pd.read_csv(
            # NOTE: might not always be Evaluation2
            os.path.join(path, "Evaluation2/PlateResults.txt"),
            skiprows=8,
            sep="\t",
        )
        df["Dilution"] = plate_mapping[plate_num]
        # add well labels to dataframe
        well_labels = []
        for row, col in df[["Row", "Column"]].itertuples(index=False):
            well_labels.append(utils.row_col_to_well(row, col))
        df["Well"] = well_labels
        df["PlateNum"] = plate_num
        dataframes.append(df)
    return pd.concat(dataframes)


### END TEMP STUFF ###
