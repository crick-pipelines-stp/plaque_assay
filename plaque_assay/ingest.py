"""
Data I/O
"""


import logging
import os
from typing import List

import pandas as pd

from . import utils
from . import consts


def read_data_from_list(plate_list: List) -> pd.DataFrame:
    """Read in data from plate list and assign dilution values by well position.

    Notes
    ------
    This will mock the data so the 4 dilutions on a single 384-well
    plate are re-labelled to appear from 4 different 96 well plates.

    Parameters
    ----------
    plate_list : list

    Returns:
    ---------
    pandas.DataFrame
    """
    dataframes = []
    barcodes = []
    for path in plate_list:
        df = pd.read_csv(
            # NOTE: might not always be Evaluation1
            os.path.join(path, "Evaluation1/PlateResults.txt"),
            skiprows=8,
            sep="\t",
        )
        plate_barcode = path.split(os.sep)[-1].split("__")[0]
        barcodes.append(plate_barcode)
        logging.info("plate barcode detected as %s", plate_barcode)
        well_labels = []
        for row, col in df[["Row", "Column"]].itertuples(index=False):
            well_labels.append(utils.row_col_to_well(row, col))
        df["Well"] = well_labels
        df["Plate_barcode"] = plate_barcode
        # Empty wells with no background produce NaNs rather than 0 in the
        # image analysis, which causes missing data for truely complete
        # inhbition. So we replace NaNs with 0 in the measurement columns we
        # use.
        fillna_cols = ["Normalised Plaque area", "Normalised Plaque intensity"]
        for colname in fillna_cols:
            df[colname] = df[colname].fillna(0)
        dataframes.append(df)
    df_concat = pd.concat(dataframes)
    # NOTE: mock barcodes before changing wells
    df_concat["Plate_barcode"] = utils.mock_384_barcode(
        existing_barcodes=df_concat["Plate_barcode"], wells=df_concat["Well"]
    )
    # mock wells
    df_concat["Well"] = [utils.well_384_to_96(i) for i in df_concat["Well"]]
    df_concat["PlateNum"] = [int(i[1]) for i in df_concat["Plate_barcode"]]
    df_concat["Dilution"] = [consts.plate_mapping[i] for i in df_concat["PlateNum"]]
    logging.debug("input data shape: %s", df_concat.shape)
    return df_concat


def get_plate_list(data_dir: str) -> List:
    """Get paths to plate directories

    Parameters
    ----------
    data_dir : str
        path to the directory containing the plate sub-directories

    Returns
    -------
    list
        list of 2 paths to plate directories
    """
    n_expected_plates = 2
    plate_list = [os.path.join(data_dir, i) for i in os.listdir(data_dir)]
    if len(plate_list) == n_expected_plates:
        logging.debug("plate list detected: %s", plate_list)
    else:
        logging.error(
            "Did not detect %s plates, detected %s :",
            n_expected_plates,
            len(plate_list),
        )
    return plate_list


def read_data_from_directory(data_dir: str) -> pd.DataFrame:
    """Read actual barcoded plate directory

    This gets a plate list from `data_dir`, and then
    reads in data from the plate list into a single
    DataFrame.

    Parameters
    -----------
    data_dir : str

    Returns
    -------
    pandas.DataFrame
    """
    plate_list = get_plate_list(data_dir)
    return read_data_from_list(plate_list)


def read_indexfiles_from_list(plate_list: List) -> pd.DataFrame:
    """Read indexfiles from a plate list

    Parameters
    ------------
    plate_list : list
        list of paths to plate directories

    Returns
    -------
    pandas.DataFrame
    """
    dataframes = []
    for path in plate_list:
        df = pd.read_csv(os.path.join(path, "indexfile.txt"), sep="\t")
        plate_barcode = path.split(os.sep)[-1].split("__")[0]
        df["Plate_barcode"] = plate_barcode
        dataframes.append(df)
    df_concat = pd.concat(dataframes)
    # remove annoying empty "Unnamed: 16" column
    to_rm = [col for col in df_concat.columns if col.startswith("Unnamed:")]
    df_concat.drop(to_rm, axis=1, inplace=True)
    if len(to_rm) > 0:
        logging.info("removed columns: %s from concatenated indexfiles", to_rm)
    logging.debug("indexfile shape: %s", df_concat.shape)
    return df_concat


def read_indexfiles_from_directory(data_dir: str) -> pd.DataFrame:
    """Return dataframe of indexfiles from a directory containing plates.

    Parameters
    -----------
    data_dir : str

    Returns
    --------
    pandas.DataFrame
    """
    plate_list = get_plate_list(data_dir)
    return read_indexfiles_from_list(plate_list)


def read_titration_data_from_list(plate_list: List[str]) -> pd.DataFrame:
    """Read in titration data from a plate_list,
    assigns dilution and sample info based on well position.

    Parameters
    -----------
    plate_list: list

    Returns
    --------
    pd.DataFrame
    """
    dataframes = []
    for path in plate_list:
        df = pd.read_csv(
            # NOTE: might not always be Evaluation1
            os.path.join(path, "Evaluation1/PlateResults.txt"),
            skiprows=8,
            sep="\t",
        )
        plate_barcode = path.split(os.sep)[-1].split("__")[0]
        logging.info("plate barcode detected as %s", plate_barcode)
        well_labels = []
        for row, col in df[["Row", "Column"]].itertuples(index=False):
            well_labels.append(utils.row_col_to_well(row, col))
        df["Well"] = well_labels
        df["Plate_barcode"] = plate_barcode
        # Empty wells with no background produce NaNs rather than 0 in the
        # image analysis, which causes missing data for truely complete
        # inhbition. So we replace NaNs with 0 in the measurement columns we
        # use.
        fillna_cols = ["Normalised Plaque area", "Normalised Plaque intensity"]
        for colname in fillna_cols:
            df[colname] = df[colname].fillna(0)
        dataframes.append(df)
    df_concat = pd.concat(dataframes)
    # remove empty wells
    df_concat = df_concat[~df_concat["Well"].isin(consts.TITRATION_EMPTY_WELLS)]
    # sample dilutions (1-4)
    dilution_int = [
        utils.titration_pos_control_dilution(well) for well in df_concat["Well"]
    ]
    # sample dilutions (40-2560)
    sample_dilution = [consts.plate_mapping.get(i) if i else None for i in dilution_int]
    df_concat["Dilution"] = sample_dilution
    # virus dilution factors (2-64)
    virus_dilutions = []
    for col_int in df_concat["Column"]:
        virus_dilution = consts.TITRATION_COLUMN_DILUTION_MAPPING.get(col_int)
        virus_dilutions.append(virus_dilution)
    df_concat["Virus_dilution_factor"] = virus_dilutions
    return df_concat
