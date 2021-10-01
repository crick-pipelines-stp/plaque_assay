import os
import logging
from typing import List

import pandas as pd

from plaque_assay import utils
from plaque_assay import consts
from plaque_assay import titration


def read_data_from_list(plate_list: List[str]) -> pd.DataFrame:
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
    df_concat = df_concat[
        ~df_concat["Well"].isin(titration.consts.TITRATION_EMPTY_WELLS)
    ]
    # sample dilutions (1-4)
    dilution_int = [
        titration.utils.pos_control_dilution(well) for well in df_concat["Well"]
    ]
    # sample dilutions (40-2560)
    sample_dilution = [consts.plate_mapping.get(i) if i else None for i in dilution_int]
    df_concat["Dilution"] = sample_dilution
    # virus dilution factors (2-64)
    virus_dilutions = []
    for col_int in df_concat["Column"]:
        virus_dilution = titration.consts.TITRATION_COLUMN_DILUTION_MAPPING.get(col_int)
        virus_dilutions.append(virus_dilution)
    df_concat["Virus_dilution_factor"] = virus_dilutions
    return df_concat
