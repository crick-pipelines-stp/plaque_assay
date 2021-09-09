import os
import string
import math
from typing import List, Union, Optional

import pandas as pd
import sqlalchemy

from .db_models import NE_available_strains
from .errors import VariantLookupError
from . import consts


RESULT_TO_INT = {
    "failed to fit model": -999,
    "no inhibition": -600,
    "weak inhibition": -400,
    "complete inhibition": -200,
}


INT_TO_RESULT = {integer: result for result, integer in RESULT_TO_INT.items()}


def result_to_int(result: str) -> int:
    """convert result string to an integer

    Parameters
    -----------
    result : str

    Returns
    -------
    int
    """
    return RESULT_TO_INT[result]


def int_to_result(integer: int) -> str:
    """convert result integer to a string

    Parameters
    -----------
    result : int

    Returns
    -------
    str
    """
    return INT_TO_RESULT[integer]


def row_col_to_well(row: int, col: int) -> str:
    """join row and column indices to well labels

    Parameters
    -----------
    row : int
        integer row label (1-indexed)
    col : int
        integer row label (1-indexed)

    Returns
    --------
    str
        well label

    Examples
    ---------
    >>> row_col_to_well(2, 8)
    "B08"
    """
    row_str = string.ascii_uppercase[row - 1]
    return f"{row_str}{col:02}"


def unpad_well(well: str) -> str:
    """
    Remove zero-padding from well labels

    Arguments
    ----------
    well: str

    Returns
    --------
    str

    Examples
    --------
    >>> unpad_well("A01")
    "A1"
    """
    row = well[0]
    col = well[1:]
    return f"{row}{int(col)}"


def unpad_well_col(well_col: Union[List, pd.Series]) -> List:
    """Remove padding from an entire column of well labels

    Parameters
    -----------
    well_col : list or pandas.Series

    Returns
    --------
    list
        list of same well labels as input but without zero-padding
    """
    return [unpad_well(i) for i in well_col]


def well_384_to_96(well: str) -> str:
    """Convert 384 well label to 96 well label.

    If a 384-well plate has been constructed by stamping out 4 96-well plates,
    this function will convert the 384-well label to the well label of the
    original 96-well plate.

    Parameters
    ------------
    well : str
        alpha-numeric well label

    Returns
    --------
    str
        well label

    Examples
    --------
    >>> well_384_to_96("P24")
    "H12"

    """
    row = math.ceil((ord(well[0].upper()) - 64) / 2) - 1
    col = math.ceil(int(well[1:]) / 2)
    return f"{string.ascii_uppercase[row]}{col:02}"


def is_odd(n: int) -> bool:
    """is odd

    Parameters
    -----------
    n: int

    Returns
    --------
    Bool
    """
    return n % 2 != 0


def is_even(n: int) -> bool:
    """is even

    Parameters
    -----------
    n: int

    Returns
    --------
    Bool
    """
    return n % 2 == 0


def get_dilution_from_384_well_label(well: str) -> int:
    """Get dilution integer given a well label

    Given a 384-well plate constructed from 4 different diltuions
    of a 96-well plate stamped in quadrants, return the dilution integer
    [1, 2, 3, 4] from the well label.

    Parameters
    -----------
    well : str
        alpha-numeric well label

    Returns
    --------
    int
        dilution number [1, 2, 3 or 4]
    """
    row_int = ord(well[0].upper()) - 64
    col_int = int(well[1:])
    if is_odd(row_int) and is_odd(col_int):
        dilution = 4
    elif is_odd(row_int) and is_even(col_int):
        dilution = 2
    elif is_even(row_int) and is_odd(col_int):
        dilution = 3
    elif is_even(row_int) and is_even(col_int):
        dilution = 1
    else:
        raise RuntimeError()
    return dilution


def mock_384_barcode(existing_barcodes: List, wells: List) -> List:
    """384 -> 96 well plate barcodes

    Create mock barcodes for the mock 96-well plates
    coming from a 384-well plate

    Parameters
    -----------
    existing_barcodes : list
    wells : list

    Returns
    --------
    list
        new barcodes
    """
    new_barcodes = []
    for well, barcode in zip(wells, existing_barcodes):
        dilution_int = get_dilution_from_384_well_label(well)
        replicate_int = barcode[2]
        workflow_id = barcode[3:]
        new_barcode = f"A{dilution_int}{replicate_int}{workflow_id}"
        new_barcodes.append(new_barcode)
    return new_barcodes


def get_prefix_from_full_path(full_path: str) -> str:
    """Get prefix from full path

    Given a full-length path to a plate, this will return the
    prefix/barcode.

    Parameters
    -----------
    full_path : str
        full-length path to a plate directory

    Returns
    --------
    str
        prefix name
    """
    basename = os.path.basename(full_path)
    barcode = basename.split("__")[0]
    prefix = barcode[:3]
    return prefix


def get_variant_from_plate_list(
    plate_list: List, session: sqlalchemy.orm.Session, titration: bool = False
) -> str:
    """
    Fetch variant name from the LIMS database.
    This uses the plate barcodes within the plate_list which contain
    prefixes which are matched to a variant name in the NE_available_strains
    table in the LIMS database.

    Parameters
    -----------
    plate_list : list
        list of 2 full-length paths to plate directories
    session : sqlalchemy.orm.session.Session
        sqlalchemy sesssion to the LIMS serology database

    Returns
    --------
    str
        proper variant name from LIMS serology database

    Raises
    -------
    plaque_assay.errors.VariantError
        if plate barcode prefixes to not mach any known variant
        in the LIMS serology database
    """
    # get both prefixes from the plate_list
    assert len(plate_list) == 2, "expected plate_list to have 2 paths"
    prefixes = [get_prefix_from_full_path(i) for i in plate_list]
    if titration:
        # titration plates start with T rather than S
        # swap T for S at beginning i.e "T01" -> "S01"
        # so they match those in NE_available_strains
        prefixes = [i.replace("T", "S") for i in prefixes]
    prefix_1, prefix_2 = sorted(prefixes)
    # query table to return variant name (mutant_strain) for entry matching
    # both the plate barcode prefixes
    return_val = (
        session.query(NE_available_strains.mutant_strain)
        .filter(
            NE_available_strains.plate_id_1 == prefix_1,
            NE_available_strains.plate_id_2 == prefix_2,
        )
        .first()
    )
    if len(return_val) != 1:
        raise VariantLookupError(
            "plate barcode prefixes do not match any known variants in the ",
            f"LIMS database: {prefixes}",
        )
    return return_val.mutant_strain


def titration_pos_control_dilution(well) -> Optional[int]:
    """
    Get dilution number from well label.
    NOTE: this is not the virus dilution factor which is positioned in pairs
          of columns, but the 4 dilutions within each dilution factor used
          to contruct the concentration-response curve.

               1 | 2
             +-------+
    H (even) | 4 | 2 |
    -------- +---+---+
    I (odd)  | 3 | 1 |
             +---+---+
    """
    if well not in consts.TITRATION_POSITIVE_CONTROL_WELLS:
        return None
    row_int = ord(well[0]) - 64
    col_int = int(well[1:])
    # positive control rows as H and I
    # as integers H = 8 = even
    #             I = 9 = odd
    if is_odd(row_int) and is_odd(col_int):
        return 3
    if is_odd(row_int) and is_even(col_int):
        return 1
    if is_even(row_int) and is_odd(col_int):
        return 4
    if is_even(row_int) and is_even(col_int):
        return 2
    return None
