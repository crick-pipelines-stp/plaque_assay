import logging
import os
import string
import math


RESULT_TO_INT = {
    "failed to fit model": -999,
    "no inhibition": -600,
    "weak inhibition": -400,
    "complete inhibition": -200,
}


INT_TO_RESULT = {integer: result for result, integer in RESULT_TO_INT.items()}


def result_to_int(result):
    """convert result string to an integer"""
    return RESULT_TO_INT[result]


def int_to_result(integer):
    """convert result integer to a string"""
    return INT_TO_RESULT[integer]


def row_col_to_well(row, col):
    """convert row-column indices (1-indexed) to well labels"""
    row_str = string.ascii_uppercase[row - 1]
    return f"{row_str}{col:02}"


def get_plate_num(plate_name):
    """get plate number from path"""
    basename = os.path.basename(plate_name)
    return int(basename.split("_")[0])


def get_dilution_from_barcode(plate_name):
    """get dilution integer [1, 2, 3, 4] from full plate path"""
    basename = os.path.basename(plate_name)
    try:
        return int(basename[1])
    except ValueError as e:
        logging.error("Failed to parse barcode from plate %s, error: %s", plate_name, e)
        raise e


def unpad_well(well):
    """
    Remove zero-padding from well labels

    Arguments:
    ----------
    well: string

    Returns:
    --------
    string
    """
    row = well[0]
    col = well[1:]
    return f"{row}{int(col)}"


def unpad_well_col(well_col):
    """
    Remove padding from an entire column of well labels
    """
    return [unpad_well(i) for i in well_col]


def well_384_to_96(well):
    """
    If a 384-well plate has been constructed by stamping out 4 96-well plates,
    this function will convert the 384-well label to the well label of the
    original 96-well plate.
    e.g well_384_to_96("P24) == "H12"
    """
    row = math.ceil((ord(well[0].upper()) - 64) / 2) - 1
    col = math.ceil(int(well[1:]) / 2)
    return f"{string.ascii_uppercase[row]}{col:02}"


def is_odd(n):
    return n % 2 != 0


def is_even(n):
    return n % 2 == 0


def get_dilution_from_384_well_label(well):
    """
    Given a 384-well plate constructed from 4 different diltuions
    of a 96-well plate stamped in quadrants, return the dilution integer
    [1, 2, 3, 4] from the well label.
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


def mock_384_barcode(existing_barcodes, wells):
    """
    create mock barcodes for the mock 96-well plates
    coming from a 384-well plate
    """
    new_barcodes = []
    for well, barcode in zip(wells, existing_barcodes):
        dilution_int = get_dilution_from_384_well_label(well)
        new_barcode = barcode[:1] + str(dilution_int) + barcode[2:]
        new_barcodes.append(new_barcode)
    return new_barcodes
