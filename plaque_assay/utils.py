import os
import string


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
    return int(basename[1])
