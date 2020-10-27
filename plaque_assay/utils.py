import os
import string


def row_col_to_well(row, col):
    """convert row-column indices (1-indexed) to well labels"""
    row_str = string.ascii_uppercase[row - 1]
    return f"{row_str}{col:02}"


def get_plate_num(plate_name):
    """get plate number from path"""
    basename = os.path.basename(plate_name)
    return int(basename.split("_")[0])
