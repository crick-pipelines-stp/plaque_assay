"""
titration-specific utility functions
"""

from typing import Optional

from plaque_assay.titration import consts
from plaque_assay.utils import is_odd, is_even


def pos_control_dilution(well) -> Optional[int]:
    """
    Get dilution number from well label.
    NOTE: this is not the virus dilution factor which is positioned in pairs
          of columns, but the 4 dilutions within each dilution factor used
          to contruct the concentration-response curve.

               1 | 2
             +-------+
    G (odd)  | 4 | 2 |
    -------- +---+---+
    H (even) | 3 | 1 |
             +---+---+
    """
    if well not in consts.TITRATION_POSITIVE_CONTROL_WELLS:
        return None
    row_int = ord(well[0]) - 64  # 1 indexed
    col_int = int(well[1:])
    # positive control rows as G and H
    # as integers G = 7 = even
    #             H = 8 = odd
    if is_odd(row_int) and is_odd(col_int):
        return 4
    if is_even(row_int) and is_odd(col_int):
        return 3
    if is_odd(row_int) and is_even(col_int):
        return 2
    if is_even(row_int) and is_even(col_int):
        return 1
    return None
