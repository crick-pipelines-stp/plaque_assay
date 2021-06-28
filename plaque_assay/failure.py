"""
Classes to hold information about plate and well failures.
"""
from typing import Dict


class PlateFailure:
    """
    Attributes
    -----------
    plate : str
    wells : list
    reason : str

    """

    def __init__(self, plate: str, wells: str, reason: str = ""):
        self.plate = plate
        self.wells = wells
        self.reason = reason

    def __str__(self):
        return str(self.to_dict())

    def __repr__(self):
        return str(self.to_dict())

    def to_dict(self) -> Dict:
        """convert to dictionary"""
        return {
            "type": "plate_failure",
            "plate": self.plate,
            "wells": self.wells,
            "reason": self.reason,
        }


class InfectionPlateFailure(PlateFailure):
    """Plate failure due to virus-only infection median outside expected range"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class CellAreaPlateFailure(PlateFailure):
    """Plate failure due to cell-image-region-area outside expected limits"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reason = "cell-image-region-area outside expected limits"


class DAPIPlateFailure(PlateFailure):
    """Plate failure due to large number of DAPI images outside expected range"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reason = "possible plate fail - check DAPI plate image"


class WellFailure:
    """Well failure class

    Attributes
    -----------
    well : str
    plate : str
    reason : str

    Methods
    -------
    to_dict()
        convert to dictionary
    """

    def __init__(self, well: str, plate: str, reason: str):
        self.well = well
        self.plate = plate
        self.reason = reason

    def __str__(self):
        return str(self.to_dict())

    def __repr__(self):
        return str(self.to_dict())

    def to_dict(self) -> Dict:
        return {
            "type": "well_failure",
            "plate": self.plate,
            "well": self.well,
            "reason": self.reason,
        }
