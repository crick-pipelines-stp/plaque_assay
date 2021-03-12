"""module docstring"""


class PlateFailure:
    """abstract plate class"""

    def __init__(self, plate, wells, reason=""):
        self.plate = plate
        self.wells = wells
        self.reason = reason

    def __str__(self):
        return str(self.to_dict())

    def __repr__(self):
        return str(self.to_dict())

    def to_dict(self):
        return {
            "type": "plate_failure",
            "plate": self.plate,
            "wells": self.wells,
            "reason": self.reason,
        }


class InfectionPlateFailure(PlateFailure):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class CellAreaPlateFailure(PlateFailure):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reason = "cell-image-region-area outside expected limits"


class DAPIPlateFailure(PlateFailure):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reason = "possible plate fail - check DAPI plate image"


class WellFailure:
    def __init__(self, well, plate, reason):
        self.well = well
        self.plate = plate
        self.reason = reason

    def __str__(self):
        return str(self.to_dict())

    def __repr__(self):
        return str(self.to_dict())

    def to_dict(self):
        return {
            "type": "well_failure",
            "plate": self.plate,
            "well": self.well,
            "reason": self.reason,
        }
