"""
Classes to hold information about plate and well failures.
"""
from typing import NamedTuple


CELL_IMAGE_AREA_FAILURE_REASON = "cell-image-region-area outside expected limits"
CELL_REGION_FAILURE_REASON = "cell-region-area outside expected range"
DAPI_PLATE_FAILURE_REASON = "possible plate fail - check DAPI plate image"


class PlateFailure(NamedTuple):
    plate: str
    well: str
    failure_reason: str
    failure_type: str = "plate_failure"


class WellFailure(NamedTuple):
    well: str
    plate: str
    failure_reason: str
    failure_type: str = "well_failure"
