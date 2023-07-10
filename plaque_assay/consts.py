VIRUS_ONLY_WELLS = ("A12", "B12", "C12", "B06", "C06", "D06", "E06", "F06", "G06")

NO_VIRUS_WELLS = ("F12", "G12", "H12")

POSITIVE_CONTROL_WELLS = ("D12", "E12", "A06", "H06")


UNWANTED_METADATA = [
    "Plane",
    "Timepoint",
    "Global Image Binning",
    "Height [µm]",
    "Time [s]",
    "Compound",
    "Concentration",
    "Cell Type",
    "Cell Count",
]

DILUTION_1 = 40
DILUTION_2 = 400
DILUTION_3 = 4000
DILUTION_4 = 40000


PLATE_MAPPING = {
    1: 1 / DILUTION_1,
    2: 1 / DILUTION_2,
    3: 1 / DILUTION_3,
    4: 1 / DILUTION_4,
}
