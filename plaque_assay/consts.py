VIRUS_ONLY_WELLS = ("A12", "B12", "C12", "B06", "C06", "D06", "E06", "F06", "G06")

NO_VIRUS_WELLS = ("F12", "G12", "H12")

POSITIVE_CONTROL_WELLS = ("D12", "E12", "A06", "H06")

plate_mapping = {1: 1 / 40, 2: 1 / 160, 3: 1 / 640, 4: 1 / 2560}

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

TITRATION_VIRUS_ONLY_ROWS = (
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
)

TITRATION_VIRUS_ONLY_WELLS = [
    f"{row}{col:02}" for row in TITRATION_VIRUS_ONLY_ROWS for col in range(1, 13)
]

TITRATION_NO_VIRUS_ROWS = ["P"]

TITRATION_NO_VIRUS_WELLS = set(
    f"{row}{col:02}" for row in TITRATION_NO_VIRUS_ROWS for col in range(1, 13)
)

TITRATION_POSITIVE_CONTROL_ROWS = ["G", "H"]

TITRATION_POSITIVE_CONTROL_WELLS = [
    f"{row}{col:02}" for row in TITRATION_POSITIVE_CONTROL_ROWS for col in range(1, 13)
]

TITRATION_COLUMN_DILUTION_MAPPING = {
    1: 2,
    2: 2,
    3: 4,
    4: 4,
    5: 8,
    6: 8,
    7: 16,
    8: 16,
    9: 32,
    10: 32,
    11: 64,
    12: 64,
}

TITRATION_EMPTY_COLUMNS = list(range(13, 25))

ALL_ROWS = [
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
]

TITRATION_EMPTY_WELLS = [
    f"{row}{col:02}" for row in ALL_ROWS for col in TITRATION_EMPTY_COLUMNS
]
