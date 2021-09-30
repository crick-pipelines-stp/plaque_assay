"""
titration-specific constants
"""


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
    1: 8,
    2: 8,
    3: 16,
    4: 16,
    5: 32,
    6: 32,
    7: 40,
    8: 40,
    9: 50,
    10: 50,
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
