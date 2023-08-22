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
    # G nanobody 1
    # H nanobody 1
    # I nanobody 2
    # J nanobody 2
    "K",
    "L",
    "M",
    "N",
    "O",
    # P no virus
)

TITRATION_VIRUS_ONLY_WELLS = [
    f"{row}{col:02}" for row in TITRATION_VIRUS_ONLY_ROWS for col in range(1, 25)
]

TITRATION_NO_VIRUS_ROWS = ["P"]

TITRATION_NO_VIRUS_WELLS = set(
    f"{row}{col:02}" for row in TITRATION_NO_VIRUS_ROWS for col in range(1, 25)
)

TITRATION_POSITIVE_CONTROL_ROWS = ["G", "H", "I", "J"]
TITRATION_NANOBODY_MAPPING = {
    "G": 1,
    "H": 1,
    "I": 2,
    "J": 2,
}

TITRATION_POSITIVE_CONTROL_WELLS = [
    f"{row}{col:02}" for row in TITRATION_POSITIVE_CONTROL_ROWS for col in range(1, 25)
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
    13: 6,
    14: 6,
    15: 12,
    16: 12,
    17: 24,
    18: 24,
    19: 48,
    20: 48,
    21: 96,
    22: 96,
    23: 192,
    24: 192,
}


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
