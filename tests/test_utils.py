from plaque_assay import utils


THRESHOLD = 50
WEAK_THRESHOLD = 60


def test_col_to_well():
    assert utils.row_col_to_well(1, 1) == "A01"
    assert utils.row_col_to_well(8, 12) == "H12"


def test_get_plate_num():
    tests = [
        ("0000017__2020-10-16T14_46_58-Measurement 1", 17),
        ("/complete/path/to/0000017__2020-10-16T14_46_58-Measurement 1", 17),
        ("../relative/path/0000017__2020-10-16T14_46_58-Measurement 1", 17),
        ("0000999__2020-10-16T14_46_58-Measurement 1", 999),
    ]
    for test_string, expected_answer in tests:
        assert utils.get_plate_num(test_string) == expected_answer


def test_get_dilution_from_barcode():
    barcode = "A110000001__2020_00_00T00_00_00-Measurement 1"
    assert utils.get_dilution_from_barcode(barcode) == 1
    barcode = "/some/other/paths/A110000001__2020_00_00T00_00_00-Measurement 1"
    assert utils.get_dilution_from_barcode(barcode) == 1
    barcode = "/some/other/paths/A410000001__2020_00_00T00_00_00-Measurement 1"
    assert utils.get_dilution_from_barcode(barcode) == 4
