from plaque_assay import utils


THRESHOLD = 50
WEAK_THRESHOLD = 60


def test_col_to_well():
    assert utils.row_col_to_well(1, 1) == "A01"
    assert utils.row_col_to_well(8, 12) == "H12"


def test_well_384_to_96():
    assert utils.well_384_to_96("A01") == "A01"
    assert utils.well_384_to_96("P24") == "H12"
    assert utils.well_384_to_96("B02") == "A01"
    assert utils.well_384_to_96("P01") == "H01"


def test_get_dilution_from_384_well_label():
    assert utils.get_dilution_from_384_well_label("A01") == 4
    assert utils.get_dilution_from_384_well_label("A02") == 2
    assert utils.get_dilution_from_384_well_label("B01") == 3
    assert utils.get_dilution_from_384_well_label("B02") == 1
    assert utils.get_dilution_from_384_well_label("P24") == 1
    assert utils.get_dilution_from_384_well_label("A03") == 4


def test_mock_384_barcode():
    existing_barcodes = ["S01900001", "S02000001", "S02000001"]
    wells = ["A01", "A02", "A03"]
    output = utils.mock_384_barcode(existing_barcodes, wells)
    assert output == ["A41900001", "A22000001", "A42000001"]


def test_get_prefix_from_full_path():
    path = "/path/to/plate/S01000001__2021_01_01T00_00_00-Measuremment 1"
    assert utils.get_prefix_from_full_path(path) == "S01"
    path = "/path/to/plate/S12000001__2021_01_01T00_00_00-Measuremment 1"
    assert utils.get_prefix_from_full_path(path) == "S12"
