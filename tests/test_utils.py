import sqlalchemy

from plaque_assay import utils
from plaque_assay import db_models


THRESHOLD = 50
WEAK_THRESHOLD = 60

VARIANT_ENGLAND = "England2"
VARIANT_INDIA = "B.1.617.2 (India)"


def setup_module():
    """create in-memory sqlite Serology database for testing"""
    global engine
    global session
    engine = sqlalchemy.create_engine("sqlite://")
    Session = sqlalchemy.orm.sessionmaker(bind=engine)
    session = Session()
    # create tables from model definitions
    db_models.Base.metadata.create_all(engine)
    # add strains to testing database
    variant_england2 = db_models.NE_available_strains(
        mutant_strain="England2", plate_id_1="S01", plate_id_2="S02"
    )
    variant_delta = db_models.NE_available_strains(
        mutant_strain="B.1.617.2 (India)", plate_id_1="S09", plate_id_2="S10"
    )
    session.add(variant_delta)
    session.add(variant_england2)
    session.commit()


def teardown_module():
    pass


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


def test_get_variant_from_plate_list():
    plate_list_analysis_england2 = [
        "/mnt/NA_raw_data/S01000100__2021_01_01T01_01_01-Measurement 1",
        "/mnt/NA_raw_data/S02000100__2021_01_01T01_01_01-Measurement 1",
    ]
    plate_list_analysis_india = [
        "/mnt/NA_raw_data/S09000100__2021_01_01T01_01_01-Measurement 1",
        "/mnt/NA_raw_data/S10000100__2021_01_01T01_01_01-Measurement 1",
    ]
    plate_list_titration_england2 = [
        "/mnt/Titration_raw_data/T01000100__2021_01_01T01_01_01-Measurement 1",
        "/mnt/Titration_raw_data/T02000100__2021_01_01T01_01_01-Measurement 1",
    ]
    plate_list_titration_inda = [
        "/mnt/Titration_raw_data/T09000100__2021_01_01T01_01_01-Measurement 1",
        "/mnt/Titration_raw_data/T10000100__2021_01_01T01_01_01-Measurement 1",
    ]
    variant_analysis_eng2 = utils.get_variant_from_plate_list(
        plate_list_analysis_england2, session, titration=False
    )
    variant_analysis_india = utils.get_variant_from_plate_list(
        plate_list_analysis_india, session, titration=False
    )
    variant_titration_eng2 = utils.get_variant_from_plate_list(
        plate_list_titration_england2, session, titration=True
    )
    variant_titration_india = utils.get_variant_from_plate_list(
        plate_list_titration_inda, session, titration=True
    )
    #
    assert variant_analysis_eng2 == VARIANT_ENGLAND
    assert variant_analysis_india == VARIANT_INDIA
    assert variant_titration_eng2 == VARIANT_ENGLAND
    assert variant_titration_india == VARIANT_INDIA


#def test_titration_pos_control_dilution():
#    examples = [
#        ("G01", 4),
#        ("H01", 3),
#        ("G02", 2),
#        ("H02", 1),
#        ("G03", 4),
#        ("H04", 1),
#        ("A01", None),
#        ("H13", None),
#    ]
#    for input, expected_output in examples:
#        assert utils.titration_pos_control_dilution(input) == expected_output
