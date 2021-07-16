import os
from datetime import datetime, timezone, timedelta

import pandas as pd
import sqlalchemy

from plaque_assay import ingest, db_uploader, db_models, utils
from plaque_assay.experiment import Experiment


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_DIR_191_ENG2 = os.path.abspath(os.path.join(CURRENT_DIR, "test_data", "NA_raw_data_191_Eng2"))
TEST_DATA_DIR_228_DELTA = os.path.abspath(os.path.join(CURRENT_DIR, "test_data", "NA_raw_data_228_delta"))
PLATE_LIST_191_ENG2 = [os.path.join(TEST_DATA_DIR_191_ENG2, i) for i in os.listdir(TEST_DATA_DIR_191_ENG2)]
PLATE_LIST_228_DELTA = [os.path.join(TEST_DATA_DIR_228_DELTA, i) for i in os.listdir(TEST_DATA_DIR_228_DELTA)]


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
    session.add(variant_england2)
    session.add(variant_delta)
    # add workflows to testing database
    # just required columns
    workflow_191 = db_models.NE_workflow_tracking(
        master_plate="master_plate_string",
        start_date=datetime.now() - timedelta(days=1),
        no_of_variants=1,
        workflow_id=191,
    )
    workflow_228 = db_models.NE_workflow_tracking(
        master_plate="master_plate_string_2",
        start_date=datetime.now() - timedelta(days=1),
        no_of_variants=5,
        workflow_id=228,
    )
    session.add(workflow_191)
    session.add(workflow_228)
    session.commit()
    run_experiment(PLATE_LIST_191_ENG2)
    run_experiment(PLATE_LIST_228_DELTA)


def teardown_module():
    # remove local testing database
    pass


def run_experiment(plate_list):
    dataset = ingest.read_data_from_list(plate_list)
    indexfiles = ingest.read_indexfiles_from_list(plate_list)
    # add variant information to dataset and indexfiles dataframes
    variant = utils.get_variant_from_plate_list(plate_list, session)
    dataset["variant"] = variant
    indexfiles["variant"] = variant
    experiment = Experiment(dataset)
    workflow_id = int(experiment.experiment_name)
    normalised_data = experiment.get_normalised_data()
    final_results = experiment.get_results_as_dataframe()
    failures = experiment.get_failures_as_dataframe()
    model_parameters = experiment.get_model_parameters()
    lims_db = db_uploader.DatabaseUploader(session)
    lims_db.upload_plate_results(dataset)
    lims_db.upload_indexfiles(indexfiles)
    lims_db.upload_normalised_results(normalised_data)
    lims_db.upload_final_results(final_results)
    lims_db.upload_failures(failures)
    lims_db.upload_model_parameters(model_parameters)
    lims_db.upload_reporter_plate_status(workflow_id, variant)
    if lims_db.is_final_upload(workflow_id):
        lims_db.update_workflow_tracking(workflow_id)
    lims_db.commit()


def test_is_final_variant_upload():
    """
    test database is only expecting a single variant for workflow 191,
    so the England2 upload in `test_integration` should be marked as the
    final upload and have a timestamp
    """
    query = (
        session.query(db_models.NE_workflow_tracking)
        .filter(db_models.NE_workflow_tracking.workflow_id == 191)
        .first()
    )
    final_variant_upload = query.final_results_upload.replace(tzinfo=timezone.utc)
    assert datetime.now(timezone.utc) - final_variant_upload < timedelta(minutes=5)


def test_final_results():
    query = session.query(db_models.NE_final_results).filter(
        db_models.NE_final_results.workflow_id == 191
    )
    df_191 = pd.read_sql(query.statement, con=engine)
    assert df_191.shape[0] == 96


def test_model_parameters():
    query = session.query(db_models.NE_model_parameters).filter(
        db_models.NE_model_parameters.workflow_id == 191
    )
    df_191 = pd.read_sql(query.statement, con=engine)
    assert df_191.shape[0] == 96


def test_reporter_plate_status():
    query = (
        session.query(db_models.NE_reporter_plate_status)
        .filter(db_models.NE_model_parameters.workflow_id == 191)
        .first()
    )
    assert query.status == "awaiting"


def test_indexfiles():
    query = session.query(db_models.NE_raw_index).filter(
        db_models.NE_raw_index.workflow_id == 191
    )
    df_191 = pd.read_sql(query.statement, con=engine)
    assert df_191.shape[0] > 0


def test_plate_results():
    query = session.query(db_models.NE_raw_results).filter(
        db_models.NE_raw_results.workflow_id == 191
    )
    df_191 = pd.read_sql(query.statement, con=engine)
    assert df_191.shape[0] > 0


def test_normalised_results():
    query = session.query(db_models.NE_normalized_results).filter(
        db_models.NE_normalized_results.workflow_id == 191
    )
    df_191 = pd.read_sql(query.statement, con=engine)
    assert df_191.shape[0] > 0


def test_already_uploaded():
    lims_db = db_uploader.DatabaseUploader(session)
    variant = "England2"
    workflow_id = 191
    assert lims_db.already_uploaded(workflow_id, variant)


def test_failed_results_228():
    """
    228 Delta has lots of failures, check that the
    failure table is as expected.
    """
    query = session.query(db_models.NE_failed_results).filter(
        db_models.NE_failed_results.workflow_id == 228,
        db_models.NE_failed_results.variant == "B.1.617.2 (India)"
    )
    df_failures = pd.read_sql(query.statement, con=engine)
    # check we have some entries
    assert df_failures.shape[0] > 0
    # check we have both plate and well failures
    assert df_failures["failure_type"].nunique() == 2
    type_of_failures = df_failures["failure_type"].unique()
    for failure_type in ["well_failure", "plate_failure"]:
        assert failure_type in type_of_failures
    # check that there are no plate and well column mix-ups
    # can do this by checking that no well labels end with "000228"
    for i in df_failures["well"].values:
        assert not i.endswith("000228")
    for i in df_failures["plate"].values:
        assert i.endswith(("000228", "DILUTION SERIES"))


def test_failed_results_191():
    """check that the failure table is as expected"""
    query = session.query(db_models.NE_failed_results).filter(
        db_models.NE_failed_results.workflow_id == 191,
        db_models.NE_failed_results.variant == "England2"
    )
    df_failures = pd.read_sql(query.statement, con=engine)
    print(df_failures)
    # check we have some entries
    assert df_failures.shape[0] > 0
    # check we have both plate and well failures
    # check that there are no plate and well column mix-ups
    # can do this by checking that no well labels end with "000191"
    for i in df_failures["well"].values:
        assert not i.endswith("000191")
    for i in df_failures["plate"].values:
        assert i.endswith(("000191", "DILUTION SERIES"))