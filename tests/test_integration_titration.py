import os
import datetime

import pandas as pd
import sqlalchemy

from plaque_assay.db_uploader import DatabaseUploader
from plaque_assay import ingest, db_models, utils, db_uploader
from plaque_assay.titration import Titration


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_MOCK = os.path.abspath(
    os.path.join(CURRENT_DIR, "test_data", "titration_data", "mock_data")
)
PLATE_LIST_MOCK = [os.path.join(TEST_DATA_MOCK, i) for i in os.listdir(TEST_DATA_MOCK)]


def setup_module():
    """create in-memory sqlite Serology database for testing"""
    global engine
    global session
    engine = sqlalchemy.create_engine("sqlite://")
    Session = sqlalchemy.orm.sessionmaker(bind=engine)
    session = Session()
    # create tables from model definitions
    db_models.Base.metadata.create_all(engine)
    # add England2 strain to testing database
    variant_england2 = db_models.NE_available_strains(
        mutant_strain="England2", plate_id_1="S01", plate_id_2="S02"
    )
    session.add(variant_england2)
    # add workflow_id 1 to workflow tracking table
    workflow_1 = db_models.NE_titration_workflow_tracking(
        plate_1="T01000001",
        plate_2="T02000001",
        variant="England2",
        start_date=datetime.datetime.utcnow() - datetime.timedelta(days=1),
        complete_operator="test operator",
        status="results",
        workflow_id=1
    )
    session.add(workflow_1)
    session.commit()
    run_titration_pipeline(PLATE_LIST_MOCK)


def teardown_module():
    pass


def run_titration_pipeline(plate_list):
    dataset = ingest.read_titration_data_from_list(plate_list)
    variant = utils.get_variant_from_plate_list(plate_list, session, titration=True)
    titration = Titration(dataset, variant=variant)
    workflow_id = titration.workflow_id
    results = titration.get_titration_results()
    lims_db = db_uploader.DatabaseUploader(session)
    lims_db.upload_titration_results(results)
    lims_db.update_titration_workflow_tracking(workflow_id, variant)
    lims_db.commit()


def test_titration_pipeline():
    query_results = (
        session
        .query(db_models.NE_virus_titration_results)
        .filter(db_models.NE_virus_titration_results.variant == "England2")
    )
    df = pd.read_sql(query_results.statement, con=session.bind)
    assert df.shape[0] > 0
    query_tracking = (
        session
        .query(db_models.NE_titration_workflow_tracking)
        .filter(db_models.NE_titration_workflow_tracking.variant == "England2")
    )
    df_tracking = pd.read_sql(query_tracking.statement, con=session.bind)
    assert df_tracking["status"].values[0] == "complete"
