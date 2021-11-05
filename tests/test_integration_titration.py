import os
import datetime

import pandas as pd
import sqlalchemy

import plaque_assay
from plaque_assay import db_models
from plaque_assay import titration
from plaque_assay import titration_db_uploader
from plaque_assay import titration_ingest


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
N_WELLS = 384  # half plate (16 rows, 12 cols)
N_REPLICATES = 2
N_DILUTIONS = 12
N_NANOBODIES = 2


def create_plate_list(dir_name):
    paths = os.path.abspath(
        os.path.join(CURRENT_DIR, "test_data", "titration_data", "real_data", dir_name)
    )
    return [os.path.join(paths, i) for i in os.listdir(paths)]


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
    variant_delta = db_models.NE_available_strains(
        mutant_strain="B.1.617.2 (India)", plate_id_1="S09", plate_id_2="S10"
    )
    variant_b117 = db_models.NE_available_strains(
        mutant_strain="B117", plate_id_1="S03", plate_id_2="S04"
    )
    variant_b1351 = db_models.NE_available_strains(
        mutant_strain="B1351", plate_id_1="S05", plate_id_2="S06"
    )
    session.add(variant_delta)
    session.add(variant_b117)
    session.add(variant_b1351)
    # add workflow_id 1 to workflow tracking table
    workflow_delta_15 = db_models.NE_titration_workflow_tracking(
        plate_1="T09000015",
        plate_2="T10000015",
        variant="B.1.617.2 (India)",
        start_date=datetime.datetime.utcnow() - datetime.timedelta(days=1),
        complete_operator="test operator",
        status="results",
        workflow_id=15,
    )
    workflow_b117_8 = db_models.NE_titration_workflow_tracking(
        plate_1="T03000008",
        plate_2="T04000008",
        variant="B117",
        start_date=datetime.datetime.utcnow() - datetime.timedelta(days=1),
        complete_operator="test operator",
        status="results",
        workflow_id=8,
    )
    workflow_b1351_7 = db_models.NE_titration_workflow_tracking(
        plate_1="T05000007",
        plate_2="T06000007",
        variant="B1351",
        start_date=datetime.datetime.utcnow() - datetime.timedelta(days=1),
        complete_operator="test operator",
        status="results",
        workflow_id=7,
    )
    session.add(workflow_delta_15)
    session.add(workflow_b117_8)
    session.add(workflow_b1351_7)
    session.commit()
    plate_list_delta_15 = create_plate_list("delta_000015")
    plate_list_b117_8 = create_plate_list("B117_000008")
    plate_list_b1351_7 = create_plate_list("B1351_000007")
    run_titration_pipeline(plate_list_delta_15)
    run_titration_pipeline(plate_list_b117_8)
    run_titration_pipeline(plate_list_b1351_7)


def teardown_module():
    pass


def run_titration_pipeline(plate_list):
    dataset = titration_ingest.read_data_from_list(plate_list)
    variant = plaque_assay.utils.get_variant_from_plate_list(
        plate_list, session, titration=True
    )
    lims_db = titration_db_uploader.TitrationDatabaseUploader(session)
    workflow_id = plaque_assay.utils.get_workflow_id_from_plate_list(plate_list)
    if lims_db.already_uploaded(workflow_id):
        print(f"workflow_id: {workflow_id} already have results in the database")
        # still exist successfully so task is marked complete
        return None
    titration_cl = titration.Titration(dataset, variant=variant)
    normalised_results = titration_cl.get_normalised_results()
    final_results = titration_cl.get_final_results()
    model_parameters = titration_cl.get_model_parameters()
    lims_db.upload_final_results(final_results)
    lims_db.upload_model_parameters(model_parameters)
    lims_db.upload_normalised_results(normalised_results)
    lims_db.update_workflow_tracking(workflow_id=workflow_id)
    lims_db.commit()


def pipeline_output(workflow_id, variant_name):
    query_normalised_results = session.query(
        db_models.NE_virus_titration_normalised_results
    ).filter(db_models.NE_virus_titration_normalised_results.workflow_id == workflow_id)
    df_normalised_results = pd.read_sql(
        query_normalised_results.statement, con=session.bind
    )
    assert df_normalised_results.shape[0] == N_WELLS * N_REPLICATES
    ###
    query_model_parameters = session.query(
        db_models.NE_virus_titration_model_parameters
    ).filter(db_models.NE_virus_titration_model_parameters.workflow_id == workflow_id)
    df_model_parameters = pd.read_sql(
        query_model_parameters.statement, con=session.bind
    )
    assert df_model_parameters.shape[0] == N_DILUTIONS * N_NANOBODIES
    ##
    query_final_results = session.query(
        db_models.NE_virus_titration_final_results
    ).filter(db_models.NE_virus_titration_final_results.workflow_id == workflow_id)
    df_final_results = pd.read_sql(query_final_results.statement, con=session.bind)
    assert df_final_results.shape[0] == N_DILUTIONS * N_NANOBODIES
    ##
    query_tracking = session.query(db_models.NE_titration_workflow_tracking).filter(
        db_models.NE_titration_workflow_tracking.workflow_id == workflow_id,
        db_models.NE_titration_workflow_tracking.variant == variant_name,
    )
    df_tracking = pd.read_sql(query_tracking.statement, con=session.bind)
    assert df_tracking["status"].values[0] == "complete"


def test_pipeline_output():
    params = [
        [15, "B.1.617.2 (India)"],
        [8, "B117"],
        [7, "B1351"],
    ]
    for workflow_id, variant_name in params:
        pipeline_output(workflow_id, variant_name)
