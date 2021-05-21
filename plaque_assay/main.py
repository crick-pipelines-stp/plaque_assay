"""
Main `run()` function to launch the analysis.
"""


import os

import sqlalchemy

from plaque_assay.experiment import Experiment
from plaque_assay.errors import AlreadyUploadedError, DatabaseCredentialError
from plaque_assay import data
from plaque_assay import utils


def create_engine(test=True):
    """Create a database engine.

    Create a database engine to the LIMS serology database.

    This requires that the database credentials are in the user's
    environment. These credentials are:
    - `NE_USER` username
    - `NE_HOST_PROD` (if using production database) host
    - `NE_HOST_TEST` (if using testing database) host
    - `NE_PASSWORD` password


    Parameters
    ----------
    test : bool
        If True, then will use the 

    Returns
    -------
    sqlalchemy.engine.Engine
    """
    user = os.environ.get("NE_USER")
    if test:
        host = os.environ.get("NE_HOST_TEST")
    else:
        host = os.environ.get("NE_HOST_PROD")
    password = os.environ.get("NE_PASSWORD")
    if None in (user, host, password):
        raise DatabaseCredentialError(
            "db credentials not found in environent.",
            "Need to set NE_USER, NE_HOST_{TEST,PROD}, NE_PASSWORD",
        )
    engine = sqlalchemy.create_engine(f"mysql://{user}:{password}@{host}/serology")
    return engine


def create_local_engine():
    """Create a local database engine.

    Create a database engine to a local sqlite database.
    """
    engine = sqlalchemy.create_engine(
        "sqlite:////home/warchas/test_variant_plaque_assay_db.sqlite"
    )
    return engine


def run(plate_list):
    """Run analysis pipeline.

    This runs the entire analysis on a pair of plates, given that the 2
    plates are replicates for a single workflow_id and variant. The analysis
    will upload the results in the LIMS database.

    Parameters
    ------------
    plate_list : list
        List of paths to the plate directories to analyse. This will be 2 plates 
        for the 2 replicates for a single workflow and variant.

    Returns
    ----------
    None

    Raises
    -------
    AlreadyUploadedError
        This exception is raised if the given workflow and variant
        are already present in the LIMS serology database.
    """
    engine = create_engine(test=False)
    Session = sqlalchemy.orm.sessionmaker(bind=engine)
    session = Session()
    dataset = data.read_data_from_list(plate_list)
    indexfiles = data.read_indexfiles_from_list(plate_list)
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
    lims_db = data.DatabaseUploader(session)
    if lims_db.already_uploaded(workflow_id, variant):
        raise AlreadyUploadedError(
            f"workflow:{workflow_id} variant:{variant} already have results in the database"
        )
    lims_db.upload_plate_results(dataset)
    lims_db.upload_indexfiles(indexfiles)
    lims_db.upload_normalised_results(normalised_data)
    lims_db.upload_final_results(final_results)
    lims_db.upload_failures(failures)
    lims_db.upload_model_parameters(model_parameters)
    if lims_db.is_final_upload(workflow_id):
        lims_db.update_workflow_tracking(workflow_id)
    lims_db.commit()
