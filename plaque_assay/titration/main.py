from typing import List

import sqlalchemy

from . import db_uploader
from . import ingest
from .. import utils
from ..main import create_engine
from .titration import Titration


def run(plate_list: List[str]) -> None:
    """Run titration analysis pipeline

    Parameters
    -----------
    plate_list : list

    Returns
    -------
    None
    """
    engine = create_engine(test=False)
    Session = sqlalchemy.orm.sessionmaker(bind=engine)
    session = Session()
    lims_db = db_uploader.TitrationDatabaseUploader(session)
    workflow_id = utils.get_workflow_id_from_plate_list(plate_list, titration=True)
    variant = utils.get_variant_from_plate_list(plate_list, titration=True)
    if lims_db.already_uploaded(workflow_id):
        print(f"workflow_id: {workflow_id} already have results in the database")
    dataset = ingest.read_data_from_list(plate_list)
    titration = Titration(dataset)
    normalised_results = titration.get_normalised_results()
    final_results = titration.get_final_results()
    model_parameters = titration.get_model_parameters()
    lims_db.upload_normalized_results(normalised_results)
    lims_db.upload_titration_final_results(final_results)
    lims_db.upload_model_parameters(model_parameters)
    lims_db.update_workflow_tracking(workflow_id=workflow_id, variant=variant)
    lims_db.commit()
