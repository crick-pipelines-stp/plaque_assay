from datetime import datetime, timezone

import pandas as pd

from ..db_uploader import BaseDatabaseUploader
from .. import db_models


class TitrationDatabaseUploader(BaseDatabaseUploader):
    """titration-specific database uploader"""

    def __init__(self, session):
        super().__init__(session)

    def already_uploaded(self, workflow_id: int) -> bool:
        """
        check if results have already been uploaded for a
        given workflow_id
        """
        raise NotImplementedError()

    def upload_normalised_results(self, normalised_results: pd.DataFrame) -> None:
        """Uploads normalised titration results to the
        LIMS database.

        Parameters
        ----------
        normalised_results: pd.DataFrame

        Returns
        -------
        None
            Uploads to the LIMS database.
        """
        # rename columns to match db
        # unpad wells
        # remove NaN/infs
        normalised_results = self.fix_for_mysql(normalised_results)
        # bulk insert mappings
        raise NotImplementedError()

    def upload_model_parameters(self, model_parameters: pd.DataFrame) -> None:
        """docstring"""
        # rename columns to match db
        # unpad wells
        # remove NaNs / infs
        model_parameters = self.fix_for_mysql(model_parameters)
        # bulk_insert_mappings
        raise NotImplementedError()

    def upload_final_results(self, final_results: pd.DataFrame) -> None:
        """docstring"""
        # rename columns to match db
        # unpad wells
        # remove NaN/infs
        final_results = self.fix_for_mysql(final_results)
        # bulk insert mappings
        raise NotImplementedError()

    def update_workflow_tracking(self, workflow_id: int) -> None:
        """Update NE_titration_workflow_tracking table to indicate the
        titration for a given workflow_id has been uploaded.

        Parameters
        ----------
        workflow_id: int
        """
        timestamp = datetime.now(timezone.utc)
        # fmt: off
        self.session\
            .query(db_models.NE_titration_workflow_tracking)\
            .filter(db_models.NE_titration_workflow_tracking.workflow_id == workflow_id)\
            .update(
                {
                    db_models.NE_titration_workflow_tracking.status: "complete",
                    db_models.NE_titration_workflow_tracking.end_date: timestamp
                }
            )
        # fmt: on
