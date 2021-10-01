from datetime import datetime, timezone

import pandas as pd

from plaque_assay.db_uploader import BaseDatabaseUploader
from plaque_assay import db_models
from plaque_assay.utils import unpad_well


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
            dataframe from `Titration.get_normalised_results()` with the
            columns:

            - plaque_area
            - background_subtracted_plaque_area
            - percentage_infected
            - dilution
            - well
            - plate_barcode
            - workflow_id

        Returns
        -------
        None
            Uploads to the LIMS database.
        """
        normalised_results.drop(columns=["Dilution", "dilution"], inplace=True)
        normalised_results = normalised_results.rename(columns={
            "Viral Plaques (global) - Area of Viral Plaques Area [µm²] - Mean per Well": "plaque_area"
        })
        # select and rename columns to match db
        new_colnames = [i.lower().replace(" ", "_") for i in normalised_results.columns]
        normalised_results.columns = new_colnames
        select_cols = [
            "plaque_area",
            "normalised_plaque_area",
            "background_subtracted_plaque_area",
            "well",
            "plate_barcode",
            "virus_dilution_factor",
            "percentage_infected",
            "workflow_id",
        ]
        normalised_results = normalised_results[select_cols].copy()
        normalised_results.rename(
            columns={"virus_dilution_factor": "dilution"}, inplace=True
        )
        # unpad wells
        normalised_results["well"] = [unpad_well(i) for i in normalised_results["well"]]
        # remove NaN/infs
        normalised_results = self.fix_for_mysql(normalised_results)
        # bulk insert mappings
        self.session.bulk_insert_mappings(
            db_models.NE_virus_titration_normalised_results,
            normalised_results.to_dict(orient="records"),
        )

    def upload_model_parameters(self, model_parameters: pd.DataFrame) -> None:
        """Upload model parameters to LIMS db

        Parameters
        -----------
        model_parameters: pd.DataFrame
            dataframe of model parameters from
            `Titration.get_model_parameters()` which has the columns:

            - dilution
            - param_top
            - param_bottom
            - param_ec50,
            - param_hillslope
            - mean_squared_error


        Returns
        -------
        None
            writes to database
        """
        model_parameters = self.fix_for_mysql(model_parameters)
        self.session.bulk_insert_mappings(
            db_models.NE_virus_titration_model_parameters,
            model_parameters.to_dict(orient="records"),
        )

    def upload_final_results(self, final_results: pd.DataFrame) -> None:
        """Upload final results to LIMS db.

        Parameters
        ----------
        final_results: pd.DataFrame
            dataframe of final results from `Titration.get_final_results()`
            with the columns :

            - dilution
            - ic50
            - status
            - workflow_id

        Returns
        -------
        None
            writes to database
        """
        # remove NaN/infs
        final_results = self.fix_for_mysql(final_results)
        # bulk insert mappings
        self.session.bulk_insert_mappings(
            db_models.NE_virus_titration_final_results,
            final_results.to_dict(orient="records"),
        )

    def update_workflow_tracking(self, workflow_id: int) -> None:
        """Update NE_titration_workflow_tracking table to indicate the
        titration for a given workflow_id has been uploaded.

        Parameters
        ----------
        workflow_id: int

        Returns
        --------
        None
            writes to database
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
