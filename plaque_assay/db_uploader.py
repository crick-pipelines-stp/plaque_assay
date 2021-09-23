import logging
from datetime import datetime, timezone

import pandas as pd
import numpy as np

from plaque_assay import db_models
from plaque_assay import utils


class DatabaseUploader:
    """
    Attributes
    ----------
    session : sql.orm.session.Session
        sqlalchemy session to LIMS serology database

    """

    def __init__(self, session):
        self.session = session

    @staticmethod
    def fix_for_mysql(df: pd.DataFrame) -> pd.DataFrame:
        """
        Replace NaN and inf values with None, as these cannot be used
        with MySQL and None values are converted to null

        Parameters
        ----------
        df : pd.DataFrame

        Returns
        --------
        pd.DataFrame
        """
        return df.replace({np.inf: None, -np.inf: None}).replace({np.nan: None})

    def commit(self) -> None:
        """commit data to LIMS serology database"""
        self.session.commit()

    def already_uploaded(
        self, workflow_id: int, variant: str, titration: bool = False
    ) -> bool:
        """
        Check if the results for a given workflow_id and variant
        have already been uploaded.

        As results are uploaded all-or-nothing, we can just check one of the
        results tables for the presence of the supplied workflow_id and
        variant, rather than each of the several results tables.

        Parameters
        -----------
        workflow_id : int
        variant : str
        titration: bool (default = False)

        Returns
        --------
        bool
        """
        if titration:
            table = db_models.NE_virus_titration_results
        else:
            table = db_models.NE_final_results
        result = (
            self.session.query(table)
            .filter(table.workflow_id == workflow_id, table.variant == variant,)
            .first()
        )
        return result is not None

    def is_final_upload(self, workflow_id: int) -> bool:
        """
        This determines if a given workflow_id is the last to be
        uploaded for that workflow_id. This is determined by the
        number of unique variants for a workflow_id currently in the
        LIMS database. If this is 1 short of the expected number, then
        we can determine the current upload is the final one.

        Parameters
        -----------
        workflow_id: int

        Returns
        ---------
        bool
            whether or not this is the final upload

        Raises
        ------
        RuntimeError
            when trying to upload more variants than specified in the
            workflow_tracking LIMS database table.
        """
        # get expected number of variants from NE_workflow_tracking
        # fmt: off
        expected = (
            self.session
            .query(db_models.NE_workflow_tracking)
            .filter(db_models.NE_workflow_tracking.workflow_id == workflow_id)
            .first()
        )
        # get number of uploaded variants from NE_final_results
        current_variants = (
            self.session
            .query(db_models.NE_final_results.variant)
            .filter(db_models.NE_final_results.workflow_id == workflow_id)
        )
        # fmt: on
        current_n_variants = current_variants.distinct().count()
        # NOTE: the sqlalchemy queries are reading from NE_final_results data
        # that includes results from this session that have not yet been
        # committed, and as is_final_upload() is called *after*
        # upload_final_results(), we pretend the results are already in the
        # database.
        is_final = int(expected.no_of_variants) == int(current_n_variants)
        if int(current_n_variants) > int(expected.no_of_variants):
            raise RuntimeError(
                f"unexpected no. of variants {current_n_variants}, expecting max of {expected.no_of_variants}"
            )
        logging.debug(f"expected no. of variants: {expected.no_of_variants}")
        logging.debug(
            f"current no. uploaded variants: {current_n_variants}: {current_variants}"
        )
        if is_final:
            logging.info(
                f"Final variant upload, marking workflow {workflow_id} as complete"
            )
        else:
            logging.info(
                f"Not final variant upload for workflow {workflow_id}, this is variant {current_n_variants}/{expected.no_of_variants}"
            )
        return is_final

    def upload_plate_results(self, plate_results_dataset: pd.DataFrame) -> None:
        """Upload raw concatenated data into database.

        This uploads the "raw" dataset into the LIMS serology database.
        The data is largely what is exported from the Phenix with columns
        renamed, some columns moved, some metadata columns added.

        Parameters
        ----------
        plate_results_dataset: pandas.DataFrame

        Returns:
        ---------
        None
        """
        # TODO: check csv matches master plate selected in NE_workflow_tracking
        plate_results_dataset = plate_results_dataset.copy()
        rename_dict = {
            "Row": "row",
            "Column": "column",
            "Viral Plaques (global) - Area of Viral Plaques Area [µm²] - Mean per Well": "VPG_area_mean",
            "Viral Plaques (global) - Intensity Viral Plaques Alexa 488 (global) Mean - Mean per Well": " VPG_intensity_mean_per_well",
            "Viral Plaques (global) - Intensity Viral Plaques Alexa 488 (global) StdDev - Mean per Well": "VPG_intensity_stddev_per_well",
            "Viral Plaques (global) - Intensity Viral Plaques Alexa 488 (global) Median - Mean per Well": "VPG_intensity_median_per_well",
            "Viral Plaques (global) - Intensity Viral Plaques Alexa 488 (global) Sum - Mean per Well": "VPG_intensity_sum_per_well",
            "Cells - Intensity Image Region DAPI (global) Mean - Mean per Well": "cells_intensity_mean_per_well",
            "Cells - Intensity Image Region DAPI (global) StdDev - Mean per Well": "cells_intensity_stddev_mean_per_well",
            "Cells - Intensity Image Region DAPI (global) Median - Mean per Well": "cells_intensity_median_mean_per_well",
            "Cells - Intensity Image Region DAPI (global) Sum - Mean per Well": "cells_intensity_sum_mean_per_well",
            "Cells - Image Region Area [µm²] - Mean per Well": "cells_image_region_area_mean_per_well",
            "Normalised Plaque area": "normalised_plaque_area",
            "Normalised Plaque intensity": "normalised_plaque_intensity",
            "Number of Analyzed Fields": "number_analyzed_fields",
            "Dilution": "dilution",
            "Well": "well",
            "PlateNum": "plate_num",
            "Plate_barcode": "plate_barcode",
            "variant": "variant",
            # "Background Subtracted Plaque Area": "background_subtracted_plaque_area",
        }
        plate_results_dataset.rename(columns=rename_dict, inplace=True)
        # filter to only desired columns
        plate_results_dataset = plate_results_dataset[list(rename_dict.values())]
        workflow_id = [int(i[3:]) for i in plate_results_dataset["plate_barcode"]]
        plate_results_dataset["workflow_id"] = workflow_id
        plate_results_dataset["well"] = utils.unpad_well_col(
            plate_results_dataset["well"]
        )
        plate_results_dataset = self.fix_for_mysql(plate_results_dataset)
        self.session.bulk_insert_mappings(
            db_models.NE_raw_results, plate_results_dataset.to_dict(orient="records")
        )

    def upload_indexfiles(self, indexfiles_dataset: pd.DataFrame) -> None:
        """Upload indexfiles from the Phenix into the database

        This uploads the IndexFile dataset into the LIMS serology database,
        which is kept because it contains the URLs to the images.

        Some columns are removed, most are renamed, and a variant metadata
        column is added.

        Parameters
        ----------
        indexfiles_dataset : pandas.DataFrame

        Returns
        -------
        None
        """
        indexfiles_dataset = indexfiles_dataset.copy()
        rename_dict = {
            "Row": "row",
            "Column": "column",
            "Field": "field",
            "Channel ID": "channel_id",
            "Channel Name": "channel_name",
            "Channel Type": "channel_type",
            "URL": "url",
            "ImageResolutionX [m]": "image_resolutionx",
            "ImageResolutionY [m]": "image_resolutiony",
            "ImageSizeX": "image_sizex",
            "ImageSizeY": "image_sizey",
            "PositionX [m]": "positionx",
            "PositionY [m]": "positiony",
            "Time Stamp": "time_stamp",
            "Plate_barcode": "plate_barcode",
            "variant": "variant",  # not renamed, just to keep it
        }
        indexfiles_dataset.rename(columns=rename_dict, inplace=True)
        # filter to only desired columns
        indexfiles_dataset = indexfiles_dataset[list(rename_dict.values())]
        # get workflow ID
        workflow_id = [int(i[3:]) for i in indexfiles_dataset["plate_barcode"]]
        indexfiles_dataset["workflow_id"] = workflow_id
        indexfiles_dataset = self.fix_for_mysql(indexfiles_dataset)
        for i in range(0, len(indexfiles_dataset), 1000):
            df_slice = indexfiles_dataset.iloc[i : i + 1000]
            self.session.bulk_insert_mappings(
                db_models.NE_raw_index, df_slice.to_dict(orient="records")
            )

    def upload_normalised_results(self, norm_results: pd.DataFrame) -> None:
        """Upload normalised results into the database.

        Uploads the normalised data, consisting of
        `background_subtracted_plaque_area`, `percentage_infected` and
        metadata, into the LIMS serology database.

        Parameters
        ----------
        norm_results: pandas.DataFrame

        Returns
        -------
        None
        """
        norm_results = norm_results.copy()
        rename_dict = {
            "Well": "well",
            "Row": "row",
            "Column": "column",
            "Dilution": "dilution",
            "Plate_barcode": "plate_barcode",
            "Background_subtracted_plaque_area": "background_subtracted_plaque_area",
            "Percentage_infected": "percentage_infected",
            "variant": "variant",  # not renamed, just to keep
        }
        norm_results.rename(columns=rename_dict, inplace=True)
        norm_results = norm_results[list(rename_dict.values())]
        workflow_id = [int(i[3:]) for i in norm_results["plate_barcode"]]
        assert len(set(workflow_id)) == 1
        norm_results["workflow_id"] = workflow_id
        norm_results["well"] = utils.unpad_well_col(norm_results["well"])
        norm_results = self.fix_for_mysql(norm_results)
        self.session.bulk_insert_mappings(
            db_models.NE_normalized_results, norm_results.to_dict(orient="records")
        )

    def upload_final_results(self, results: pd.DataFrame) -> None:
        """Upload final results to database

        Final results are mainly IC50 and metadata.

        Parameters
        -----------
        results : pandas.DataFrame

        Returns
        -------
        None
        """
        results = results.copy()
        # don't have master_plate details from accessible tables, set as None
        results["master_plate"] = None
        # get workflow_id
        assert results["experiment"].nunique() == 1
        assert results["variant"].nunique() == 1
        results["workflow_id"] = results["experiment"].astype(int)
        results["well"] = utils.unpad_well_col(results["well"])
        results = self.fix_for_mysql(results)
        self.session.bulk_insert_mappings(
            db_models.NE_final_results, results.to_dict(orient="records")
        )

    def upload_failures(self, failures: pd.DataFrame) -> None:
        """Upload failure information to database

        Parameters
        ----------
        failures : pandas.DataFrame

        Returns
        -------
        None
        """
        failures = failures.copy()
        if failures.shape[0] > 0:
            assert failures["experiment"].nunique() == 1
            failures["workflow_id"] = failures["experiment"].astype(int)
            self.session.bulk_insert_mappings(
                db_models.NE_failed_results, failures.to_dict(orient="records")
            )

    def upload_model_parameters(self, model_parameters: pd.DataFrame) -> None:
        """Upload model parameters to database

        Parameters
        -----------
        model_parameters: pandas.DataFrame

        Returns
        --------
        None
        """
        model_parameters = model_parameters.copy()
        model_parameters.rename(columns={"experiment": "workflow_id"}, inplace=True)
        model_parameters = self.fix_for_mysql(model_parameters)
        model_parameters["well"] = utils.unpad_well_col(model_parameters["well"])
        self.session.bulk_insert_mappings(
            db_models.NE_model_parameters, model_parameters.to_dict(orient="records")
        )

    def update_workflow_tracking(self, workflow_id: int) -> None:
        """Update workflow_tracking table to indicate all variants for
        a workflow have been uploaded.

        This doesn't check that a workflow is complete, that is handled
        by `DatabaseUploader.is_final_upload()`.

        Parameters
        -----------
        workflow_id: int

        Returns
        -------
        None
        """
        # set status to "complete"
        # set final_results_upload to current datetime
        # set end_date to current datetime
        timestamp = datetime.now(timezone.utc)
        # fmt: off
        self.session\
            .query(db_models.NE_workflow_tracking)\
            .filter(db_models.NE_workflow_tracking.workflow_id == workflow_id)\
            .update(
                {
                    db_models.NE_workflow_tracking.status: "complete",
                    db_models.NE_workflow_tracking.end_date: timestamp,
                    db_models.NE_workflow_tracking.final_results_upload: timestamp,
                }
            )
        # fmt: on

    def upload_reporter_plate_status(self, workflow_id: int, variant: str) -> None:
        """Inserts new row in NE_reporter_plate_status to indicate
        current plate (workflow_id & variant) is awaiting
        reporter decision (pass, fail).

        Parameters
        ----------
        workflow_id: int
        variant: str

        Returns
        -------
        None
        """
        plate_entry = db_models.NE_reporter_plate_status(
            workflow_id=int(workflow_id), variant=variant, status="awaiting"
        )
        self.session.add(plate_entry)

    def upload_titration_results(self, titration_results: pd.DataFrame) -> None:
        """
        Parameters
        ----------
        titration_results: pd.DataFrame

        Returns
        -------
        None
            Uploads results to the LIMS database.
        """
        # can't store NaNs
        titration_results = self.fix_for_mysql(titration_results)
        self.session.bulk_insert_mappings(
            db_models.NE_virus_titration_results,
            titration_results.to_dict(orient="records"),
        )

    def update_titration_workflow_tracking(self, workflow_id: int, variant: str):
        """Update NE_titration_workflow_tracking table to indicate the
        titration has been uploaded.

        Parameters
        ----------
        workflow_id: int
        variant: str
        """
        timestamp = datetime.now(timezone.utc)
        # fmt: off
        self.session\
            .query(db_models.NE_titration_workflow_tracking)\
            .filter(
                db_models.NE_titration_workflow_tracking.workflow_id == workflow_id,
                db_models.NE_titration_workflow_tracking.variant == variant
            )\
            .update(
                {
                    db_models.NE_titration_workflow_tracking.status: "complete",
                    db_models.NE_titration_workflow_tracking.end_date: timestamp
                }
            )
        # fmt: on
