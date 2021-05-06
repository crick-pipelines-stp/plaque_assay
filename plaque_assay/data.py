from datetime import datetime, timezone
import logging
import os

import pandas as pd
import numpy as np

from . import utils
from . import consts
from . import db_models


def read_data_from_list(plate_list):
    """
    Read in data from plate list and assign dilution values by well position.
    NOTE: this will mock the data so the 4 dilutions on a single 384-well plate
          are re-labelled to appear from 4 different 96 well plates.
    """
    dataframes = []
    barcodes = []
    for path in plate_list:
        df = pd.read_csv(
            # NOTE: might not always be Evaluation1
            os.path.join(path, "Evaluation1/PlateResults.txt"),
            skiprows=8,
            sep="\t",
        )
        plate_barcode = path.split(os.sep)[-1].split("__")[0]
        barcodes.append(plate_barcode)
        logging.info("plate barcode detected as %s", plate_barcode)
        well_labels = []
        for row, col in df[["Row", "Column"]].itertuples(index=False):
            well_labels.append(utils.row_col_to_well(row, col))
        df["Well"] = well_labels
        df["Plate_barcode"] = plate_barcode
        dataframes.append(df)
    df_concat = pd.concat(dataframes)
    # NOTE: mock barcodes before changing wells
    df_concat["Plate_barcode"] = utils.mock_384_barcode(
        existing_barcodes=df_concat["Plate_barcode"], wells=df_concat["Well"]
    )
    # mock wells
    df_concat["Well"] = [utils.well_384_to_96(i) for i in df_concat["Well"]]
    df_concat["PlateNum"] = [int(i[1]) for i in df_concat["Plate_barcode"]]
    df_concat["Dilution"] = [consts.plate_mapping[i] for i in df_concat["PlateNum"]]
    logging.debug("input data shape: %s", df_concat.shape)
    return df_concat


def get_plate_list(data_dir):
    n_expected_plates = 2
    plate_list = [os.path.join(data_dir, i) for i in os.listdir(data_dir)]
    if len(plate_list) == n_expected_plates:
        logging.debug("plate list detected: %s", plate_list)
    else:
        logging.error(
            "Did not detect %s plates, detected %s :",
            n_expected_plates,
            len(plate_list),
        )
    return plate_list


def read_data_from_directory(data_dir):
    """
    read actual barcoded plate directory
    """
    plate_list = get_plate_list(data_dir)
    return read_data_from_list(plate_list)


def read_indexfiles_from_list(plate_list):
    dataframes = []
    for path in plate_list:
        df = pd.read_csv(os.path.join(path, "indexfile.txt"), sep="\t")
        plate_barcode = path.split(os.sep)[-1].split("__")[0]
        df["Plate_barcode"] = plate_barcode
        dataframes.append(df)
    df_concat = pd.concat(dataframes)
    # remove annoying empty "Unnamed: 16" column
    to_rm = [col for col in df_concat.columns if col.startswith("Unnamed:")]
    df_concat.drop(to_rm, axis=1, inplace=True)
    if len(to_rm) > 0:
        logging.info("removed columns: %s from concatenated indexfiles", to_rm)
    logging.debug("indexfile shape: %s", df_concat.shape)
    return df_concat


def read_indexfiles_from_directory(data_dir):
    plate_list = get_plate_list(data_dir)
    return read_indexfiles_from_list(plate_list)


def get_awaiting_raw(session, master_plate):
    """
    docstring
    """
    # FIXME: this table doesn't actually exist
    awaiting_raw = (
        session.query(db_models.NE_workflow_tracking)
        .filter(db_models.NE_workflow_tracking.master_plate == master_plate)
        .filter(db_models.NE_workflow_tracking.status == "raw_results")
        .first()
    )
    return awaiting_raw


class DatabaseUploader:
    def __init__(self, session):
        self.session = session

    def commit(self):
        self.session.commit()

    def already_uploaded(self, workflow_id, variant):
        """
        Check if the results for a given workflow_id and variant
        have already been uploaded.

        As results are uploaded all-or-nothing, we can just check one of the
        results tables for the presence of the supplied workflow_id and
        variant, rather than each of the several results tables.

        Parameters:
        -----------
        workflow_id: int
        variant: string

        Returns:
        --------
        bool
        """
        result = (
            self.session.query(db_models.NE_final_results)
            .filter(
                db_models.NE_final_results.workflow_id == workflow_id,
                db_models.NE_final_results.variant == variant,
            )
            .first()
        )
        return result is not None

    def is_final_upload(self, workflow_id):
        """
        This determines if a given workflow_id and variant are the last
        to be uploaded for that workflow_id.

        More description here.

        Parameters:
        -----------
        workflow_id: int

        Returns:
        ---------
        bool
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

    def upload_plate_results(self, plate_results_dataset):
        """upload raw concatenated data into database"""
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
        # can't store NaN in mysql, to convert to None which are stored as null
        plate_results_dataset = plate_results_dataset.replace({np.nan: None})
        self.session.bulk_insert_mappings(
            db_models.NE_raw_results, plate_results_dataset.to_dict(orient="records")
        )

    def upload_indexfiles(self, indexfiles_dataset):
        """docstring"""
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
        # can't store NaN in mysql, to convert to None which are stored as null
        indexfiles_dataset = indexfiles_dataset.replace({np.nan: None})
        for i in range(0, len(indexfiles_dataset), 1000):
            df_slice = indexfiles_dataset.iloc[i : i + 1000]
            self.session.bulk_insert_mappings(
                db_models.NE_raw_index, df_slice.to_dict(orient="records")
            )

    def upload_normalised_results(self, norm_results):
        """docstring"""
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
        # can't store NaN in mysql, to convert to None which are stored as null
        norm_results = norm_results.replace({np.nan: None})
        self.session.bulk_insert_mappings(
            db_models.NE_normalized_results, norm_results.to_dict(orient="records")
        )

    def upload_final_results(self, results):
        """docstring"""
        results = results.copy()
        # TODO: double-check what master_plate is??
        results["master_plate"] = None
        # get workflow_id
        assert results["experiment"].nunique() == 1
        assert results["variant"].nunique() == 1
        results["workflow_id"] = results["experiment"].astype(int)
        # can't store NaN in mysql, to convert to None which are stored as null
        results = results.replace({np.nan: None})
        results["well"] = utils.unpad_well_col(results["well"])
        self.session.bulk_insert_mappings(
            db_models.NE_final_results, results.to_dict(orient="records")
        )

    def upload_failures(self, failures):
        """docsring"""
        failures = failures.copy()
        # FIXME: get workflow_id
        if failures.shape[0] > 0:
            assert failures["experiment"].nunique() == 1
            failures["workflow_id"] = failures["experiment"].astype(int)
            # FIXME: doesn't work with multiple well in plate failures
            # failures["well"] = utils.unpad_well_col(failures["well"])
            self.session.bulk_insert_mappings(
                db_models.NE_failed_results, failures.to_dict(orient="records")
            )

    def upload_model_parameters(self, model_parameters):
        """docstring"""
        model_parameters = model_parameters.copy()
        model_parameters.rename(columns={"experiment": "workflow_id"}, inplace=True)
        # can't store NaNs
        model_parameters = model_parameters.replace({np.nan: None})
        model_parameters["well"] = utils.unpad_well_col(model_parameters["well"])
        self.session.bulk_insert_mappings(
            db_models.NE_model_parameters, model_parameters.to_dict(orient="records")
        )

    def upload_barcode_changes_384(self, workflow_id):
        """docstring"""
        replicates = [1, 2]
        # NOTE: dilution numbers match the LIMS values rather than assay values
        assay_plates = []
        ap_10 = []
        ap_40 = []
        ap_160 = []
        ap_640 = []
        workflow_ids = []
        for replicate in replicates:
            assay_plate = f"AA{replicate}{workflow_id}"
            assay_plates.append(assay_plate)
            ap_10.append(f"A1{replicate}{workflow_id}")
            ap_40.append(f"A2{replicate}{workflow_id}")
            ap_160.append(f"A3{replicate}{workflow_id}")
            ap_640.append(f"A4{replicate}{workflow_id}")
            workflow_ids.append(workflow_id)
        df = pd.DataFrame(
            {
                "assay_plate_384": assay_plates,
                "ap_10": ap_10,
                "ap_40": ap_40,
                "ap_160": ap_160,
                "ap_640": ap_640,
                "workflow_id": workflow_ids,
            }
        )
        self.session.bulk_insert_mappings(
            db_models.NE_assay_plate_tracker_384, df.to_dict(orient="records")
        )

    def is_awaiting_results(self, workflow_id):
        """
        docstring
        """
        raise NotImplementedError()

    def update_workflow_tracking(self, workflow_id):
        """
        Update the status in NE_workflow_tracking indicating
        the analysis for a workflow is complete.
        """
        # set status to "complete"
        # set final_results_upload to current datetime
        # set end_date to current datetime
        timestamp = datetime.now(timezone.utc)
        # fmt: off
        self.session\
            .query(db_models.NE_workflow_tracking)\
            .filter( db_models.NE_workflow_tracking.workflow_id == workflow_id)\
            .update(
                {
                    db_models.NE_workflow_tracking.status: "complete",
                    db_models.NE_workflow_tracking.end_date: timestamp,
                    db_models.NE_workflow_tracking.final_results_upload: timestamp,
                }
            )
        # fmt: on
