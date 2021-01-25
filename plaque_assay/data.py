import logging
import os
import time

import pandas as pd
import numpy as np

from . import utils
from . import consts
from . import db_models


def read_data_from_list(plate_list):
    plate_name_dict = {
        os.path.abspath(i): utils.get_dilution_from_barcode(i) for i in plate_list
    }
    dataframes = []
    for path, plate_num in plate_name_dict.items():
        df = pd.read_csv(
            # NOTE: might not always be Evaluation1
            os.path.join(path, "Evaluation1/PlateResults.txt"),
            skiprows=8,
            sep="\t",
        )
        plate_barcode = path.split(os.sep)[-1].split("__")[0]
        logging.info("plate barcode detected as %s", plate_barcode)
        df["Dilution"] = consts.plate_mapping[plate_num]
        # add well labels to dataframe
        well_labels = []
        for row, col in df[["Row", "Column"]].itertuples(index=False):
            well_labels.append(utils.row_col_to_well(row, col))
        df["Well"] = well_labels
        df["PlateNum"] = plate_num
        df["Plate_barcode"] = plate_barcode
        dataframes.append(df)
    df_concat = pd.concat(dataframes)
    logging.debug("input data shape: %s", df_concat.shape)
    return df_concat


def get_plate_list(data_dir):
    plate_list = [os.path.join(data_dir, i) for i in os.listdir(data_dir)]
    if len(plate_list) == 8:
        logging.debug("plate list detected: %s", plate_list)
    else:
        logging.error(
            "Did not detect 8 plates, detected %s :", len(plate_list), plate_list
        )
    return plate_list


def read_data_from_directory(data_dir):
    """
    read actual barcoded plate directory
    """
    plate_list = get_plate_list(data_dir)
    return read_data_from_list(plate_list)


def read_indexfiles_from_list(plate_list):
    plate_name_dict = {
        os.path.abspath(i): utils.get_dilution_from_barcode(i) for i in plate_list
    }
    dataframes = []
    for path, plate_num in plate_name_dict.items():
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


def upload_plate_results(session, plate_results_dataset):
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
        # "Background Subtracted Plaque Area": "background_subtracted_plaque_area",
    }
    plate_results_dataset.rename(columns=rename_dict, inplace=True)
    # filter to only desired columns
    plate_results_dataset = plate_results_dataset[list(rename_dict.values())]
    workflow_id = [int(i[3:]) for i in plate_results_dataset["plate_barcode"]]
    plate_results_dataset["workflow_id"] = workflow_id
    plate_results_dataset["well"] = utils.unpad_well_col(plate_results_dataset["well"])
    session.bulk_insert_mappings(
        db_models.NE_raw_results, plate_results_dataset.to_dict(orient="records")
    )
    # FIXME: update status and result_upload_time for NE_workflow_tracking
    session.commit()


def upload_indexfiles(session, indexfiles_dataset):
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
    }
    indexfiles_dataset.rename(columns=rename_dict, inplace=True)
    # filter to only desired columns
    indexfiles_dataset = indexfiles_dataset[list(rename_dict.values())]
    # get workflow ID
    workflow_id = [int(i[3:]) for i in indexfiles_dataset["plate_barcode"]]
    indexfiles_dataset["workflow_id"] = workflow_id
    indexfiles_dataset["well"] = utils.unpad_well_col(indexfiles_dataset["well"])
    for i in range(0, len(indexfiles_dataset), 1000):
        df_slice = indexfiles_dataset.iloc[i : i + 1000]
        session.bulk_insert_mappings(
            db_models.NE_raw_index, df_slice.to_dict(orient="records")
        )
    # FIXME: update status and result_upload_time for NE_workflow_tracking
    session.commit()


def upload_normalised_results(session, norm_results):
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
    }
    norm_results.rename(columns=rename_dict, inplace=True)
    norm_results = norm_results[list(rename_dict.values())]
    workflow_id = [int(i[3:]) for i in norm_results["plate_barcode"]]
    assert len(set(workflow_id)) == 1
    norm_results["workflow_id"] = workflow_id
    norm_results["well"] = utils.unpad_well_col(norm_results["well"])
    session.bulk_insert_mappings(
        db_models.NE_normalized_results, norm_results.to_dict(orient="records")
    )
    # FIXME: update status and result_upload time for NE_workflow_tracking
    session.commit()


def upload_final_results(session, results):
    """docstring"""
    results = results.copy()
    # TODO: double-check what master_plate is??
    results["master_plate"] = None
    # get workflow_id
    assert results["experiment"].nunique() == 1
    results["workflow_id"] = results["experiment"].astype(int)
    # can't store NaN in mysql, to convert to None which are stored as null
    results = results.replace({np.nan: None})
    results["well"] = utils.unpad_well_col(results["well"])
    session.bulk_insert_mappings(
        db_models.NE_final_results, results.to_dict(orient="records")
    )
    # FIXME: update status and results upload time for NE_workflow tracking
    session.commit()


def upload_failures(session, failures):
    """docsring"""
    failures = failures.copy()
    # FIXME: get workflow_id
    assert failures["experiment"].nunique() == 1
    failures["workflow_id"] = failures["experiment"].astype(int)
    failures["well"] = utils.unpad_well_col(failures["well"])
    session.bulk_insert_mappings(
        db_models.NE_failed_results, failures.to_dict(orient="records")
    )
    session.commit()


def upload_model_parameters(session, model_parameters):
    """docstring"""
    model_parameters = model_parameters.copy()
    model_parameters.rename(columns={"experiment": "workflow_id"}, inplace=True)
    # can't store NaNs
    model_parameters = model_parameters.replace({np.nan: None})
    model_parameters["well"] = utils.unpad_well_col(model_parameters["well"])
    session.bulk_insert_mappings(
        db_models.NE_model_parameters, model_parameters.to_dict(orient="records")
    )
    session.commit()
