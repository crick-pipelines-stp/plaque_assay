"""
module docstring
"""
import os
import json

import pandas as pd

from . import utils
from .consts import NO_VIRUS_WELLS
from .plate import Plate
from .sample import Sample


class Experiment:
    """Experiment class, holds plates"""

    def __init__(self, df):
        self.df = self.subtract_plaque_area_background(df)
        self.experiment_name = df["Plate_barcode"].values[0][3:]
        self.plate_store = {name: Plate(df) for name, df in df.groupby("Plate_barcode")}
        self.df = pd.concat([plate.df for plate in self.plate_store.values()])
        self.sample_store = self.make_samples()

    @property
    def samples(self):
        """docstring"""
        return self.sample_store.items()

    @property
    def plates(self):
        """docstring"""
        return self.plate_store.items()

    def subtract_plaque_area_background(self, df):
        """
        This is done on an experiment-level rather than plate-level, so belongs here rather
        than the Plate class.
        - Calculate the median of "Normalised Plaque area" fo no virus wells.
        - Subtract median from "Normalised Plaque area" for each well and save
          as "Background Subtracted Plaque Area"
        """
        feature = "Normalised Plaque area"
        new_colname = "Background Subtracted Plaque Area"
        no_virus_bool = df.Well.isin(NO_VIRUS_WELLS)
        background = df[no_virus_bool][feature].median()
        df[new_colname] = df[feature] - background
        return df

    def make_samples(self):
        """
        create a dictionary of {sample_name: Sample}
        """
        sample_dict = dict()
        # NOTE: "Well" might change to be a sample name in the future
        for name, group in self.df.groupby("Well"):
            sample_df = group[["Dilution", "Percentage Infected"]]
            sample_dict[name] = Sample(name, sample_df)
        return sample_dict

    def get_failures_as_json(self):
        """collect all failures into a dictionary ready for JSON output"""
        failure_dict = {}
        failure_dict["plate_failures"] = {}
        failure_dict["well_failures"] = []
        for plate_name, plate_object in self.plates:
            # plate failures
            if plate_object.plate_failed:
                plate_failures_as_dict = [
                    i.to_dict() for i in plate_object.plate_failures
                ]
                failure_dict["plate_failures"][plate_name] = plate_failures_as_dict
            # well failures
            # convert WellFailure objects to dictionaries
            well_failures_as_dict = [i.to_dict() for i in plate_object.well_failures]
            failure_dict["well_failures"].extend(well_failures_as_dict)
        # failures in the Sample object are actually failures of the
        # positive control wells
        for sample_name, sample_object in self.samples:
            sample_failures_as_dict = [i.to_dict() for i in sample_object.failures]
            failure_dict["well_failures"].extend(sample_failures_as_dict)
        return failure_dict

    def get_failures_as_dataframe(self):
        """convert failure JSON file to a dataframe"""
        failures_dict = self.get_failures_as_json()
        types = []
        plates = []
        wells = []
        reasons = []
        # go through plate failures first as have to join multiple wells into a string
        plate_failures = failures_dict["plate_failures"]
        for _, plate_list in plate_failures.items():
            plate = plate_list[0]  # it's a single-element list
            types.append(plate["type"])
            plates.append(plate["plate"])
            wells.append(";".join(plate["wells"]))
            reasons.append(plate["reason"])
        # go through well failures
        well_failures = failures_dict["well_failures"]
        for well_failure in well_failures:
            types.append(well_failure["type"])
            plates.append(well_failure["plate"])
            wells.append(well_failure["well"])
            reasons.append(well_failure["reason"])
        df = pd.DataFrame(
            {
                "failure_type": types,
                "plate": plates,
                "well": wells,
                "failure_reason": reasons,
            }
        )
        df["experiment"] = self.experiment_name
        return df

    def save_failures_as_dataframe(self, output_dir):
        failures_df = self.get_failures_as_dataframe()
        failure_output_path = os.path.join(
            output_dir, f"failures_{self.experiment_name}.csv"
        )
        failures_df.to_csv(failure_output_path, index=False)

    def save_failures_as_json(self, output_dir):
        """docstring"""
        failures = self.get_failures_as_json()
        failure_output_path = os.path.join(
            output_dir, f"failures_{self.experiment_name}.json"
        )
        with open(failure_output_path, "w") as f:
            json.dump(failures, f, indent=4)

    def get_results_as_json(self):
        """collect all results into a dictionary read for JSON output"""
        results = {}
        results["result_mapping"] = utils.INT_TO_RESULT
        results["results"] = {}
        for sample_name, sample in self.samples:
            results["results"][sample_name] = sample.ic50
        return results

    def get_results_as_dataframe(self):
        """convert results JSON file to a dataframe"""
        results_dict = self.get_results_as_json()
        result_mapping = results_dict["result_mapping"]
        wells = []
        ic50s = []
        statuses = []
        for well, value in results_dict["results"].items():
            wells.append(well)
            if value < 0:
                # is an error code
                ic50s.append(None)
                status_string = result_mapping[value]
                statuses.append(status_string)
            else:
                # is an ic50 value
                ic50s.append(value)
                statuses.append(None)
        df = pd.DataFrame({"well": wells, "ic50": ic50s, "status": statuses})
        df["experiment"] = self.experiment_name
        return df

    def save_results_as_dataframe(self, output_dir):
        results_df = self.get_results_as_dataframe()
        result_output_path = os.path.join(
            output_dir, f"results_{self.experiment_name}.csv"
        )
        results_df.to_csv(result_output_path, index=False)

    def save_results_as_json(self, output_dir):
        """docstring"""
        results = self.get_results_as_json()
        result_output_path = os.path.join(
            output_dir, f"results_{self.experiment_name}.json"
        )
        with open(result_output_path, "w") as f:
            json.dump(results, f, indent=4)

    def save_normalised_data(self, output_dir, concatenate=True):
        """docstring"""
        if concatenate:
            # concatenated all dataframes together
            dataframes = []
            for _, plate_object in self.plates:
                df = plate_object.get_normalised_data()
                dataframes.append(df)
            df_concat = pd.concat(dataframes)
            save_path = os.path.join(
                output_dir, f"normalised_{self.experiment_name}.csv"
            )
            df_concat.to_csv(save_path, index=False)
        else:
            # save as individual dataframe per plate
            for _, plate_object in self.plates:
                plate_object.save_normalised_data(output_dir)
