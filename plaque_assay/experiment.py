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
        self.plate_store = {name: Plate(df) for name, df in df.groupby("PlateNum")}
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

    def collect_failures(self):
        """collect all failures into a dictionary ready for JSON output"""
        failure_dict = {}
        failure_dict["plate_failures"] = {}
        failure_dict["well_failures"] = []
        for plate_name, plate_object in self.plates:
            # plate failures
            if plate_object.plate_failed:
                failure_dict["plate_failures"][plate_name] = plate_object.plate_failures
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

    def save_failures(self, output_dir):
        """docstring"""
        failures = self.collect_failures()
        failure_output_path = os.path.join(
            output_dir, f"failures_{self.experiment_name}.json"
        )
        with open(failure_output_path, "w") as f:
            json.dump(failures, f, indent=4)

    def collect_results(self):
        """collect all results into a dictionary read for JSON output"""
        results = {}
        results["result_mapping"] = utils.INT_TO_RESULT
        results["results"] = {}
        for sample_name, sample in self.samples:
            results["results"][sample_name] = sample.ic50
        return results

    def save_results(self, output_dir):
        """docstring"""
        results = self.collect_results()
        result_output_path = os.path.join(
            output_dir, f"results_{self.experiment_name}.json"
        )
        with open(result_output_path, "w") as f:
            json.dump(results, f, indent=4)

    def save_normalised_data(self, output_dir):
        """docstring"""
        for _, plate_object in self.plates:
            plate_object.save_normalised_data(output_dir)
