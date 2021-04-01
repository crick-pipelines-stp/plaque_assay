"""
module docstring
"""

import os

import pandas as pd

from . import failure
from . import qc_criteria
from .consts import VIRUS_ONLY_WELLS, NO_VIRUS_WELLS


class Plate:
    """Plate class"""

    def __init__(self, df):
        self.df = self.subtract_plaque_area_background(df)
        assert df["PlateNum"].nunique() == 1
        self.barcode = df["Plate_barcode"].values[0]
        assert df["Dilution"].nunique() == 1
        self.dilution = df["Dilution"].values[0]
        self.plate_failed = False
        self.well_failures = []
        self.plate_failures = []
        self.calc_percentage_infected()
        self.outside_image_area()

    def __len__(self):
        return self.df.shape[0]

    def __str__(self):
        return f"Plate {self.barcode}"

    def outside_image_area(self):
        """docstring"""
        feature = "Cells - Image Region Area [µm²] - Mean per Well"
        experiment_median = self.df[feature].median()
        ratio = self.df[feature] / experiment_median
        self.df["ratio"] = ratio
        lower_limit = qc_criteria.low_cells_image_region_area_low
        upper_limit = qc_criteria.low_cells_image_region_area_high
        low = self.df[self.df["ratio"] < lower_limit]
        high = self.df[self.df["ratio"] > upper_limit]
        outliers = pd.concat([low, high])
        control_outliers = outliers[outliers["Well"].str.endswith("12")]
        if control_outliers.shape[0] > 0:
            # plate failure due to control well failure
            self.plate_failed = True
            self.plate_failures.append(
                failure.CellAreaPlateFailure(
                    plate=self.barcode, wells=control_outliers["Well"].tolist()
                )
            )
        if outliers.shape[0]:
            for _, row in outliers.iterrows():
                self.well_failures.append(
                    failure.WellFailure(
                        well=row["Well"],
                        plate=self.barcode,
                        reason="cell region area outside expected range",
                    )
                )
            # if there's more than 8 failures for DAPI wells, then flag
            # as a possible plate failure
            if outliers.shape[0] > 8:
                # flag possible plate fail
                self.plate_failures.append(
                    failure.DAPIPlateFailure(
                        plate=self.barcode, wells=outliers["Well"].tolist()
                    )
                )

    def subtract_plaque_area_background(self, df):
        """
        This is now done on a plate-by-plate basis.
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

    def calc_percentage_infected(self):
        """docstring"""
        # determine if infection rate of virus-only-wells
        # is within accetable limts, flag the plate if this is false
        lower_limit = qc_criteria.infection_rate_low
        upper_limit = qc_criteria.infection_rate_high
        feature = "Background Subtracted Plaque Area"
        virus_only_bool = self.df.Well.isin(VIRUS_ONLY_WELLS)
        infection = self.df[virus_only_bool][feature].median()
        if infection < lower_limit or infection > upper_limit:
            reason = f"virus-only infection median ({infection:3f}) outside range: ({lower_limit}, {upper_limit})"
            self.plate_failed = True
            self.plate_failures.append(
                failure.InfectionPlateFailure(
                    plate=self.barcode, wells=VIRUS_ONLY_WELLS, reason=reason
                )
            )
        self.df["Percentage Infected"] = self.df[feature] / infection * 100

    def get_normalised_data(self):
        """
        return a simplified dataframe of just the normalised data
        that is saved alongside the final IC50 results
        """
        wanted_cols = [
            "Well",
            "Row",
            "Column",
            "Dilution",
            "Plate_barcode",
            "Background Subtracted Plaque Area",
            "Percentage Infected",
            "variant",
        ]
        temp_df = self.df.copy()
        df_wanted = temp_df[wanted_cols]
        df_wanted = df_wanted.rename(
            columns={
                "Background Subtracted Plaque Area": "Background_subtracted_plaque_area",
                "Percentage Infected": "Percentage_infected",
            }
        )
        return df_wanted

    def save_normalised_data(self, output_dir):
        """
        save csv of the normalised data
        """
        norm_data = self.get_normalised_data()
        save_path = os.path.join(output_dir, f"{self.barcode}.csv")
        norm_data.to_csv(save_path, index=False)
