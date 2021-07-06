"""
Mock 96-well plates of a single dilution, taken from the physical
384 well plates.

The assay originated in 96-well plates, with a dilution per plate
totalling 4 plates.  This was then taken forward when the assay moved
to 384-well plates, where all 4 dilutions were present on a single
384-well plate, with samples dilution in quadrants, so a sample was
deposited in 4 adjacent wells, each well a different dilution.

e.g in the physical 384-well plate where `A` and `B` are the samples,
and `1`, `2`, `3`, `4` are the dilution integers.

```txt
    01   02
  +----+----+
A | A4 | A3 |
  +----+----+
B | A2 | A1 |
  +----+----+
C | B4 | B3 |
  +----+----+
D | B2 | B1 |
  +----+----+
```

This also means that the well labels in this class have been converted from their
384-well plate quadrant to the 96 well equivalent. This is carried out in
`plaque_assay.utils.well_384_to_96()`.
"""

import os
from typing import Set

import pandas as pd

from . import failure
from . import qc_criteria
from .consts import VIRUS_ONLY_WELLS, NO_VIRUS_WELLS


class Plate:
    """Plate class

    This is not the physical 384-well plate, but the mock 96-well
    plate of a single dilution.


    Parameters
    -----------
    df : pandas.DataFrame

    Attributes
    -----------
    df : pandas.DataFrame
        Dataframe, similar to input `df` but with additional columns from
        `Plate.subtract_plaque_area_background()` and
        `Plate.calc_percentage_infected()`.
    barcode : string
        Mock plate barcode. This is not the scanned barcode,
        but the mock 96-well plate barcode determined from the
        actual 384-well barcode.
    dilution : int
        Dilution integer of the plate, (1, 2, 3, 4).
    plate_failed : bool
        `True` if the entire plate is a QC failure, otherwise `False`.
    well_failures : list
        List containing WellFailure classes if any wells have failed.
    plate_failures : list
        List containing PlateFailure classes if the plate has failed
        for one of more QC reasons.
    """

    def __init__(self, df: pd.DataFrame):
        self.df = self.subtract_plaque_area_background(df)
        assert df["PlateNum"].nunique() == 1
        self.barcode = df["Plate_barcode"].values[0]
        assert df["Dilution"].nunique() == 1
        self.dilution = df["Dilution"].values[0]
        self.variant = df["variant"].values[0]
        self.plate_failed = False
        self.well_failures: Set[failure.WellFailure] = set()
        self.plate_failures: Set[failure.PlateFailure] = set()
        self.calc_percentage_infected()
        self.outside_image_area()

    def __len__(self):
        return self.df.shape[0]

    def __str__(self):
        return f"Plate {self.barcode}"

    def outside_image_area(self) -> None:
        """QC check for `cell_region_area`

        Determines if `cell_region_area` is outside the expected range and
        adds any failures to `well_failures`.
        """
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
            failed_plate = failure.PlateFailure(
                plate=self.barcode,
                well=";".join(control_outliers["Well"]),
                failure_reason=failure.CELL_IMAGE_AREA_FAILURE_REASON,
            )
            self.plate_failures.add(failed_plate)
        if outliers.shape[0]:
            for _, row in outliers.iterrows():
                well_failure = failure.WellFailure(
                    well=row["Well"],
                    plate=self.barcode,
                    failure_reason=failure.CELL_REGION_FAILURE_REASON,
                )
                self.well_failures.add(well_failure)
            # if there's more than 8 failures for DAPI wells, then flag
            # as a possible plate failure
            if outliers.shape[0] > 8:
                # flag possible plate fail
                self.plate_failed = True
                failed_plate = failure.PlateFailure(
                    plate=self.barcode,
                    well=";".join(outliers["Well"]),
                    failure_reason=failure.DAPI_PLATE_FAILURE_REASON,
                )
                self.plate_failures.add(failed_plate)

    def subtract_plaque_area_background(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove background from `plaque_area`.

        1. Calculate the median of "Normalised Plaque area" fo no virus wells.

        2. Subtract median from "Normalised Plaque area" for each well and save
          as "Background Subtracted Plaque Area"

        This is now done on a plate-by-plate basis.

        Parameters
        -----------
        df : pandas.DataFrame

        Returns
        -------
        pandas.DataFrame
        """
        feature = "Normalised Plaque area"
        new_colname = "Background Subtracted Plaque Area"
        no_virus_bool = df.Well.isin(NO_VIRUS_WELLS)
        background = df[no_virus_bool][feature].median()
        df[new_colname] = df[feature] - background
        return df

    def check_infection(self, infection: float) -> None:
        """
        Determines if infection rate of virus-only-wells is within
        acceptable limts, flag as a plate failure if this is false.

        Parameters
        -----------
        infection: float
            infection rate of virus-only wells

        Returns
        -------
        None
            Adds a `plaque_assay.failure.PlateFailure` to
            `self.plate_failures` if plate has failed.
        """
        infection_limits = qc_criteria.infection_rate[self.variant]
        lower_limit = infection_limits["low"]
        upper_limit = infection_limits["high"]
        if infection < lower_limit or infection > upper_limit:
            reason = f"virus-only infection median ({infection:3f}) outside range: ({lower_limit}, {upper_limit})"
            self.plate_failed = True
            failed_plate = failure.PlateFailure(
                plate=self.barcode,
                well=";".join(VIRUS_ONLY_WELLS),
                failure_reason=reason,
            )
            self.plate_failures.add(failed_plate)

    def calc_percentage_infected(self) -> None:
        """Calculate percentage infected.

        `Percentage Infected` is the `Background Subtracted Plaque Area`
        divided by the median of the virus only wells x 100.

        Notes
        ------
        Also carries out a QC check. Determines if infection rate of
        virus-only-wells is within acceptable limts, flag the plate if
        this is false.
        """
        feature = "Background Subtracted Plaque Area"
        virus_only_bool = self.df.Well.isin(VIRUS_ONLY_WELLS)
        infection = self.df[virus_only_bool][feature].median()
        self.check_infection(infection)
        self.df["Percentage Infected"] = self.df[feature] / infection * 100

    def get_normalised_data(self) -> pd.DataFrame:
        """Return a simplified dataframe of just the normalised data

        This subsets and renames some columns.

        Returns
        -------
        pandas.DataFrame
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

    def save_normalised_data(self, output_dir: str) -> None:
        """Save csv of the normalised data

        Parameters
        -----------
        output_dir : str
            directory in which to save the csv file
        """
        norm_data = self.get_normalised_data()
        save_path = os.path.join(output_dir, f"{self.barcode}.csv")
        norm_data.to_csv(save_path, index=False)
