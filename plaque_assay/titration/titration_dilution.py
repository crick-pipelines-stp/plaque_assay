"""
similar to the Plate class for the analysis, but does not
have QC checking.
"""

import pandas as pd

from . import consts


class TitrationDilution:
    """
    A TitrationDilution covers 2 columns relating to a single virus dilution
    factor on a titration plate.

    Parameters
    -----------
    df: pd.DataFrame
    """

    def __init__(self, df: pd.DataFrame):
        assert (
            df["Virus_dilution_factor"].nunique() == 1
        ), "more than 1 dilution factor present"
        self.df = df
        self.df_background_subtracted = self.subtract_plaque_area_background(df)
        self.dilution = df["Dilution"].values[0]
        self.calc_percentage_infected()

    @property
    def infection_rate(self) -> float:
        """alias for self.median_virus_only_plaque_area"""
        return self.median_virus_only_plaque_area

    @property
    def median_virus_only_plaque_area(self) -> float:
        """
        median of the virus only well's plaque area is equivalent
        to the infection rate.
        """
        feature = "Background Subtracted Plaque Area"
        virus_only_bool = self.df_background_subtracted.Well.isin(
            consts.TITRATION_VIRUS_ONLY_WELLS
        )
        return self.df_background_subtracted[virus_only_bool][feature].median()

    def subtract_plaque_area_background(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove background from `plaque_area`.

        1. Calculate the median of "Normalised Plaque area" fo no virus wells.

        2. Subtract median from "Normalised Plaque area" for each well and save
          as "Background Subtracted Plaque Area"

        This is done individually for each pair of dilution columns.

        Parameters
        -----------
        df : pandas.DataFrame

        Returns
        -------
        pandas.DataFrame
        """
        feature = "Normalised Plaque area"
        new_colname = "Background Subtracted Plaque Area"
        # self.df is a subset of the plate only containing 2 columns for a
        # dilution
        no_virus_bool = df.Well.isin(consts.TITRATION_NO_VIRUS_WELLS)
        background = df[no_virus_bool][feature].median()
        df[new_colname] = df[feature] - background
        return df

    def calc_percentage_infected(self) -> None:
        """Calculate percentage infected. Adds "Percentage Infected" column
        to self.df.

        `Percentage Infected` is the `Background Subtracted Plaque Area`
        divided by the median of the virus only wells x 100.
        """
        feature = "Background Subtracted Plaque Area"
        infection_rate = self.median_virus_only_plaque_area
        self.df["Percentage Infected"] = (
            self.df_background_subtracted[feature] / infection_rate * 100
        )
