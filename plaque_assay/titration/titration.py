"""
Titration class is similar to the Experiment class but used
with titration plates.
"""

from typing import Dict, List, Optional

import pandas as pd

from .sample import Sample
from .titration_dilution import TitrationDilution
from . import utils


class Titration:
    """Titration class, holds dilution samples

    Parameters
    -----------
    titration_dataset: pd.DataFrame
    """

    def __init__(self, titration_dataset: pd.DataFrame, variant: str):
        self.dataset = titration_dataset
        self.variant = variant
        self.workflow_id = self.dataset["Plate_barcode"].values[0][3:]
        dilution_store = dict()
        for dilution, df in titration_dataset.groupby("Virus_dilution_factor"):
            dilution_store[dilution] = TitrationDilution(df)
        self.dilution_store = dilution_store
        self.df = pd.concat([dilution.df for dilution in self.dilution_store.values()])
        self.sample_store = self.make_samples()

    @property
    def samples(self):
        """key-value store of all samples in the experiment"""
        return self.sample_store.items()

    def make_samples(self) -> Dict[str, Sample]:
        """
        Returns samples in the form of a dictionary

        Returns
        --------
        dict
            `{sample_name: Sample`}
        """
        sample_dict: Dict[str, Sample] = dict()
        for dilution, group in self.df.groupby("Virus_dilution_factor"):
            sample_df = group[["Dilution", "Percentage Infected"]]
            sample_dict[dilution] = Sample(str(dilution), sample_df, self.variant)
        return sample_dict

    def get_titration_results(self) -> pd.DataFrame:
        """
        return dataframe of titration results suitable for uploading to
        the LIMS database. e.g:
        +----------+------+--------+--------------------+-------------------------------+---------+-------------+
        | dilution | ic50 | status | mean_squared_error | median_virus_only_plaque_area | variant | workflow_id |
        +----------+------+--------+--------------------+-------------------------------+---------+-------------+
        |          |      |        |                    |                               |         |             |

        Arguments
        ---------
        None

        Returns
        -------
        pd.DataFrame
        """
        dilutions: List[int] = []
        ic50s: List[Optional[float]] = []
        statuses: List[Optional[str]] = []
        mean_squared_errors: List[Optional[float]] = []
        infection_rates: List[float] = []
        for dilution, dilution_sample in self.samples:
            dilutions.append(dilution)
            titration_dilution_object = self.dilution_store[dilution]
            infection_rate = titration_dilution_object.infection_rate
            infection_rates.append(infection_rate)
            mean_squared_errors.append(dilution_sample.mean_squared_error)
            if dilution_sample.ic50 < 0:
                # is an error code
                ic50s.append(None)
                status_str = utils.int_to_result(dilution_sample.ic50)
                statuses.append(status_str)
            else:
                # is a valid ic50 value
                ic50s.append(dilution_sample.ic50)
                statuses.append(None)
        results_df = pd.DataFrame(
            {
                "dilution": dilutions,
                "ic50": ic50s,
                "status": statuses,
                "mean_squared_error": mean_squared_errors,
                "median_virus_only_plaque_area": infection_rates,
            }
        )
        results_df["variant"] = self.variant
        results_df["workflow_id"] = self.workflow_id
        return results_df
