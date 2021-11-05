"""
Titration class is similar to the Experiment class but used
with titration plates.
"""

from collections import defaultdict
from typing import Dict, List, Optional

import pandas as pd

from plaque_assay.sample import Sample
from plaque_assay.titration_dilution import TitrationDilution
from plaque_assay import utils


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
        for (dilution, nanobody), group in self.df.groupby(
            ["Virus_dilution_factor", "nanobody"]
        ):
            sample_df = group[["Dilution", "Percentage Infected"]]
            sample_name = f"{dilution}-{int(nanobody)}"
            sample_dict[sample_name] = Sample(sample_name, sample_df, self.variant)
        return sample_dict

    def get_final_results(self) -> pd.DataFrame:
        """
        return dataframe of titration results suitable for uploading to
        the LIMS database. e.g:
        ```
        +----------+----------+------+--------+-------------+
        | dilution | nanobody | ic50 | status | workflow_id |
        +----------+----------+------+--------+-------------+
        |          |          |      |        |             |
        ```

        Arguments
        ---------
        None

        Returns
        -------
        pd.DataFrame
        """
        dilutions: List[int] = []
        nanobodies: List[int] = []
        ic50s: List[Optional[float]] = []
        statuses: List[Optional[str]] = []
        for sample_name, dilution_sample in self.samples:
            dilution, nanobody = sample_name.split("-")
            dilution = int(dilution)
            nanobody = int(nanobody)
            dilutions.append(dilution)
            nanobodies.append(nanobody)
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
                "nanobody": nanobodies,
                "ic50": ic50s,
                "status": statuses,
            }
        )
        results_df["workflow_id"] = self.workflow_id
        return results_df

    def get_model_parameters(self) -> pd.DataFrame:
        """
        return a dataframe of model parameters

        Arguments
        ----------
        None

        Returns
        --------
        pd.DataFrame
        """
        param_dict = defaultdict(list)
        for sample_name, dilution_sample in self.samples:
            dilution, nanobody = sample_name.split("-")
            dilution = int(dilution)
            nanobody = int(nanobody)
            model_params = dilution_sample.model_params
            mean_squared_error = dilution_sample.mean_squared_error
            if model_params:
                top, bottom, ec50, hillslope = model_params
            else:
                top, bottom, ec50, hillslope = None, None, None, None
            param_dict["dilution"].append(dilution)
            param_dict["nanobody"].append(nanobody)
            param_dict["param_top"].append(top)
            param_dict["param_bottom"].append(bottom)
            param_dict["param_ec50"].append(ec50)
            param_dict["param_hillslope"].append(hillslope)
            param_dict["mean_squared_error"].append(mean_squared_error)
        df = pd.DataFrame(param_dict)
        df["workflow_id"] = self.workflow_id
        return df

    def get_normalised_results(self) -> pd.DataFrame:
        """
        return a dataframe of normalised results

        Arguments
        ----------
        None

        Returns
        --------
        pd.DataFrame
        """
        dilution_dataframe_list = []
        for dilution_name, dilution_sample in self.dilution_store.items():
            tmp_df = dilution_sample.df
            tmp_df["dilution"] = dilution_name
            dilution_dataframe_list.append(tmp_df)
        df_all = pd.concat(dilution_dataframe_list)
        df_all["workflow_id"] = self.workflow_id
        return df_all
