"""
Titration class is similar to the Experiment class but used
with titration plates.
"""

from typing import Dict

import pandas as pd

from .sample import Sample


class Titration:
    """Titration class, holds dilution samples

    Parameters
    -----------
    titration_dataset: pd.DataFrame
    """

    def __init__(self, titration_dataset: pd.DataFrame):
        self.dataset = titration_dataset
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
        sample_dict = dict()
        raise NotImplementedError()
        # TODO: loop through and generate a dataframe for each sample with two
        # columns ["Dilution", "Percentage Infected"]
        return sample_dict

    def get_titration_results(self) -> pd.DataFrame:
        """
        return dataframe of titration results suitable for uploading to
        the LIMS database.

        Arguments
        ---------
        None

        Returns
        -------
        pd.DataFrame
        """
        raise NotImplementedError()

    @property
    def workflow_id(self):
        """get workflow_id from self.dataset"""
        raise NotImplementedError()
        pass

    @property
    def variant(self):
        """get variant name from self.dataset"""
        raise NotImplementedError()
        pass
