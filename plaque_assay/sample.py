"""
module docstring
"""
from typing import List, Union, Set

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plaque_assay import stats
from plaque_assay import failure
from plaque_assay import utils
from plaque_assay import qc_criteria
from plaque_assay.consts import POSITIVE_CONTROL_WELLS


class Sample:
    """Sample object holds the data for a sample across 4 concentrations
    including replicates.

    The key attribute here is a dataframe of 2 columns: `dilution` and `value`.
    Which is used to calculate the IC50 and detect and failures for the sample.
    `value` is the percentage infected metric.

    Parameters
    -----------
    sample_name : string
        sample name, typically well-label
    data : pandas.DataFrame
        2 column dataframe: [dilution, value]
    variant : string
        virus variant name, used for variant-specific QC checks

    Attributes
    ----------
    sample_name : str
        Typically a well label e.g "A01".
    data : pd.DataFrame
        Dataframe containing results filtered to a single sample. Containing
        2 columns: `dilution` and `value`.
    failures : set
        set containing potential `failures.WellFailures`. These include reasons
        such as:
        - invalid IC50 value if the sample is a positive control
        - failing to fit a model
        - model with mean squared error > 150
        - discordant replicate values
    is_positive_control : bool
        Whether or not the sample is a positive control or not.
    ic50 : float
        Calculated ic50 values.
    ic50_pretty : float or str
        Calculated IC50 value, but negative integer values that indicate
        categories such as "complete inhbition" are converted to their
        string representation.
    model_params : array
        Fitted model parameters from `scipy.optimize.curve_fit()` for
        the 4-parameter dose-response model. These are
        `(top, bottom, ic50, hillslope)`
    fit_method : str
        How the model was fitted, either through an heuristic or
        via fitting a curve.
    """

    def __init__(self, sample_name: str, data: pd.DataFrame, variant: str):
        self.sample_name = sample_name
        self.data = data
        self.variant = variant
        self.failures: Set[Union[failure.WellFailure, failure.PlateFailure]] = set()
        self.calc_ic50()
        self.is_positive_control = sample_name in POSITIVE_CONTROL_WELLS
        self.check_positive_control()
        self.check_duplicate_differences()
        self.check_for_model_fit_failure()

    def calc_ic50(self):
        """calculate IC50 value

        This calculates an IC50 value if possible, though this could also
        be an heuristic such as `complete inhibition`, `no inhibition`,
        `weak inhibition` or `model fit failure`.

        The `ic50` attribute is always numeric, if it is in heuristic it is
        stored as a negative integer (see `utils.int_to_result`).
        the `ic50_pretty` attribute is used for display, where it will
        show either the numeric ic50 value, or the heuristic text.

        This method also calculates the `mean_squared_error` of the model
        if a model has been used to calculate the IC50.

        Parameters
        ----------

        Returns
        --------
        None
        """
        model_results = stats.calc_model_results(self.sample_name, self.data)
        self.fit_method = model_results.fit_method
        self.ic50 = model_results.result
        self.ic50_pretty = (
            self.ic50 if self.ic50 > 0 else utils.int_to_result(int(self.ic50))
        )
        self.model_params = model_results.model_params
        self.mean_squared_error = model_results.mean_squared_error
        self.check_for_model_mse_failure(qc_criteria.mse_upper_limit)

    def check_positive_control(self) -> None:
        """
        If this sample is a postitive control, then determine
        if the IC50 value is between a variant-specific range.
        If it's outside this range it's flagged as a well failure.

        Appends `failure.WellFailure`s to the `failures` attribute.

        Parameters
        ----------
        None

        Returns
        --------
        None
        """
        if not self.is_positive_control:
            return None
        lower_limit = qc_criteria.positive_control_ic50[self.variant]["low"]
        upper_limit = qc_criteria.positive_control_ic50[self.variant]["high"]
        if self.ic50 < lower_limit or self.ic50 > upper_limit:
            reason = f"positive control failure. IC50 = {self.ic50_pretty} not in range ({lower_limit}, {upper_limit})"
            positive_control_failure = failure.WellFailure(
                plate="DILUTION SERIES", well=self.sample_name, failure_reason=reason,
            )
            self.failures.add(positive_control_failure)

    def check_duplicate_differences(self) -> None:
        """
        Want to calculate the differences between duplicates and determine
        if the difference is larger than some threshold. If this is the case
        for more than 2 duplicates in a sample then this should be flagged as
        a QC well failure.

        Appends `failure.WellFailure`s to the `failures` attribute.

        Parameters
        ----------

        Returns
        --------
        None
        """
        failed_count = 0
        difference_threshold = qc_criteria.duplicate_difference
        if self.ic50_pretty in ("no inhibition", "failed to fit model"):
            # don't flag for bad replicates if there's no inhbition
            # or aleady failed model fit
            return None
        for _, group in self.data.groupby("Dilution"):
            if group.shape[0] == 2:
                x, y = group["Percentage Infected"].values
                difference = abs(x - y)
                if difference >= difference_threshold:
                    failed_count += 1
        if failed_count >= 2:
            # is a well failure
            duplicate_failure = failure.WellFailure(
                plate="DILUTION SERIES",
                well=self.sample_name,
                failure_reason=f"2 or more duplicates differ by >= {difference_threshold} % infected",
            )
            self.failures.add(duplicate_failure)

    def check_for_model_fit_failure(self) -> None:
        """
        If the model has failed to fit this should also be flagged as a QC failure

        Appends `failure.WellFailure`s to the `failures` attribute.

        Parameters
        ----------

        Returns
        --------
        None
        """
        if self.ic50_pretty == "failed to fit model" or self.ic50 == -999:
            model_fit_failure = failure.WellFailure(
                plate="DILUTION SERIES",
                well=self.sample_name,
                failure_reason="failed to fit model to data points",
            )
            self.failures.add(model_fit_failure)

    def check_for_model_mse_failure(self, mse_upper_limit: Union[int, float]) -> None:
        """
        If mean squared error of the model is greater than `limit`,
        then flag as a well failure

        Parameters
        -----------
        mse_upper_limit: float or int

        Returns
        --------
        None
            Possibly adds `plaque_assay.failure.WellFailure`
            to `self.failures`.
        """
        if self.mean_squared_error is not None:
            if self.mean_squared_error > mse_upper_limit:
                model_mse_failure = failure.WellFailure(
                    plate="DILUTION SERIES",
                    well=self.sample_name,
                    failure_reason=f"model MSE > {mse_upper_limit} ({self.mean_squared_error:.3f})",
                )
                self.failures.add(model_mse_failure)

    def plot(self) -> List:
        """
        Simple static plot of points and fitted curve.

        Notes
        ------
        This is not used within the automated analysis, but is
        here for convenience when using `plaque_assay` as a separate module.

        Parameters
        -----------

        Returns
        --------
        matplotlib.pyplot.plot
        """
        plt.figure(figsize=[10, 6])
        plt.axhline(y=50, linestyle="--", color="grey")
        plt.scatter(1 / self.data["Dilution"], self.data["Percentage Infected"])
        x = np.logspace(
            np.log10(self.data["Dilution"].min()),
            np.log10(self.data["Dilution"].max()),
            10000,
        )
        if self.model_params is not None:
            curve = stats.dr_4(x, *self.model_params)
            plt.plot(1 / x, curve, linestyle="--", label="4 param dose-response")
            plt.legend(loc="upper left")
            try:
                intersect = stats.find_y_intercept(*self.model_params)
                if intersect:
                    plt.plot(
                        1.0 / intersect, 50, marker="P", color="black", zorder=999,
                    )
            except RuntimeError:
                pass
        plt.xscale("log")
        plt.grid(linestyle=":")
        plt.ylim([0, 110])
        plt.xlim([30, 3000])
        plt.xlabel("Dilution")
        plt.ylabel("Percentage Infected")
        plt.title(self.sample_name)
        return plt
