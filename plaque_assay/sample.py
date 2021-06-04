"""
module docstring
"""
from . import stats
from . import failure
from . import utils
from . import qc_criteria
from .consts import POSITIVE_CONTROL_WELLS

import numpy as np
import matplotlib.pyplot as plt


class Sample:
    """Sample object holds the data for a sample across 4 concentrations
    including replicates.

    Parameters
    -----------
    sample_name : string
        sample name, typically well-label
    data : pandas.DataFrame
        2 column dataframe: [dilution, value]

    Attributes
    ----------
    sample_name : str
        Typically a well label e.g "A01".
    data : pd.DataFrame
        Dataframe containing results filtered to a single sample. Containing
        2 columns: `dilution` and `value`.
    failures : list
        List containing potential `failures.WellFailures`. These include reasons
        such as:

            - invalid IC50 value if the sample is a positive control

            - failing to fit a model

            - discordant replicate values

    is_positive_control : Bool
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

    Methods
    --------
    plot()
    """

    def __init__(self, sample_name, data):
        """
        """
        self.sample_name = sample_name
        self.data = data
        self.failures = []
        self.calc_ic50()
        self.is_positive_control = sample_name in POSITIVE_CONTROL_WELLS
        self.check_positive_control()
        self.check_duplicate_differences()
        self.check_for_model_fit_failure()

    def calc_ic50(self):
        """calculate IC50 value"""
        model_results = stats.calc_results_model(self.sample_name, self.data)
        self.fit_method = model_results.fit_method
        self.ic50 = model_results.result
        self.ic50_pretty = self.ic50 if self.ic50 > 0 else utils.int_to_result(self.ic50)
        self.model_params = model_results.model_params
        self.mean_squared_error = model_results.mean_square_error

    def check_positive_control(self):
        """
        If this sample is a postitive control, then determine
        if the IC50 value is between 500-800. Otherwise it's an
        failure.
        """
        lower_limit = qc_criteria.positive_control_low
        upper_limit = qc_criteria.positive_control_high
        if self.is_positive_control:
            if self.ic50 < lower_limit or self.ic50 > upper_limit:
                positive_control_failure = failure.WellFailure(
                    well=self.sample_name,
                    plate="DILUTION SERIES",
                    reason=f"positive control failure. IC50 = {self.ic50_pretty}",
                )
                self.failures.append(positive_control_failure)

    def check_duplicate_differences(self):
        """
        Want to calculate the differences between duplicates and determine
        if the difference is larger than some threshold. If this is the case
        for more than 2 duplicates in a sample then this should be flagged as
        a QC well failure.
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
                well=self.sample_name,
                plate="DILUTION SERIES",
                reason=f"2 or more duplicates differ by >= {difference_threshold} % infected",
            )
            self.failures.append(duplicate_failure)

    def check_for_model_fit_failure(self):
        """
        If the model has failed to fit this should also be flagged as a QC failure
        """
        if self.ic50_pretty == "failed to fit model" or self.ic50 == -999:
            model_fit_failure = failure.WellFailure(
                well=self.sample_name,
                plate="DILUTION SERIES",
                reason="failed to fit model to data points",
            )
            self.failures.append(model_fit_failure)

    def plot(self):
        """Simple static plot of points and fitted curve

        Returns
        --------
        matplotlib.pyplot.plot
        """
        plt.figure(figsize=[10, 6])
        plt.axhline(y=50, linestyle="--", color="grey")
        plt.scatter(1 / self.data["Dilution"], self.data["Percentage Infected"])
        x = np.linspace(self.data["Dilution"].min(), self.data["Dilution"].max(), 1000)
        if self.model_params is not None:
            x_min = self.data["Dilution"].min()
            x_max = self.data["Dilution"].max()
            curve = stats.dr_4(x, *self.model_params)
            plt.plot(1 / x, curve, linestyle="--", label="4 param dose-response")
            plt.legend(loc="upper left")
            try:
                intersect_x, intersect_y = stats.find_intersect_on_curve(
                    x_min, x_max, curve
                )
                plt.plot(
                    1 / intersect_x, intersect_y, marker="P", color="black", zorder=999
                )
            except RuntimeError:
                # caused by stats.find_intersect_on_curve finding more than
                # 1 intersect point, so a model fit failure anyway, so dont
                # plot the curve
                pass
        plt.xscale("log")
        plt.grid(linestyle=":")
        plt.ylim([0, 110])
        plt.xlim([30, 3000])
        plt.xlabel("Dilution")
        plt.ylabel("Percentage Infected")
        plt.title(self.sample_name)
        return plt
