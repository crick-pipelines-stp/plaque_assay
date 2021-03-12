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
    """
    A Sample object holds the data for a sample across 4 concentrations
    including replicates.
    """

    def __init__(self, sample_name, data):
        """
        data: pandas.DataFrame
            2 column dataframe: [dilution, value]
        """
        self.sample_name = sample_name
        self.data = data
        self.failures = []
        self.calc_ic50()
        self.is_positive_control = sample_name in POSITIVE_CONTROL_WELLS
        self.check_positive_control()

    def calc_ic50(self):
        """calculate IC50 value"""
        fit_method, ic50, model_params = stats.calc_results_model(
            self.sample_name, self.data
        )
        self.fit_method = fit_method
        self.ic50 = ic50
        self.ic50_pretty = ic50 if ic50 > 0 else utils.int_to_result(ic50)
        self.model_params = model_params

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

    def plot(self):
        """plot dilution with points and curve"""
        plt.figure(figsize=[10, 6])
        plt.axhline(y=50, linestyle="--", color="grey")
        plt.scatter(1 / self.data["Dilution"], self.data["Percentage Infected"])
        x = np.linspace(self.data["Dilution"].min(), self.data["Dilution"].max(), 1000)
        if self.model_params is not None:
            x_min = self.data["Dilution"].min()
            x_max = self.data["Dilution"].max()
            curve = stats.dr_4(x, *self.model_params)
            plt.plot(1 / x, curve, linestyle="--", label="3 param dose-response")
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
