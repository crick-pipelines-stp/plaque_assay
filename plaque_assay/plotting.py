import os

import numpy as np
import matplotlib.pyplot as plt

from . import stats


def test_plot(df):
    if not os.path.isdir("../output"):
        os.makedirs("../output")
    for name, group in df.groupby("Well"):
        plt.figure(figsize=[10, 6])
        group = group.copy()
        group.sort_values("Dilution", inplace=True)
        try:
            model_params, model_cov = stats.non_linear_model(
                group["Dilution"], group["Percentage Infected"]
            )
        except RuntimeError:
            model_params = None
        plt.axhline(y=50, linestyle="--", color="grey")
        plt.scatter(1 / group["Dilution"], group["Percentage Infected"])
        x = np.linspace(group["Dilution"].min(), group["Dilution"].max(), 1000)
        if model_params is not None:
            x_min = group["Dilution"].min()
            x_max = group["Dilution"].max()
            curve =stats.exp_decay(x, *model_params)
            intersect_x, intersect_y = stats.find_intersect_on_curve(
                x_min, x_max, curve
            )
            plt.plot(
                1 / x,
                stats.exp_decay(x, *model_params),
                linestyle="--",
                label="exponential function (1st choice)",
            )
            plt.plot(
                1 / intersect_x, intersect_y, marker="P", color="black", zorder=999
            )
        poly_fit = stats.polynomial_model(
            group["Dilution"], group["Percentage Infected"]
        )
        plt.plot(1 / x, poly_fit(x), label="quadratic function (2nd choice)")
        plt.xscale("log")
        plt.grid(linestyle=":")
        plt.ylim([0, 110])
        plt.xlim([30, 3000])
        plt.xlabel("Dilution")
        plt.ylabel("% infected")
        plt.legend(loc="upper left")
        plt.title(name)
        plt.savefig(f"../output/plot_{name}.pdf")
        plt.close()
