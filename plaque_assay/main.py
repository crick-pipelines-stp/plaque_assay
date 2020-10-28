import os
import json
from glob import glob
from string import ascii_uppercase
from collections import namedtuple

import pandas as pd
import numpy as np

import data
import qc
import stats


def test_plot(
    df, wells=["E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09", "E10", "E11"]
):
    import matplotlib.pyplot as plt

    df = df[df["Well"].isin(wells)]
    plt.figure(figsize=[10, 6])
    for name, group in df.groupby("Well"):
        group = group.copy()
        group.sort_values("Dilution", inplace=True)
        try:
            model_params = stats.non_linear_model(
                group["Dilution"], group["Percentage Infected"]
            )
        except RuntimeError:
            print(f"Model failure: {name}")
            model_params = None
        plt.axhline(y=50, linestyle="--", color="grey")
        if model_params is not None:
            x_min = group["Dilution"].min()
            x_max = group["Dilution"].max()
            intersect_x, intersect_y = stats.find_intersect_on_curve(
                x_min, x_max, stats.exp_decay, model_params
            )
            x = np.linspace(group["Dilution"].min(), group["Dilution"].max(), 1000)
            plt.plot(
                1 / x, stats.exp_decay(x, *model_params), linestyle="-", label=name
            )
            plt.plot(
                1 / intersect_x, intersect_y, marker="P", color="black", zorder=999
            )
            plt.scatter(1 / group["Dilution"], group["Percentage Infected"])
    plt.xscale("log")
    plt.grid(linestyle=":")
    plt.xlabel("Dilution")
    plt.ylabel("% infected")
    plt.legend(loc="upper left")
    plt.savefig("plot.pdf")


def main():
    df = data.read_data_19()
    failures = qc.detect_failures(df, 0.7)
    df = stats.subtract_background(df)
    df = stats.calc_percentage_infected(df)
    failures_high_background = qc.detect_high_background(df)
    results = stats.calc_results_model(df)
    all_results = {
        "failures": {
            "high_background": failures_high_background,
            "low_cells_image_region_area": {
                "failed_entire_plates": failures.failed_plates,
                "failed_wells": failures.failed_wells,
            },
        },
        "results": results,
    }
    with open("results.json", "w") as f:
        json.dump(all_results, f, indent=4)
    test_plot(df)


if __name__ == "__main__":
    main()
