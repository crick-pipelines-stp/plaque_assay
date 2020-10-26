import os
import json
from glob import glob
from string import ascii_uppercase
from collections import namedtuple

import pandas as pd
import numpy as np


THIS_DIR = os.path.dirname(os.path.abspath(__file__))

DATA_DIR = os.path.join(THIS_DIR, "data/RuthExpt20201014")

# NOTE: this is temporary until the barcodes are standardised
plate_mapping = {
    17: 1/40,
    19: 1/160,
    21: 1/640,
    23: 1/2560
}

VIRUS_ONLY_WELLS = ("A12", "B12", "C12")
NO_VIRUS_WELLS = ("F12", "G12", "H12")


def row_col_to_well(row, col):
    """convert row-column indices (1-indexed) to well labels"""
    row_str = ascii_uppercase[row-1]
    return f"{row_str}{col:02}"


def get_plate_num(plate_name):
    """get plate number from path"""
    basename = os.path.basename(plate_name)
    return int(basename.split("_")[0])


def read_data(data_dir=DATA_DIR):
    """docstring"""
    # don't yet know how plates will be barcoded to identify the dilution series
    # for now just hard-code this information
    experiments = [os.path.join(DATA_DIR, i) for i in os.listdir(data_dir)]
    plate_name_dict = {os.path.abspath(i): get_plate_num(i) for i in experiments}
    dataframes = []
    for path, plate_num in plate_name_dict.items():
        df = pd.read_csv(
            # NOTE: might not always be Evaluation2
            os.path.join(path, "Evaluation2/PlateResults.txt"),
            skiprows=8,
            sep="\t"
        )
        df["Dilution"] = plate_mapping[plate_num]
        # add well labels to dataframe
        well_labels = []
        for row, col in df[["Row", "Column"]].itertuples(index=False):
            well_labels.append(row_col_to_well(row, col))
        df["Well"] = well_labels
        df["PlateNum"] = plate_num
        dataframes.append(df)
    return pd.concat(dataframes)


def calc_median_all_plates(df):
    """
    calculate the median of "Cells - Image Region Area - Mean per Well"
    for all wells of all plates
    """
    return df["Cells - Image Region Area [µm²] - Mean per Well"].median()


def detect_failures(df, threshold=0.7):
    """
    Identify wells where "Cells - Image Region Area [µm²] - Mean per Well"
    is < 70% of the median calculated in "calc_median_all_plates()"
    - return these wells as a list
    - if any are A-H12, fail entire plate
    """
    output = namedtuple("Failures", ["failed_plates", "failed_wells"])
    all_median = calc_median_all_plates(df)
    col_name = "Cells - Image Region Area [µm²] - Mean per Well"
    thresholded_df = df[(df[col_name] / all_median) < threshold]
    # detect failed wells per plate, store in a dictionary
    failed_well_dict = {}
    for name, group in thresholded_df.groupby("PlateNum"):
        failed_well_dict[name] = group["Well"].tolist()
    # detect failed plates
    failure_wells = [f"{row_letter}12" for row_letter in ascii_uppercase[:8]]
    failed_plates = []
    for plate_num, failed_wells in failed_well_dict.items():
        if any(well in failure_wells for well in failed_wells):
            failed_plates.append(plate_num)
    return output(failed_plates, failed_well_dict)


# calculate median of "Normalised Plaque area" for no virus wells (F12, G12, H12)
# subtract median from the "Normalised Plaque area" for each well and save as "Background subtracted Plaque Area"
def subtract_background(
    df,
    colname="Normalised Plaque area",
    new_colname = "Background subtracted Plaque Area",
    no_virus_wells=NO_VIRUS_WELLS
):
    background = df[df["Well"].isin(NO_VIRUS_WELLS)][colname].median()
    df[new_colname] = df[colname] - background
    return df


# calculate median of "Background subtracted Plaque Area" for virus only wells (A12, B12, C12)
# Divide each well's "Background subtracted Plaque Area" by the median of "Background subtracted ..." to get "Percentage Infected"
def calc_percentage_infected(df):
    colname = "Background subtracted Plaque Area"
    virus_only_median = df[df["Well"].isin(VIRUS_ONLY_WELLS)][colname].median()
    # TODO: if virus_only_median < 0.3, indicate plate fail
    if virus_only_median < 0.3:
        print("plate fail: virus_only_median < 0.3")
    df["Percentage Infected"] = (df[colname] / virus_only_median) * 100
    return df


def simple_threshold(x, y, ec):
    """
    A *really* simple EC50 sort of thing.
    Determine first dilution value at which 'Percentage Infected' is < 50%
    """
    if min(y) > ec:
        # if y never crosses below ec threshold
        return "no inhibition"
    elif y[0] <= ec:
        # if already starts below ec threshold
        return "complete inhibition"
    else:
        # return first dilution at which y crosses below threshold
        return x[np.argmax(y <= ec)]


def exp_decay(x, a, b, c):
    """exponential decay function"""
    return a * np.exp(-b * x) + c


def find_intersect_on_curve(x_min, x_max, func, params, intersect=50):
    """
    Really hacky way of finding dilution where percentage infection is 50%.
    Works by finding the intersection between two curves, in this case
    one of the curves is just a horizontal line where y = 50.
    TODO: solve this mathematically
    """
    x = np.linspace(x_min, x_max, 1000)
    curve = func(x, *params)
    line = np.full(x.shape, intersect)
    idx = np.argwhere(np.diff(np.sign(line - curve))).flatten()
    if len(idx) > 1:
        print(f"Found more than 1 intersect. len = {len(idx)}")
    return x[idx], curve[idx]


def non_linear_model(x, y):
    """
    fit non-linear least squares to the data
    """
    from scipy.optimize import curve_fit
    try:
        popt, pcov = curve_fit(exp_decay, x, y)
        return(popt)
    except RuntimeError as e:
        print(f"Failed to fit model: {e}")
        return None


def calc_results_simple(df, threshold=50):
    """
    simple threshold
    for each well
      iterate through dilutions
      first dilution where "Percentage infected" > 0.5
    """
    output = dict()
    for name, group in df.groupby("Well"):
        group = group.copy()
        group.sort_values("Dilution", inplace=True)
        x = group["Dilution"].values
        y = group["Percentage Infected"].values
        result = simple_threshold(x, y, threshold)
        output[name] = result
    return output


def calc_results_model(df, threshold=50):
    output = dict()
    for name, group in df.groupby("Well"):
        group = group.copy()
        group.sort_values("Dilution", inplace=True)
        x = group["Dilution"].values
        x_min, x_max = x.min(), x.max()
        y = group["Percentage Infected"].values
        if min(y) > threshold:
            # if y never crosses below ec threshold
            result = "no inhibition"
        elif y[0] <= threshold:
            # if already starts below ec threshold
            result = "complete inhibition"
        else:
            # fit non-linear_model
            model_params = non_linear_model(x, y)
            if model_params is not None:
                # model successfully fit
                intersect_x, intersect_y = find_intersect_on_curve(
                    x_min, x_max, exp_decay, model_params
                )
                result = 1 / intersect_x[0]
            else:
                # model failed to fit
                result = "Failed to fit model"
        output[name] = result
    return output


def test_plot(
    df,
    wells=["C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11"]
):
    import matplotlib.pyplot as plt
    df = df[df["Well"].isin(wells)]
    plt.figure(figsize=[10, 6])
    for name, group in df.groupby("Well"):
        group = group.copy()
        group.sort_values("Dilution", inplace=True)
        model_params = non_linear_model(group["Dilution"], group["Percentage Infected"])
        plt.axhline(y=50, linestyle="--", color="grey")
        if model_params is not None:
            x_min = group["Dilution"].min()
            x_max = group["Dilution"].max()
            intersect_x, intersect_y = find_intersect_on_curve(x_min, x_max, exp_decay, model_params)
            x = np.linspace(group["Dilution"].min(), group["Dilution"].max(), 1000)
            plt.plot(1/x, exp_decay(x, *model_params), linestyle="-", label=name)
            plt.plot(1/intersect_x, intersect_y, marker="P", color="black", zorder=999)
            plt.scatter(1/group["Dilution"], group["Percentage Infected"])
    plt.xscale("log")
    plt.grid(linestyle=":")
    plt.xlabel("Dilution")
    plt.ylabel("% infected")
    plt.legend(loc="upper left")
    plt.savefig("plot.pdf")


def main():
    df = read_data()
    failures = detect_failures(df, 0.7)
    df = subtract_background(df)
    df = calc_percentage_infected(df)
    results = calc_results_model(df)
    all_results = {
        "failures": {

            "plates": failures.failed_plates,
            "wells": failures.failed_wells
        },
        "results": results,
    }
    with open("results.json", "w") as f:
        json.dump(all_results, f, indent=4)
    test_plot(df)



if __name__ == "__main__":
    main()
