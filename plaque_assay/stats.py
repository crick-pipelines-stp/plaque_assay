import numpy as np
import numpy.polynomial.polynomial as poly
import scipy.optimize

from . import data


def exp_decay(x, a, b, c):
    """exponential decay function"""
    return a * np.exp(-b * x) + c


def find_intersect_on_curve(x_min, x_max, func, params, intersect=50):
    """
    Really hacky way of finding intersect of two curves,
    used for finding dilution where percentage infection is 50%.
    In this case one of the curves is just a horizontal line where y = 50.
    TODO: solve this mathematically
    """
    x = np.linspace(x_min, x_max, 1000)
    curve = func(x, *params)
    line = np.full(x.shape, intersect)
    idx = np.argwhere(np.diff(np.sign(line - curve))).flatten()
    if len(idx) > 1:
        print(f"Found more than 1 intersect. len = {len(idx)}")
    return x[idx], curve[idx]


def find_intersect_on_curve_polynomial(x_min, x_max, poly_fit, intersect=50):
    x = np.linspace(x_min, x_max, 1000)
    line = np.full(x.shape, intersect)
    curve = poly_fit(x)
    idx = np.argwhere(np.diff(np.sign(line - curve))).flatten()
    if len(idx) > 1:
        print(f"Found more than 1 intersect. len = {len(idx)}")
    return x[idx], curve[idx]


def non_linear_model(x, y):
    """
    fit non-linear least squares to the data
    """
    # initial guess at sensible parameters
    p0 = [100, 100, 0]
    popt, pcov = scipy.optimize.curve_fit(
        exp_decay, x, y, p0=p0, method="lm", maxfev=1000,
    )
    return popt, pcov


def polynomial_model(x, y, d=2):
    """
    fit polynomial model, return function which can be used on new data.
    """
    popt = poly.polyfit(x, y, d)
    return poly.Polynomial(popt)


def confidence_bands(func, model_params, model_cov, x=None):
    """
    returns upper and lower confidence bands
    """
    sigma = np.sqrt(np.diagonal(model_cov))
    if x is None:
        x = np.linspace(40, 2560, 1000)
    bound_upper = func(x, *(model_params + sigma))
    bound_lower = func(x, *(model_params - sigma))
    return bound_lower, bound_upper


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


def calc_results_model(df, threshold=50, weak_threshold=60):
    output = dict()
    for name, group in df.groupby("Well"):
        group = group.copy()
        group.sort_values("Dilution", inplace=True)
        x = group["Dilution"].values
        x_min, x_max = x.min(), x.max()
        y = group["Percentage Infected"].values
        if min(y) > weak_threshold:
            # if y never crosses below ec threshold
            result = "no inhibition"
        elif min(y) < weak_threshold and min(y) > threshold:
            # if 1/40 dilution reaches 60% percent infected
            result = "weak inhibition"
        elif y[0] <= threshold:
            # if already starts below ec threshold
            result = "complete inhibition"
        else:
            # fit non-linear_model
            try:
                model_params, _ = non_linear_model(x, y)
            except RuntimeError:
                print(f"Model fit failure: {name}")
                model_params = None
            if model_params is not None:
                # model successfully fit
                intersect_x, intersect_y = find_intersect_on_curve(
                    x_min, x_max, exp_decay, model_params
                )
                try:
                    result = 1 / intersect_x[0]
                except IndexError:
                    result = "Failed to fit model"
            else:
                # model failed to fit
                poly_fit = polynomial_model(x, y)
                result = "Failed to fit model"
                intersect_x, intersect_y = find_intersect_on_curve_polynomial(
                    x_min, x_max, poly_fit
                )
                try:
                    result = 1 / intersect_x[0]
                except IndexError:
                    result = "Failed to fit model"
        output[name] = result
    return output


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


def calc_percentage_infected(df):
    colname = "Background subtracted Plaque Area"
    virus_only_median = df[df["Well"].isin(data.VIRUS_ONLY_WELLS)][colname].median()
    # TODO: if virus_only_median < 0.3, indicate plate fail
    if virus_only_median < 0.3:
        print("plate fail: virus_only_median < 0.3")
    df["Percentage Infected"] = (df[colname] / virus_only_median) * 100
    return df


def subtract_background(
    df,
    colname="Normalised Plaque area",
    new_colname="Background subtracted Plaque Area",
    no_virus_wells=data.NO_VIRUS_WELLS,
):
    background = df[df["Well"].isin(no_virus_wells)][colname].median()
    df[new_colname] = df[colname] - background
    return df


def calc_median_all_plates(df):
    """
    calculate the median of "Cells - Image Region Area - Mean per Well"
    for all wells of all plates
    """
    return df["Cells - Image Region Area [µm²] - Mean per Well"].median()
