import logging

import numpy as np
import scipy.optimize

from . import utils


def dr_3(x, top, bottom, ec50):
    """3 parameter dose response curve"""
    return bottom + x * (top - bottom) / (ec50 + x)


def find_intersect_on_curve(x_min, x_max, curve, intersect=50):
    """
    Really hacky way of finding intersect of two curves,
    used for finding dilution where percentage infection is 50%.
    In this case one of the curves is just a horizontal line where y = 50.
    TODO: solve this mathematically
    """
    x = np.linspace(x_min, x_max, 1000)
    line = np.full(x.shape, intersect)
    idx = np.argwhere(np.diff(np.sign(line - curve))).flatten()
    if len(idx) > 1:
        logging.error(f"Found more than 1 intersect. len = {len(idx)}")
        return None
    return x[idx], curve[idx]


def non_linear_model(x, y, func=dr_3):
    """
    fit non-linear least squares to the data
    """
    # initial guess at sensible parameters
    p0 = [0, 100, 0.015]
    bounds = ((-0.01, 0, -10), (100, 120, 10))
    popt, pcov = scipy.optimize.curve_fit(
        func, x, y, p0=p0, method="trf", bounds=bounds, maxfev=500
    )
    return popt, pcov


def calc_heuristics_dilutions(group, threshold, weak_threshold):
    """
    simple heuristics based on the values without model fitting
    """
    result = None
    avg = group.groupby("Dilution")["Percentage Infected"].mean()
    # convert dilutions into 40 -> 2560
    avg.index = 1 / avg.index
    avg.index = avg.index.astype(int)
    # round to nearest 10
    avg.index = [round(i, -1) for i in avg.index]
    # for complete inhibition
    if all(avg.values <= threshold):
        result = "complete inhibition"
    try:
        if avg[2560] <= threshold and avg[2560] <= threshold:
            result = "complete inhibition"
    except KeyError:
        # missing this dilution, possibly removed due to high-background
        try:
            # try the next dilution
            if avg[640] <= threshold and avg[640] <= threshold:
                result = "complete inhibition"
        except KeyError:
            # if this returns a KeyError aswell we're missing 2 dilutions
            # so we can't get anthing meaningful from the model
            result = "failed to fit model"
    # check for weak inhibition
    try:
        if avg[40] > threshold and avg[40] < weak_threshold:
            result = "weak inhibition"
    except KeyError:
        # missing this dilution, possibly removed due to high-background
        try:
            # try the next dilution
            if avg[160] > threshold and avg[160] < weak_threshold:
                result = "weak inhibition"
        except KeyError:
            # if this returns a KeyError aswell we're missing 2 dilutions
            # so we can't get anthing meaningful from the model
            result = "failed to fit model"
    # check for no inhibition
    if all(avg.values > weak_threshold):
        result = "no inhibition"
    if result:
        return utils.result_to_int(result)


def calc_heuristics_curve(name, x, y, threshold, weak_threshold):
    """
    heuristics based on the model fit where we cannot calculate the intercept
    at the threshold
    """
    result = None
    # look for sharp changes in the curve shape indicating a bad fit
    outliers = hampel(y, 5)
    if outliers:
        result = "failed to fit model"
        logging.warning("well %s model failed due to hampel outliers on curve", name)
    # look for times when the curve doesn't reach below threshold but
    # drops below weak_threshold indicated "weak inhibition"
    if min(y) > threshold and min(y) < weak_threshold:
        # determine minimum is on the side we would expect (1:40)
        idx_min = np.argmin(y)
        # checking for greater than as the actual values are inverted because
        # we're using 1/dilution
        if idx_min > len(x) / 4:
            result = "weak inhibition"
        else:
            result = "failed to fit model"
    if result:
        return utils.result_to_int(result)


def calc_results_model(name, df, threshold=50, weak_threshold=60):
    """
    Try simple heuristics first without model fitting.
    Once a model is fitted try heuristics based on the fitted curve.
    Then calculate the value based on the intercept where the curve = threshold.
    """
    df = df.sort_values("Dilution")
    x = df["Dilution"].values
    x_min = 0.0000390625
    x_max = 0.25
    x_interpolated = np.linspace(x_min, x_max, 1000)
    y = df["Percentage Infected"].values
    model_params = None
    heuristic = calc_heuristics_dilutions(df, threshold, weak_threshold)
    if heuristic is not None:
        result = heuristic
        fit_method = "heuristic"
    else:
        # fit non-linear_model
        try:
            model_params, _ = non_linear_model(x, y)
        except RuntimeError:
            model_params = None
            result = utils.result_to_int("failed to fit model")
        fit_method = "model fit"
        if model_params is not None:
            # model successfully fit
            dr_curve = dr_3(x_interpolated, *model_params)
            curve_heuristics = calc_heuristics_curve(
                name, x_interpolated, dr_curve, threshold, weak_threshold
            )
            if curve_heuristics is not None:
                result = curve_heuristics
            else:
                try:
                    intersect_x, intersect_y = find_intersect_on_curve(
                        x_min, x_max, dr_curve
                    )
                    result = 1 / intersect_x[0]
                    if result < 1 / x.max():
                        logging.info(
                            "%s IC50 of %s less than lowest dilution, weak inhibition",
                            name,
                            result,
                        )
                        result = utils.result_to_int("weak inhibition")
                except (IndexError, RuntimeError) as e:
                    logging.error("during model fitting: %s", e)
                    result = utils.result_to_int("failed to fit model")
    logging.debug("well %s fitted with method %s", name, fit_method)
    return fit_method, result, model_params


def hampel(x, k, t0=3):
    """
    adapted from hampel function in R package pracma
        x: 1-d numpy array of numbers to be filtered
        k: number of items in window/2 (# forward and backward wanted to capture in median filter)
        t0: number of standard deviations to use; 3 is default
    """
    n = len(x)
    L = 1.4826
    indices = []
    for i in range((k + 1), (n - k)):
        if np.isnan(x[(i - k) : (i + k + 1)]).all():
            continue
        x0 = np.nanmedian(x[(i - k) : (i + k + 1)])
        S0 = L * np.nanmedian(np.abs(x[(i - k) : (i + k + 1)] - x0))
        if np.abs(x[i] - x0) > t0 * S0:
            indices.append(i)
    return indices


def calc_median_all_plates(df):
    """
    calculate the median of "Cells - Image Region Area - Mean per Well"
    for all wells of all plates
    """
    median = df["Cells - Image Region Area [µm²] - Mean per Well"].median()
    logging.info("median 'Cells = Image Region Area' for all plates %s", median)
    return median
