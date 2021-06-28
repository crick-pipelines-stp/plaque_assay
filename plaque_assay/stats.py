"""
Stats and number crunching functions.
"""

import logging
from numbers import Number
from typing import NamedTuple, List, Callable, Optional

import numpy as np
import pandas as pd
import scipy.optimize

from . import utils


class Intersect(NamedTuple):
    x: Number
    y: Number


class ModelParams(NamedTuple):
    top: float
    bottom: float
    ec50: float
    hill_slope: float


class ModelResults(NamedTuple):
    fit_method: str
    result: Number
    model_params: Optional[ModelParams]
    mean_squared_error: float


def dr_3(x: np.array, top: Number, bottom: Number, ec50: Number) -> np.array:
    """3 parameter dose response curve

    Parameters
    -----------
    x : array-like
    top : numeric
    bottom : numeric
    ec50 : numeric

    Returns
    --------
    array-like
        same shape as input `x`
    """
    return bottom + x * (top - bottom) / (ec50 + x)


def dr_4(
    x: np.array, top: Number, bottom: Number, ec50: Number, hill_slope: Number
) -> np.array:
    """4 parameter dose response curve

    Parameters
    -----------
    x : array-like
    top : numeric
    bottom : numeric
    ec50 : numeric
    hill_slope : numeric

    Returns
    --------
    array-like
        same shape as input `x`
    """
    # during optimisation numpy doesn't like raising negative numbers
    # to a fractional power, this stops it spamming up the logs
    with np.errstate(invalid="ignore"):
        return (bottom - top) / (1 + (x / ec50) ** hill_slope)


def find_intersect_on_curve(
    x_min: Number, x_max: Number, curve: np.array, intersect: Number = 50
) -> Optional[Intersect]:
    """Find intersect of two curves.

    Really hacky way of finding intersect of two curves,
    used for finding dilution where percentage infection is 50%.
    In this case one of the curves is just a horizontal line where y = 50.
    TODO: solve this mathematically

    Parameters
    ----------
    x_min : numeric
    x_max : numeric
    curve : array-like
    intersect : numeric

    Returns
    --------
    Intersect or None
    """
    x = np.logspace(np.log10(x_min), np.log10(x_max), 10000)
    line = np.full(x.shape, intersect)
    idx = np.argwhere(np.diff(np.sign(line - curve))).flatten()
    if len(idx) > 1:
        logging.error(f"Found more than 1 intersect. len = {len(idx)}")
        return None
    # we have a numpy array of length 1, so just get the value
    idx = idx[0]
    return Intersect(x[idx], curve[idx])


def non_linear_model(x: Number, y: Number, func: Callable = dr_4) -> ModelParams:
    """
    fit non-linear least squares to the data

    Parameters
    ----------
    x : numeric
    y : numeric
    func : Callable

    Returns
    --------
    `plaque_assay.stats.ModelParams`
    """
    # initial guess at sensible parameters
    p0 = [0, 100, 0.015, 1]
    bounds = ((-0.01, 0, -10, -10), (100, 120, 10, 10))
    popt, *_ = scipy.optimize.curve_fit(
        func, x, y, p0=p0, method="trf", bounds=bounds, maxfev=500
    )
    return ModelParams(*popt)


def model_mse(y_observed: np.array, y_fitted: np.array) -> float:
    """
    Mean squared error between the observed
    and the estimated values.

    Parameters
    -----------
    y_observed : numeric
        observed y-values
    y_estimated : numeric
        estimated y-values from model fit

    Returns
    -------
    float
    """
    assert y_observed.shape == y_fitted.shape
    return np.nanmean((y_observed - y_fitted) ** 2)


def calc_heuristics_dilutions(
    group: pd.DataFrame, threshold: Number, weak_threshold: Number
) -> Number:
    """Simple heuristics based on the values without model fitting

    Parameters
    ----------
    group
    threshold
    weak_threshold

    Returns
    --------
    numeric
        IC50 value, or negative integer indicating an error code
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
        # if 2 most dilute values are below threshold, then
        # label it as complete inhibition
        if avg[2560] <= threshold and avg[640] <= threshold:
            result = "complete inhibition"
    except KeyError:
        # missing this dilution, possibly removed due to high-background
        try:
            # try the next dilution
            if avg[640] <= threshold and avg[160] <= threshold:
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


def calc_heuristics_curve(
    name: str, x: np.array, y: np.array, threshold: Number, weak_threshold: Number
) -> Optional[int]:
    """
    heuristics based on the model fit where we cannot calculate the intercept
    at the threshold

    Parameters
    -----------
    name : str
        name of sample or well
    x : array-like
    y : array-like
    threshold : numeric
    weak_threshold : numeric

    Returns
    -------
    int or None
        If `None`, then an heuristic was not found.
        Otherwise returns a negative integer code.
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


def calc_model_results(
    name: str, df: pd.DataFrame, threshold: Number = 50, weak_threshold: Number = 60
) -> ModelResults:
    """
    Try simple heuristics first without model fitting.
    Once a model is fitted try heuristics based on the fitted curve.
    Then calculate the value based on the intercept where the curve = threshold.

    Parameters
    -----------
    name : str
    df : pandas.DataFrame
    threshold : numeric
    weak_threshold : numeric

    Returns
    --------
    `plaque_assay.stats.ModelResults`
    """
    df = df.dropna().sort_values("Dilution")
    x = df["Dilution"].values
    x_min = 0.0000390625
    x_max = 0.25
    x_interpolated = np.logspace(np.log10(x_min), np.log10(x_max), 10000)
    y = df["Percentage Infected"].values
    model_params = None
    mean_squared_error = None
    heuristic = calc_heuristics_dilutions(df, threshold, weak_threshold)
    if heuristic is not None:
        result = heuristic
        fit_method = "heuristic"
    else:
        # fit non-linear_model
        try:
            model_params = non_linear_model(x, y)
        except RuntimeError:
            model_params = None
            result = utils.result_to_int("failed to fit model")
        fit_method = "model fit"
        if model_params is not None:
            # predicted y-values for interpolated x-values, useful to generate curve
            y_fitted = dr_4(x_interpolated, *model_params)
            # predicted y-values only for dilution x-values, useful for MSE
            # calculation
            y_hat = dr_4(x, *model_params)
            mean_squared_error = model_mse(y_hat, y)
            if mean_squared_error > 99999:
                logging.warning("MSE > 99999, clipped to 99999 to fit in database")
                mean_squared_error = 99999
            curve_heuristics = calc_heuristics_curve(
                name, x_interpolated, y_fitted, threshold, weak_threshold
            )
            if curve_heuristics is not None:
                result = curve_heuristics
            else:
                try:
                    intersect = find_intersect_on_curve(x_min, x_max, y_fitted)
                    result = 1 / intersect.x
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
    return ModelResults(fit_method, result, model_params, mean_squared_error)


def hampel(x: np.array, k: int, t0: Number = 3) -> List:
    """Hampel's outlier test

    Adapted from hampel function in R package pracma.

    Parameters
    -----------
    x : 1-d array
        numbers to be filtered
    k : int
        number of items in window/2 (# forward and backward wanted to capture in median filter)
    t0 : numeric (float or int)
        number of standard deviations to use; 3 is default

    Returns
    -------
    array-like
        indices in x of outliers
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
