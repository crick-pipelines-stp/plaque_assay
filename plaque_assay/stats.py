"""
Stats and number crunching functions.
"""

import logging
from typing import NamedTuple, List, Callable, Optional, Union

import numpy as np
import pandas as pd
import scipy.optimize
from numba import jit

from plaque_assay import utils
from plaque_assay import consts


Numeric = Union[int, float]


class Intersect(NamedTuple):
    x: Numeric
    y: Numeric
    error: bool


class ModelParams(NamedTuple):
    top: float
    bottom: float
    ec50: float
    hill_slope: float


class ModelResults(NamedTuple):
    fit_method: str
    result: Numeric
    model_params: Optional[ModelParams]
    mean_squared_error: Optional[float]


def dr_3(x: np.ndarray, top: Numeric, bottom: Numeric, ec50: Numeric) -> np.ndarray:
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
    x: np.ndarray, top: Numeric, bottom: Numeric, ec50: Numeric, hill_slope: Numeric
) -> np.ndarray:
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


def intersect_between_curves(
    x_min: Numeric, x_max: Numeric, curve: np.ndarray, intersect: Numeric = 50
) -> Intersect:
    """Find intersect of two curves.
    Really hacky way of finding intersect of two curves,
    used for QC purposes.
    In this case one of the curves is just a horizontal line where y = 50.
    Parameters
    ----------
    x_min : numeric
    x_max : numeric
    curve : array-like
    intersect : numeric
    Returns
    --------
    Intersect
    """
    x = np.logspace(np.log10(x_min), np.log10(x_max), 10000)
    line = np.full(x.shape, intersect)
    idx_arr = np.argwhere(np.diff(np.sign(line - curve))).flatten()
    error = False
    if len(idx_arr) != 1:
        error = True
        x_intersect = np.nan
        y_intersect = np.nan
    else:
        try:
            idx = int(idx_arr)
            x_intersect = float(x[idx])
            y_intersect = float(curve[idx])
        except (IndexError, ValueError):
            x_intersect = np.nan
            y_intersect = np.nan
            error = True
    result = Intersect(x_intersect, y_intersect, error)
    return result


def find_y_intercept(
    top: float, bottom: float, ec50: float, hillslope: float, y: float = 50.0
) -> float:
    """
    Psuedo-IC50 value.

    Used to find the dilution where percentage infection is 50%.
    NOTE: top and bottom parameters are named incorrectly on purpose
    to mirror the incorrectly named database columns
    """
    return ec50 * (((bottom - top) / (y - top)) - 1.0) ** (1.0 / hillslope)


def non_linear_model(x: Numeric, y: Numeric, func: Callable = dr_4) -> ModelParams:
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
    bounds = ((-0.01, 0, -10, 0), (100, 120, 10, 5))
    popt, *_ = scipy.optimize.curve_fit(
        func, x, y, p0=p0, method="trf", bounds=bounds, maxfev=500
    )
    return ModelParams(*popt)


def model_mse(y_observed: np.ndarray, y_fitted: np.ndarray) -> float:
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
    group: pd.DataFrame, threshold: Numeric, weak_threshold: Numeric
) -> Optional[Numeric]:
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
    # convert dilutions into 40 -> 40_000
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
        if avg[consts.DILUTION_4] <= threshold and avg[consts.DILUTION_3] <= threshold:
            result = "complete inhibition"
    except KeyError:
        # missing this dilution, possibly removed due to high-background
        try:
            # try the next dilution
            if (
                avg[consts.DILUTION_3] <= threshold
                and avg[consts.DILUTION_2] <= threshold
            ):
                result = "complete inhibition"
        except KeyError:
            # if this returns a KeyError aswell we're missing 2 dilutions
            # so we can't get anthing meaningful from the model
            result = "failed to fit model"
    # check for weak inhibition
    try:
        if (
            avg[consts.DILUTION_1] > threshold
            and avg[consts.DILUTION_1] < weak_threshold
        ):
            result = "weak inhibition"
    except KeyError:
        # missing this dilution, possibly removed due to high-background
        try:
            # try the next dilution
            if (
                avg[consts.DILUTION_2] > threshold
                and avg[consts.DILUTION_2] < weak_threshold
            ):
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
    else:
        return None


def calc_heuristics_curve(
    name: str, x: np.ndarray, y: np.ndarray, threshold: Numeric, weak_threshold: Numeric
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
    else:
        return None


def calc_model_results(
    name: str, df: pd.DataFrame, threshold: int = 50, weak_threshold: int = 60
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
    # TODO: fix this monstrosity
    df = df.dropna().sort_values("Dilution")
    x = df["Dilution"].values
    x_min = (1 / consts.DILUTION_4) / 10
    x_max = (1 / consts.DILUTION_1) * 10
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
                intersect = intersect_between_curves(x_min, x_max, y_fitted)
                if intersect.error:
                    logging.error("error caused when finding intersect at y=50")
                    result = utils.result_to_int("failed to fit model")
                    model_params = None
                else:
                    result = 1.0 / intersect.x
                    result = recast_if_out_of_bounds_ic50(result, x, name)
    logging.debug("well %s fitted with method %s", name, fit_method)
    return ModelResults(fit_method, result, model_params, mean_squared_error)


def recast_if_out_of_bounds_ic50(
    result: float, x: np.ndarray, name: str
) -> Union[float, int]:
    """
    if IC50 value falls beyond the tested dilutions then recast as either
    "weak inhibition" or "complete inhibition"
    """
    # IC50 < lowest dilution (1:40)
    if result < 1.0 / x.max():
        logging.info(
            "%s IC50 of %s less than lowest dilution, weak inhibition",
            name,
            result,
        )
        result = utils.result_to_int("weak inhibition")
    # IC50 > highest dilution (1:40_000)
    if result > 1.0 / x.min():
        logging.info(
            "%s IC50 of %s greater than highest dilution, complete inhibition",
            name,
            result,
        )
        result = utils.result_to_int("complete inhibition")
    return result


@jit(nopython=True)
def hampel(x: np.ndarray, k: int, t0: int = 3) -> List:
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
