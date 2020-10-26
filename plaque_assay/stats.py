import numpy as np
import scipy.optimize


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


def non_linear_model(x, y):
    """
    fit non-linear least squares to the data
    """
    try:
        popt, pcov = scipy.optimize.curve_fit(exp_decay, x, y)
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


