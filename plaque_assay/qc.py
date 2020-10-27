from collections import namedtuple
from string import ascii_uppercase

import stats


def detect_failures(df, threshold=0.7):
    """
    Identify wells where "Cells - Image Region Area [µm²] - Mean per Well"
    is < 70% of the median calculated in "calc_median_all_plates()"
    - return these wells as a list
    - if any are A-H12, fail entire plate
    """
    output = namedtuple("Failures", ["failed_plates", "failed_wells"])
    all_median = stats.calc_median_all_plates(df)
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


def detect_high_background(df):
    """
    Identify wells with high fluorescent background.
    This is detected by wells with over 110% of the experiment median
    in the DAPI channel.
    """
    colname = "Cells - Intensity Image Region DAPI (global) Mean - Mean per Well"
    experiment_median = df[colname].median()
    threshold = experiment_median * 1.1
    failures = []
    columns = ["Well", "PlateNum", colname]
    for idx, well, plate, val in df[columns].itertuples():
        if val > threshold:
            print(f"QC failure (high background): plate:{plate} well:{well}")
            failures.append({"well": well, "plate": plate, "value": val})
    return failures
