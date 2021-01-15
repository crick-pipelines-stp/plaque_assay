import os
import textwrap

from . import consts

# Identify wells where "Cells - Image Region Area [µm²] - Mean per Well"
# is < 70% or > 125% of the median calculated in "calc_median_all_plates()"
# - return these wells as a list
# - if any are A-H12, fail entire plate
low_cells_image_region_area_low = 0.7
low_cells_image_region_area_high = 1.25


# Identify wells with high fluorescent background.
#    This is detected by wells with over 110% of the experiment median
#    in the DAPI channel.
high_background_proportion = 1.1

# If this sample is a postitive control, then determine
# if the IC50 value is between 500-800. Otherwise it's an
# failure.
positive_control_low = 40
positive_control_high = 300


# Accetable infection rate, which is calculated as the median
# "Background Subtracted Plaque Area" of the virus-only-wells
# The plate get's flagged is this is out of the expected limits.
infection_rate_low = 0.15
infection_rate_high = 0.8


def print_criteria():
    output = textwrap.dedent(
        f"""
        Quality control criteria used in this analysis
        ===============================================

        Wells identified where "Cells - Image Region Area [µm²] - Mean per Well"
        is < {low_cells_image_region_area_low*100:.1f}% or > {low_cells_image_region_area_high*100:.1f} % of the plate median.
        If any of these wells are in column 12 then flag the entire plate for failure.

        Wells identified with high fluorescent background. These are wells with over
        {high_background_proportion*100:.1f}% of the plate median in feature "Cells - Intensity Image Region DAPI (global) Mean - Mean per Well".

        The IC50 values of the positive control wells have to be between {positive_control_low:.1f} and {positive_control_high:.1f}.

        The infection rate of the virus only wells has to be between {infection_rate_low*100:.1f}% and {infection_rate_high*100:.1f}%.
        This is calculated from "Background Subtracted Plaque Area", and the entire plate is failed if the median
        of the virus-only-wells is outside of these limits.

        virus only wells: {consts.VIRUS_ONLY_WELLS}

        no virus wells: {consts.NO_VIRUS_WELLS}

        positive control wells: {consts.POSITIVE_CONTROL_WELLS}
        """
    )
    return output


def save_qc_criteria(output_dir, experiment):
    text = print_criteria()
    save_path = os.path.join(output_dir, f"qc_criteria_{experiment}.txt")
    with open(save_path, "w") as f:
        f.write(text)
