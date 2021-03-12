import os
import textwrap

from . import consts

"""
1. Calculate median of “Cells - Image Region Area [¬µm¬≤] - Mean per Well”
   for all wells of all plates

2. Identify wells where “Cells - Image Region Area [¬µm¬≤] - Mean per Well”
   is < 65% or > 150% of the median calculated in step 1.
       a. Flag these wells as a list for ‘well fail’
       b. If any of these wells are A6:H6 and A12:H12, flag ‘control well fail’
       c. If >8 wells fail, flag ‘possible plate fail – check DAPI plate image’
"""
# Identify wells where "Cells - Image Region Area [µm²] - Mean per Well"
# is < 65% or > 150% of the median calculated in "calc_median_all_plates()"
# - return these wells as a list
# - if any are A-H12, fail entire plate
low_cells_image_region_area_low = 0.6
low_cells_image_region_area_high = 1.5

"""
3. Calculate median of “Normalised Plaque Area” for no virus wells which are
   F12, G12, H12.

4. Subtract median calculated in step 3 from the “Normalised Plaque Area” for
   each well and save as “Background Subtracted Plaque Area”

5. Calculate median of “Background Subtracted Plaque Area” for virus only wells
   which are B6:G6 and A12:C12.
       a. If median calculated is < 0.4 or >0.8 then display this value and
          flag ‘plate fail due to infection outside optimal range’.
"""
# Accetable infection rate, which is calculated as the median
# "Background Subtracted Plaque Area" of the virus-only-wells
# The plate get's flagged is this is out of the expected limits.
infection_rate_low = 0.4
infection_rate_high = 0.8

"""
6. Divide each well’s “Background Subtracted Plaque Area” by the median calculated
   in step 5 to get “Percentage Infected”.

7. Calculated dilution where “Percentage Infected” for a particular well would be 0.5 or 50% (IC50).
    a. If no inhibition at any dilution, then return ‘no inhibition’
    b. If > 50% inhibition at 1:2560 dilution, return ‘complete inhibition’
    c. If between 40-50% inhibition at 1:40, return ‘weak inhibition’

8. Check IC50 for control wells A6, H6, D12, and E12.
    a. If IC50 is not between 300-800, then flag “positive control IC50 outside expected range”.
"""
# If this sample is a postitive control, then determine
# if the IC50 value is between specified values. Otherwise it's an
# failure.
positive_control_low = 300
positive_control_high = 800


def print_criteria():
    output = textwrap.dedent(
        f"""
        Quality control criteria used in this analysis
        ===============================================

        Wells identified where "Cells - Image Region Area [µm²] - Mean per Well"
        is < {low_cells_image_region_area_low*100:.1f}% or > {low_cells_image_region_area_high*100:.1f} % of the plate median.
        If any of these wells are in column 12 then flag the entire plate for failure.

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
