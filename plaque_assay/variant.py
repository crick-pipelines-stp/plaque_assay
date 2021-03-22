import math
import os
import string


class MultipleVariantError(Exception):
    """class docstring"""
    def __init__(self, message="Multiple variants detected."):
        self.message = message
        super().__init__(self.message)


class Variant:
    """class docstring"""

    def __init__(self):
        self.mapper = self.create_mapper()

    def create_mapper(self):
        """docstring"""
        variant_dict = dict()
        for i in range(1, 27):
            letter_int = math.ceil(i / 2) - 1
            variant_dict[i] = string.ascii_lowercase[letter_int]
        return variant_dict

    def get_variant_from_barcode(self, barcode):
        """docstring"""
        variant_int = int(barcode[1:3])
        return self.mapper[variant_int]

    def get_variant_from_barcodes(self, barcodes):
        """docstring"""
        variant_set = set()
        for barcode in barcodes:
            variant = self.get_variant_from_barcode(barcode)
            variant_set.add(variant)
        if len(variant_set) > 1:
            raise MultipleVariantError
        return next(iter(variant_set))

    def get_variant_from_experiment_df(self, experiment_df):
        """docstring"""
        barcodes = experiment_df["Plate_barcode"].values
        variant = self.get_variant_from_barcodes(barcodes)
        return variant

    def get_variant_from_plate_path_list(self, plate_list):
        """docstring"""
        just_plates = [os.path.basename(i) for i in plate_list]
        barcodes = [i.split("__")[0] for i in just_plates]
        return self.get_variant_from_barcodes(barcodes)
