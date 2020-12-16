import argparse

from plaque_assay.experiment import Experiment
from plaque_assay import data


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="input directory containing the Phenix-created subdirectories",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="output directory where to store the results",
    )
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    dataset = data.read_data_from_directory(args.input)
    experiment = Experiment(dataset)
    experiment.save_results(args.output)
    experiment.save_failures(args.output)
    experiment.save_normalised_data(args.output)


if __name__ == "__main__":
    main()
