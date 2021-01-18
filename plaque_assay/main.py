import argparse
import logging
import os
import sqlite3

from plaque_assay.experiment import Experiment
from plaque_assay.consts import UNWANTED_METADATA
from plaque_assay import data
from plaque_assay import qc_criteria


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
    logging.basicConfig(
        filename=f"{os.path.expanduser('~')}/plaque_analysis.log",
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
    )
    logging.info("Started analysis with args %s", args)
    # create database for dashboard if it doesn't already exist
    db_path = "./plaque_assay_results.sqlite"
    conn = sqlite3.connect(db_path)
    dataset = data.read_data_from_directory(args.input)
    experiment = Experiment(dataset)
    logging.info("Experiment named as: %s", experiment.experiment_name)
    # save concatenated "raw" data
    dataset.drop(UNWANTED_METADATA, axis=1).to_csv(
        os.path.join(args.output, f"plateResults_{experiment.experiment_name}.csv"),
        index=False,
    )
    indexfiles = data.read_indexfiles_from_directory(args.input)
    indexfiles.to_csv(
        os.path.join(args.output, f"indexfiles_{experiment.experiment_name}.csv"),
        index=False,
    )
    results = experiment.get_results_as_dataframe()
    results["experiment"] = experiment.experiment_name
    results.to_sql("results", con=conn, if_exists="append", index=False)
    experiment.save_results_as_dataframe(args.output)
    experiment.save_failures_as_dataframe(args.output)
    percentage_infected_df = experiment.get_percentage_infected_dataframe()
    percentage_infected_df.to_sql("percentage_infected", con=conn, if_exists="append")
    failures = experiment.get_failures_as_dataframe()
    failures.to_sql("failures", con=conn, if_exists="append", index=False)
    experiment.save_normalised_data(args.output)
    model_parameters = experiment.get_model_parameters_as_dataframe()
    model_parameters.to_sql("parameters", con=conn, if_exists="append")


def run(plate_list, output_dir, plot=True):
    db_path = "/usr/src/app/test_output/plaque_assay_results.sqlite"
    conn = sqlite3.connect(db_path)
    dataset = data.read_data_from_list(plate_list)
    experiment = Experiment(dataset)
    # save concatenated "raw" data
    # save concatenated "raw" data
    dataset.drop(UNWANTED_METADATA, axis=1).to_csv(
        os.path.join(output_dir, f"plateResults_{experiment.experiment_name}.csv"),
        index=False,
    )
    data.read_indexfiles_from_list(plate_list).to_csv(
        os.path.join(output_dir, f"indexfile_{experiment.experiment_name}.csv"),
        index=False,
    )
    results = experiment.get_results_as_dataframe()
    results["experiment"] = experiment.experiment_name
    results.to_sql("results", con=conn, if_exists="append", index=False)
    experiment.save_results_as_dataframe(output_dir)
    experiment.save_failures_as_dataframe(output_dir)
    percentage_infected_df = experiment.get_percentage_infected_dataframe()
    percentage_infected_df.to_sql("percentage_infected", con=conn, if_exists="append")
    failures = experiment.get_failures_as_dataframe()
    failures.to_sql("failures", con=conn, if_exists="append", index=False)
    experiment.save_normalised_data(output_dir)
    model_parameters = experiment.get_model_parameters_as_dataframe()
    model_parameters.to_sql("parameters", con=conn, if_exists="append")
    qc_criteria.save_qc_criteria(output_dir, experiment.experiment_name)


if __name__ == "__main__":
    main()
