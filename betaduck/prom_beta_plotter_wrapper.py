#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import logging

from betaduck.prom_beta_plotter_gen import plot_data
from betaduck.prom_beta_plotter_reader import get_summary_files
from betaduck.prom_beta_plotter_reader import get_fastq_files
from betaduck.prom_beta_plotter_reader import read_summary_datasets, read_fastq_datasets

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

"""
Sequencing summary file columns
filename	read_id	run_id	channel	start_time	duration	num_events	template_start	num_events_template
template_duration	sequence_length_template	mean_qscore_template	strand_score_template
"""

"""
fastq file columns
"fastq_id", "sample_id", "read", "channel", "start_time_utc"
"""


def get_args():
    """Get arguments from commandline"""
    """
    Two simple arguments.
    1. Path to Sequencing Summary Directory
    """
    parser = argparse.ArgumentParser(description="Plot a run as it is going")
    parser.add_argument("--summary_dir", type=str, required=True,
                        help="Contains the txt files")
    parser.add_argument("--fastq_dir", type=str, required=True,
                        help="Where are the fastq files")
    parser.add_argument("--plots_dir", type=str, required=True,
                        help="Where do the plots go")
    parser.add_argument("--name", type=str, required=True,
                        help="Titles for plots")
    args = parser.parse_args()
    return args


def main(args):

    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)

    # Get summary files
    summary_files = get_summary_files([summary_dir
                                       for summary_dir in args.summary_dir.split(",")])

    # Get fastq files
    fastq_files = get_fastq_files([fastq_dir
                                   for fastq_dir in args.fastq_dir.split(",")])

    # Reduce thread count unless already 1.
    threads = 1 if args.threads == 1 else args.threads - 1
    logging.info("Given we need to take of the parent script, running %d jobs in parallel" % threads)

    # Read in summary datasets
    logging.info("Reading in summary datasets")
    summary_datasets = read_summary_datasets(summary_files, args.threads)

    # Read in fastq_datasets
    logging.info("Reading in fastq datasets")
    fastq_datasets = read_fastq_datasets(fastq_files, args.threads)

    # Merge summary and fastq datasets
    logging.info("Merging datasets")
    dataset = pd.merge(summary_datasets, fastq_datasets, on=['read_id', 'run_id', 'channel'])

    # Check plots_dir exists
    if not os.path.isdir(args.plots_dir):
        os.mkdir(args.plots_dir)

    # Plot yields and histograms
    logging.info("Generating plots")
    plot_data(dataset, args.name, args.plots_dir)


if __name__ == "__main__":
    main()
