#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import dask.dataframe as dd
import logging

from betaduck.prom_beta_plotter_gen import plot_data, print_stats
from betaduck.prom_beta_plotter_reader import get_summary_files, get_fastq_files
from betaduck.prom_beta_plotter_reader import convert_sample_time_columns, trim_dataset
from betaduck.prom_beta_plotter_reader import read_summary_datasets, read_fastq_datasets
from betaduck.prom_beta_plotter_reader import get_read_count, get_channel_yield
from betaduck.prom_beta_plotter_reader import get_quality_yield, get_yield, get_quality_count

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M')

logger = logging.getLogger()

"""
Sequencing summary file columns
filename    read_id run_id  channel start_time  duration    num_events  template_start  num_events_template
template_duration   sequence_length_template    mean_qscore_template    strand_score_template
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
    parser.add_argument("--summary-dir", type=str, required=True,
                        help="Contains the txt files")
    parser.add_argument("--fastq-dir", type=str, required=True,
                        help="Where are the fastq files")
    parser.add_argument("--plots-dir", type=str, required=True,
                        help="Where do the plots go")
    parser.add_argument("--name", type=str, required=True,
                        help="Titles for plots")
    args = parser.parse_args()
    return args


def main(args):

    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)

    # Check plots_dir exists
    if not os.path.isdir(args.plots_dir):
        os.mkdir(args.plots_dir)

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
    print(summary_datasets.head())

    # Read in fastq_datasets
    logging.info("Reading in fastq datasets")
    fastq_datasets = read_fastq_datasets(fastq_files, args.threads)
    # Rename run_id column
    fastq_datasets.rename(columns={"run_id": "run_id_orig"}, inplace=True)
    print(fastq_datasets.head())

    # Merge summary and fastq datasets
    logging.info("Merging datasets")

    # Merging works faster on indexes
    summary_datasets.set_index(['read_id', 'run_id_orig', 'channel'], inplace=True)
    fastq_datasets.set_index(['read_id', 'run_id_orig', 'channel'], inplace=True)

    # Merge on indexes
    dataset = dd.merge(summary_datasets, fastq_datasets, left_index=True, right_index=True)

    # Now drop the indexes
    dataset.reset_index(['read_id', 'run_id_orig', 'channel'], inplace=True)

    # Drop summary and fastq datasets which will lower the memory requirements of the system.
    del summary_datasets
    del fastq_datasets

    # Add in the start_time_float_by_sample (allows us to later iterate through plots by sample.
    dataset = convert_sample_time_columns(dataset)

    print(dataset)
    print(dataset.columns)

    # Print the non-trimmed stats before proceeding
    print_stats(dataset, args.name+".unfiltered", args.plots_dir)

    # Trim the dataset
    dataset = trim_dataset(dataset, args.plots_dir, args.name)

    # Create a venn diagram here of the reads lost.

    # Reprint the filtered stats
    print_stats(dataset, args.name+".filtered", args.plots_dir)

    # Re-grab the fastq times
    dataset = convert_sample_time_columns(dataset)

    # Get read_count column
    dataset['read_count'] = get_read_count(dataset)

    # Get yield column
    dataset['yield'] = get_yield(dataset)

    # Get the cumulative channel yield
    dataset['channel_yield'] = get_channel_yield(dataset)

    # Get the cumulative quality yield
    dataset['quality_yield'] = get_quality_yield(dataset)

    # Get the cumulative  quality count
    dataset['quality_count'] = get_quality_count(dataset)

    # Plot yields and histograms
    logging.info("Generating plots")
    plot_data(dataset, args.name, args.plots_dir)


if __name__ == "__main__":
    main()
