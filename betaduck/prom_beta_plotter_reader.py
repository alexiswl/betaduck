#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import gzip
from Bio import SeqIO
import concurrent.futures
from datetime import timedelta
import dask.dataframe as dd

import logging

def get_summary_files(summary_dirs):
    summary_files = [os.path.join(summary_dir, summary_file)
                     for summary_dir in summary_dirs
                     for summary_file in os.listdir(summary_dir)
                     if summary_file.endswith(".txt")
                     and "sequencing_summary" in summary_file]

    return summary_files


def get_fastq_files(fastq_dirs):
    fastq_files = [os.path.join(fastq_dir, fastq)
                   for fastq_dir in fastq_dirs
                   for fastq in os.listdir(fastq_dir)
                   if fastq.endswith(".fastq.gz")
                   ]

    return fastq_files


def get_series_from_seq(record):
    """Takes in a seq record and returns dataframe with index"""
    # Series index
    index = ["read_id", "run_id", "sample_id", "read", "channel", "start_time_utc"]
    # Get metadata
    fastq_id = record.id.split()[0].lstrip("@")
    row_as_dict = dict(x.split("=") for x in record.description.split()[1:])

    return pd.Series([fastq_id, row_as_dict['runid'],
                      row_as_dict['sampleid'], row_as_dict['read'],
                      row_as_dict['ch'], row_as_dict['start_time']],
                     index=index)


def get_fastq_dataframe(fastq_file, is_gzipped=True):
    """Use get_series_from_seq in list comprehension to generate dataframe then transpose"""
    # Open fastq file and return list comprehension
    # PD concat merges series as columns, we then transpose.
    try:
        if not is_gzipped:
            with open(fastq_file, "rt") as handle:
                fastq_df = pd.concat([get_series_from_seq(record)
                                      for record in SeqIO.parse(handle, "fastq")],
                                     sort=True,
                                     axis='columns').transpose()
        else:
            with gzip.open(fastq_file, "rt") as handle:
                fastq_df = pd.concat([get_series_from_seq(record)
                                      for record in SeqIO.parse(handle, "fastq")],
                                     sort=True,
                                     axis='columns').transpose()
        # Specify types for each line
        numeric_cols = ["read", "channel"]
        fastq_df[numeric_cols] = fastq_df[numeric_cols].apply(pd.to_numeric, axis='columns')
        # Convert StartTime to date
        fastq_df['start_time_utc'] = pd.to_datetime(fastq_df['start_time_utc'])
        return fastq_df
    except ValueError:
        print("Value error when generating dataframe for %s. Unknown cause of issue." % fastq_file)
        return pd.DataFrame(columns=["read_id", "run_id", "sample_id", "read", "channel", "start_time_utc"])


def read_fastq_datasets(fastq_files, threads):
    # Run in parallel
    datasets = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        iterator = {executor.submit(get_fastq_dataframe, fastq_file, is_gzipped=True):
                    fastq_file for fastq_file in fastq_files}
        for item in concurrent.futures.as_completed(iterator):
            pandas_input = iterator[item]
            datasets.append(item.result())
    # Read in each of the fastq files and retrieve header info
    dataset = pd.concat(datasets, sort=True, ignore_index=True)

    # Return dataset
    return dataset


def read_summary_datasets(sequencing_summary_files, threads):
    # Run in parallel using the dask API
    dataset = dd.read_csv(sequencing_summary_files, sep="\t", header=0)

    # Convert to pandas object
    dataset = dataset.compute()

    # Reset the dtypes for the time columns
    dataset = set_summary_time_dtypes(dataset)

    # Get pass column
    dataset['pass'] = get_pass(dataset)

    # Get qualitative pass
    dataset['qualitative_pass'] = get_qualitative_pass(dataset)

    # Get duration ratio
    dataset['pore_speed'] = get_duration_ratio(dataset)

    # Get events ratio
    dataset['events_ratio'] = get_events_ratio(dataset)

    # Return the object
    return dataset


def set_summary_time_dtypes(dataset):
    """
    :rtype: pd.DataFrame
    """
    # Return the time datasets as appropriate
    dataset.rename(columns={"start_time": "start_time_float",
                            "template_start": "template_time_float"}, inplace=True)

    # Add start_time_float template_time_float
    dataset['start_time_timedelta'] = pd.to_timedelta(dataset['start_time_float'], unit='s')
    dataset['template_start_time_timedelta'] = pd.to_timedelta(dataset['template_time_float'], unit='s')

    return dataset


def sort_summary_dataset(dataset):
    """
    :rtype: pd.DataFrame
    """
    # Sort the dataset by start_time (in seconds)
    return dataset.sort_values(['start_time_timedelta'])


def get_channel_yield(dataset):
    # Get the yield per channel 
    return dataset.groupby(['channel'])['sequence_length_template'].cumsum()


def get_yield(dataset):
    # Get the yield datset
    return dataset['sequence_length_template'].cumsum()


def get_pass(dataset):
    # Determine if sequence passed quality
    return dataset['mean_qscore_template'].apply(lambda x: True if x > 9 else False)


def get_qualitative_pass(dataset):
    # Describe the pass (Passed / Failed)
    return dataset['pass'].apply(lambda x: 'Passed' if x is True else "Failed")


def get_duration_ratio(dataset):
    # Return the length in bases over time in seconds.
    return dataset.apply(
        lambda x: np.nan if x.template_duration == 0 else x.sequence_length_template / x.template_duration,
        axis='columns')


def get_events_ratio(dataset):
    # Return the events-per-base ratio
    return dataset.apply(
        lambda x: np.nan if x.sequence_length_template == 0 else x.num_events / x.sequence_length_template,
        axis='columns')


def trim_dataset(dataset):
    logging.info("Filtering dataset to make plots nice")
    logging.info("Starting with %d reads" % dataset.shape[0])
    # Restrict extremely long reads
    read_length_max_quantile = 0.995
    read_length_query = "sequence_length_template < %d" % dataset['sequence_length_template'].quantile(
        read_length_max_quantile)
    # Events thresold (there's some extreme fail reads up there)
    events_ratio_threshold = 10
    events_ratio_query = "events_ratio < %d" % events_ratio_threshold
    # Some of the times of the fastq are a little whacked.
    time_min_quantile = 0.001
    time_max_quantile = 0.999
    time_min_query = "start_time_float_by_sample > %d" % dataset['start_time_float_by_sample'].quantile(
        time_min_quantile)
    time_max_query = "start_time_float_by_sample < %d" % dataset['start_time_float_by_sample'].quantile(
        time_max_quantile)
    dataset = dataset.query(' & '.join([read_length_query, events_ratio_query, time_min_query, time_max_query]))
    logging.info("Finished filtering with %d reads" % dataset.shape[0])

    return dataset


def convert_sample_time_columns(dataset):
    # Use the utc in the fastq file to work around restarts
    min_start_time = dataset['start_time_utc'].min()
    dataset['start_time_timedelta_by_sample'] = dataset['start_time_utc'].apply(lambda x: x - min_start_time)

    # Convert to float because matplotlib doesn't seem to do timedelta on the x axis well.
    # Need to divide by another timedelta object in order to get float
    dataset['start_time_float_by_sample'] = dataset['start_time_timedelta_by_sample'].apply(lambda x:
                                                                                            x / timedelta(seconds=1))

    # Sort values to start_time_float_by_sample (to assist yield plotting)
    dataset.sort_values(['start_time_float_by_sample'], inplace=True)

    # Reset index values to match
    dataset.reset_index(drop=True, inplace=True)

    # Return
    return dataset


def get_quality_yield(dataset):
    # Get the yield per quality
    return dataset.groupby(['pass'])['sequence_length_template'].cumsum()


def get_quality_count(dataset):
    # Get the yield per quality
    return dataset.groupby(['pass']).cumcount() + 1


def get_read_count(dataset):
    # Get the read count dataset
    return dataset.reset_index()['index'] + 1
