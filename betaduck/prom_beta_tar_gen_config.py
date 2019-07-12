#!/usr/bin/env python3

import argparse
import os
import yaml
import re
import pandas as pd
import json
import logging
import subprocess
import h5py
import sys
import shutil

"""
Generate a yaml file used to run each of the alpha_light python commands
Re-distribute the sequencing summary file as a smaller subset.
"""

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

zfill = 5

config_columns = ['zfill_num', 'rand_id',
                  'fast5_pass_file', 'fast5_pass_file_renamed',
                  'fast5_fail_file', 'fast5_fail_file_renamed',
                  'fastq_pass_file', 'fastq_pass_file_renamed',
                  'fastq_fail_file', 'fastq_fail_file_renamed',
                  'md5_fast5_pass', 'md5_fastq_pass',
                  'md5_fast5_fail', 'md5_fastq_fail']


def get_all_files(rand_id, summary_df):
    """
    :rtype: pd.DataFrame
    :param rand_id: string
    :param summary_df: pd.DataFrame

    Run directory looks like this
    ./fast5_pass
        PAD29258_b799457dcb4fece6795e4bca6ae64348f2687c97_0.fast5
        PAD29258_b799457dcb4fece6795e4bca6ae64348f2687c97_1000.fast5
    ./fast5_fail
        PAD29258_b799457dcb4fece6795e4bca6ae64348f2687c97_0.fast5
        PAD29258_b799457dcb4fece6795e4bca6ae64348f2687c97_1000.fast5
    ./fastq_fail
        PAD29258_b799457dcb4fece6795e4bca6ae64348f2687c97_0.fastq
        PAD29258_b799457dcb4fece6795e4bca6ae64348f2687c97_1000.fastq
    ./fast5_skip
    ./fastq_pass
        PAD29258_b799457dcb4fece6795e4bca6ae64348f2687c97_0.fastq
        PAD29258_b799457dcb4fece6795e4bca6ae64348f2687c97_1000.fastq
    ./sequencing_summary
        PCT0035_20190226_0004A30B00245CD6_2_A9_D9_sequencing_run_hcc33t_s5_34419_sequencing_summary.txt
    """

    logging.info("Grabbing run information from directory")

    config_list = []

    # Grab fastq and fast5 values
    if 'filename_fastq' in summary_df.columns.tolist() and 'filename_fast5' in summary_df.columns.tolist():
        # only if has been basecalled directly
        unique_files = summary_df[['filename_fastq', 'filename_fast5']].drop_duplicates()
        basecall_direct = True
    elif 'filename' in summary_df.columns.tolist():
        unique_files = summary_df[['filename']]
        basecall_direct = False
    else:
        logging.error("Could not find right filename columns in %s" % ','.join(summary_df.columns.tolist()))
        sys.exit(1)
        

    # Add to config dict
    for row in unique_files.itertuples():
        flowcell, run_id, num = row.filename_fastq.split(".", 1)[0].split("_", 3)
        zfill_num = str(num).zfill(zfill)
        if basecall_direct:
            fast5_pass_file = os.path.join('fast5_pass', row.filename_fast5)
            fast5_pass_file_renamed = os.path.join('fast5_pass',
                                                   '_'.join(map(str, [flowcell, rand_id, 'pass', zfill_num]))
                                                   + ".fast5.gz"
                                                   )
            fast5_fail_file = os.path.join('fast5_fail', row.filename_fast5)
            fast5_fail_file_renamed = os.path.join('fast5_fail',
                                                   '_'.join(map(str, [flowcell, rand_id, 'fail', zfill_num]))
                                                   + ".fast5.gz"
                                                   ) 
        else:
            fast5_file = os.path.join('fast5', row.filename)
            fast5_file_renamed = os.path.join('fast5', 
                                              '_'.join(map(str, [flowcell, rand_id, zfill_num])) 
                                               + ".fast5.gz"
                                              )
        fastq_pass_file = os.path.join('fastq_pass', row.filename_fastq)
        fastq_pass_file_renamed = os.path.join('fastq_pass',
                                               '_'.join(map(str, [flowcell, rand_id, 'pass', zfill_num]))
                                               + ".fastq.gz"
                                               )
        fastq_fail_file = os.path.join('fastq_fail', row.filename_fastq)
        fastq_fail_file_renamed = os.path.join('fastq_fail',
                                               '_'.join(map(str, [flowcell, rand_id, 'fail', zfill_num]))
                                               + ".fastq.gz"
                                               )

        # Add md5sum outputs
        output_md5sum_fast5_pass = os.path.join("fast5_pass", "checksum.fast5.pass.md5")
        output_md5sum_fastq_pass = os.path.join("fastq_pass", "checksum.fastq.pass.md5")
        output_md5sum_fast5_fail = os.path.join("fast5_fail", "checksum.fast5.fail.md5")
        output_md5sum_fastq_fail = os.path.join("fast5_fail", "checksum.fastq.fail.md5")

        config_list.append(pd.Series(data=[zfill_num, rand_id,
                                           fast5_pass_file, fast5_pass_file_renamed,
                                           fast5_fail_file, fast5_fail_file_renamed,
                                           fastq_pass_file, fastq_pass_file_renamed,
                                           fastq_fail_file, fastq_fail_file_renamed,
                                           output_md5sum_fast5_pass, output_md5sum_fastq_pass,
                                           output_md5sum_fast5_fail, output_md5sum_fastq_fail],
                                     index=config_columns))
    return pd.concat(config_list, axis='columns', ignore_index=True, sort=True).transpose()


def output_yaml(yaml_file, dataset):
    logging.info("Printing yaml file to %s" % yaml_file)
    with open(yaml_file, 'w') as file:
        yaml.dump(json.loads(dataset.set_index('zfill_num', drop=False).to_json(orient='index')), file, default_flow_style=False)


def sanitise_fastq_files(fastq_path):
    sanitise_command = ['fix_concatenated_fastqs', '--input=%s' % fastq_path, "--non_recursive"]
    logging.info("Running sanitiser on fastq files")
    sanitise_proc = subprocess.run(sanitise_command, capture_output=True)
    logging.info("Completed sanitation")

    if sanitise_proc.returncode != 0:
        logging.warning("Sanitiser produced non-zero exit code")
        logging.warning("Stderr %s" % sanitise_proc.stderr.decode())


def tidy_summary_df(summary_df, config_df):
    """
    Tidy and add columns to summary dataframe from config dataframe
    :param summary_df: pd.DataFrame
    :param config_df: pd.DataFrame
    :return: pd.DataFrame
    """
    # Adjust values of fastq and fast5 filenames
    # Grab number value to bind to config_df
    summary_df['zfill_num'] = summary_df['filename_fastq'].apply(
        lambda x: str(x.split(".")[0].rsplit("_")[-1]).zfill(zfill))


    # Rename fastq file to path
    summary_df = summary_df.merge(config_df, how='inner',
                                  left_on='zfill_num', right_on='zfill_num',
                                  suffixes=("", "_config"))

    # Rename old columns
    rename_columns = ['filename_fast5', 'filename_fastq', 'run_id']
    rename_columns_dict = {key: "%s_orig" % key for key in rename_columns}
    summary_df.rename(columns=rename_columns_dict, inplace=True)

    # Adjust the files (replace columns)
    summary_df['filename_fast5'] = summary_df.apply(lambda x: x.fast5_pass_file_renamed
                                                                  if x.passes_filtering
                                                                  else x.fast5_fail_file_renamed,
                                                    axis='columns')
    summary_df['filename_fastq'] = summary_df.apply(lambda x: x.fastq_pass_file_renamed
                                                                 if x.passes_filtering
                                                                 else x.fastq_fail_file_renamed,
                                                    axis='columns')

    # Adjust the run_id
    summary_df['run_id'] = summary_df.apply(lambda x: x.rand_id, axis='columns')

    return summary_df


def output_mini_dfs(summary_df, summary_dir, fcid, rand_id, active=False):
    """
    Output mini df, a subset of the bulk df
    :param summary_df: pd.DataFrame
    :param summary_dir: path
    :param fcid: string
    :param rand_id: int
    :param active: bool
    :return: None
    """
    for zfill_num in summary_df['zfill_num'].unique().tolist():
        # Don't play with the last fast5 file if last row is still active
        if active and zfill_num == max(summary_df['zfill_num'].tolist()):
            # Don't play with the last file
            break
        # Grab summary df, drop config columns and any that may have been renamed
        # to the _config suffix during the merge
        mini_summary_df = summary_df.query("zfill_num=='%s'" % zfill_num).\
            drop(columns=config_columns, errors='ignore').\
            filter(regex='.*(?<!_config)$', axis='columns')
        # Set output file
        mini_output_file = os.path.join(summary_dir,
                                        '_'.join(map(str, [fcid, rand_id, zfill_num, "sequencing_summary"])) + ".txt")
        mini_summary_df.to_csv(mini_output_file, sep="\t", index=False, header=True)


def main(args):
    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)

    # First sanitise fastq files
    if args.sanitiser:
        sanitise_fastq_files(os.path.join(args.run_dir, 'fastq_pass'))
        sanitise_fastq_files(os.path.join(args.run_dir, 'fastq_fail'))

    # Split out components of the run to grab
    # Date, UTCTIME, PORT, FCIDD, RID
    args.run_dir = os.path.normpath(args.run_dir)
    try:
        date, utc_time, port, fcid, rand_id = os.path.basename(args.run_dir).split("_", 5)
    except ValueError:
        logger.error("Could not unpack date, utctime port fcid and rand_id from %s"
                     % os.path.basename(args.run_dir))
        sys.exit(1)

    summary_dir = os.path.join(args.run_dir, 'sequencing_summary')
    try:
        sequencing_summary_file = next(
            filter(lambda x:
                   x.endswith("_sequencing_summary.txt")
                   and not x.endswith("_sequencing_summary.bulk.txt")
                   and 'sequencing_run' in x,
                   iter(os.listdir(summary_dir))
                   )
        )
        sequencing_summary_file = os.path.join(summary_dir, sequencing_summary_file)

    except (StopIteration, FileNotFoundError):
        logger.error("Error, could not find summary file in %s" % summary_dir)
        sys.exit(1)

    # Read in summary dataframe
    summary_df = pd.read_csv(sequencing_summary_file, sep='\t', low_memory=False)

    # Grab all files and create new dataframe
    config_df = get_all_files(rand_id, summary_df)

    # Copy sequencing_summary_file to bulk file
    sequencing_summary_file_renamed = re.sub("_sequencing_summary.txt$",
                                             "_sequencing_summary.bulk.txt", sequencing_summary_file)

    # Tidy summary df
    summary_df = tidy_summary_df(summary_df, config_df)

    # Output smaller dfs
    output_mini_dfs(summary_df, summary_dir, fcid, rand_id, active=args.active)

    # Copy across the bulk file
    if not os.path.isfile(sequencing_summary_file_renamed):
        shutil.copy(sequencing_summary_file, sequencing_summary_file_renamed)

    # Output the yaml file
    output_yaml(args.output_yaml_file, config_df)


if __name__ == "__main__":
    main()

