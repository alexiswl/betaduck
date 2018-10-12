#!/usr/bin/env python3

import argparse
import yaml
import json
import pandas as pd
import subprocess
import logging
import os

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def read_config(config):
    with open(config) as f:
        config_data = yaml.load(f)
    return config_data


def run_process(config_data, keep=False, overwrite=False, dry_run=False):
    # Get path to this github repo
    here = os.path.dirname(os.path.realpath(__file__))

    # Generate tar command
    tar_command = ["python", os.path.join(here, "prom_beta_tar_runner.py"),
                   # Then come the options.
                   "--sequencing_summary_path=%s" % config_data.sequencing_summary_file,
                   "--fastq_path=%s" % config_data.fastq_file,
                   "--fast5_path=%s" % config_data.fast5_dir,
                   "--flowcellID=%s" % config_data.FlowcellID,
                   "--rnumber=%s" % config_data.rnumber,
                   "--md5_fast5=%s" % config_data.md5_fast5,
                   "--md5_fastq=%s" % config_data.md5_fastq]
    # Do we want to keep the data
    if not keep:
        tar_command.append("--inplace")

    # Do we want to test
    if dry_run:
        tar_command.append("--dry-run")

    # Do we want to overwrite
    if overwrite:
        tar_command.append("--overwrite")

    tar_proc = subprocess.run(tar_command, capture_output=True)

    if tar_proc.returncode == 0:
        logging.info("Process completed successfully")
        logging.info("Stdout = %s" % tar_proc.stdout.decode())
        logging.info("Stderr = %s" % tar_proc.stderr.decode())
    else:
        logging.warning("Process returned non-zero exit code.")
        logging.warning("Stdout = %s" % tar_proc.stdout.decode())
        logging.warning("Stderr = %s" % tar_proc.stderr.decode())


def main(args):
    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)

    # Read in pandas dataframe
    dataframe = pd.DataFrame(read_config(args.config))

    # Iterate through each row of the configuration file.
    for row in dataframe.itertuples():
        run_process(row, keep=args.keep, overwrite=args.overwrite, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
