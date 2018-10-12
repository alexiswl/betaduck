#!/usr/bin/env python3

import argparse
import yaml
import json
import pandas as pd
import subprocess
import logging
import os
import concurrent.futures
import numpy as np
import time

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
    time.sleep(2)

    if tar_proc.returncode == 0:
        logging.info("Process completed successfully")
        logging.info("Stdout = %s" % tar_proc.stdout.decode())
        logging.info("Stderr = %s" % tar_proc.stderr.decode())
        return True
    else:
        logging.warning("Process returned non-zero exit code.")
        logging.warning("Stdout = %s" % tar_proc.stdout.decode())
        logging.warning("Stderr = %s" % tar_proc.stderr.decode())
        return False


def main(args):
    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)

    # Read in pandas dataframe
    dataframe = pd.DataFrame(read_config(args.config))

    # Reduce thread count unless already 1.
    threads = 1 if args.threads == 1 else args.threads - 1
    logging.info("Given we need to take of the parent script, running %d jobs in parallel" % threads)
   
    # Run in parallel 
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        iterator = {executor.submit(run_process, config, args.keep, args.overwrite, args.dry_run): 
                    config for config in dataframe.itertuples()}
        for item in concurrent.futures.as_completed(iterator):
            pandas_input = iterator[item]
            success = item.result()
            

if __name__ == "__main__":
    main()
