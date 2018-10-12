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

# """
# Generate a yaml file used to run each of the alpha_light python commands
# """

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def get_flowcell_id(fast5_file):
    with h5py.File(fast5_file) as f:
        # Get flowcell ID from UniqueGlobalKey/tracking_id
        try:
            read_attributes = dict(f['UniqueGlobalKey/tracking_id'].attrs.items())
        except KeyError:
            logging.warning("Could not find flowcell ID from %s" % fast5_file)
            return None

        return read_attributes['flow_cell_id']


def get_random_number(fast5_file):
    _, rnumber, _, read, _, channel, _ = fast5_file.rsplit("_", 6)
    try:
        rnumber = int(rnumber)
        return rnumber
    except TypeError:
        logging.warning("Tried to get rnumber from fast5 file. Got %s from %s" % (fast5_file, rnumber))
        return None


def get_all_files(sequencing_summary_dir, fastq_dir, fast5_dir):
    """
    :rtype: pd.DataFrame
    :param sequencing_summary_dir: string
    :param fastq_dir: string
    :param fast5_dir: string
    """

    logging.info("Grabbing sequencing summary files")
    sequencing_summary_files = [os.path.join(sequencing_summary_dir, sequencing_summary_file)
                                for sequencing_summary_file in os.listdir(sequencing_summary_dir)
                                if re.match('sequencing_summary_\d+.txt', sequencing_summary_file)]

    logging.info("Grabbing fastq files")
    fastq_files = [os.path.join(fastq_dir, fastq_file)
                   for fastq_file in os.listdir(sequencing_summary_dir)
                   if re.match('fastq_\d+.fastq', fastq_file)]

    logging.info("Grabbig fast5 directories")
    fast5_dirs = [os.path.join(fast5_dir, fast5_folder)
                  for fast5_folder in os.listdir(fast5_dir)
                  if os.path.isdir(os.path.join(fast5_dir, fast5_folder))
                  and re.match("^\d+$", fast5_folder)]

    # Get rnumber and flowcell id
    logging.info("Grabbing a flowcell ID from the fast5 attributes")
    flowcell_id = None
    rnumber = None
    while flowcell_id is None or rnumber is None:
        for fast5_dir in fast5_dirs:
            fast5_files = [os.path.join(fast5_dir, fast5_file) 
                           for fast5_file in os.listdir(fast5_dir)
                           if fast5_file.endswith('.fast5')
                           and os.path.isfile(os.path.join(fast5_dir, fast5_file))]
            for fast5_file in fast5_files:
                flowcell_id = get_flowcell_id(fast5_file)
                rnumber = get_random_number(fast5_file)
                if rnumber is not None and flowcell_id is not None:
                    break
        if rnumber is not None and flowcell_id is not None:
            break

    logging.info("Got flowcell ID as %s" % flowcell_id) 
    logging.info("Got rnumber as %s" % rnumber)

    sequencing_summary_df = pd.DataFrame(sequencing_summary_files, 
                                         columns=["sequencing_summary_file"])
    fastq_df = pd.DataFrame(fastq_files, columns=["fastq_file"])
    fast5_df = pd.DataFrame(fast5_dirs, columns=["fast5_dir"])

    # Append number onto each dataframe
    sequencing_summary_df['number'] = sequencing_summary_df['sequencing_summary_file'].apply(
        lambda x: int(re.match("sequencing_summary_(\d+).txt", os.path.basename(x)).group(1)))
    fastq_df['number'] = fastq_df['fastq_file'].apply(
        lambda x: int(re.match("fastq_(\d+).fastq", os.path.basename(x)).group(1)))
    fast5_df['number'] = fast5_df['fast5_dir'].apply(
        lambda x: int(re.match('(\d+)', os.path.basename(x)).group(1)))

    # Sort dataframes by number
    sequencing_summary_df.sort_values(by=['number'], inplace=True)
    fastq_df.sort_values(by=['number'], inplace=True)
    fast5_df.sort_values(by=['number'], inplace=True)

    # Set number as index for each dataframe
    sequencing_summary_df.set_index("number", inplace=True)
    fastq_df.set_index("number", inplace=True)
    fast5_df.set_index("number", inplace=True)

    # Now merge each using pd.concat
    dataset = pd.concat([sequencing_summary_df, fastq_df, fast5_df], axis='columns', join='inner', sort=True)

    # Add on flowcellIDs and rnumbers
    dataset['FlowcellID'] = flowcell_id
    dataset['rnumber'] = rnumber

    return dataset


def output_yaml(yaml_file, dataset):
    logging.info("Printing yaml file to %s" % yaml_file)
    with open(yaml_file, 'w') as file:
        yaml.dump(json.loads(dataset.to_json(orient='records')), file, default_flow_style=False)


def sanitise_fastq_files(fastq_path):
    sanitise_command = ['fix_concatenated_fastqs', '--input=%s' % fastq_path, "--non_recursive"]
    logging.info("Running sanitiser on fastq files")
    sanitise_proc = subprocess.run(sanitise_command, capture_output=True)
    logging.info("Completed sanitation")

    if sanitise_proc.returncode != 0:
        logging.warning("Sanitiser produced non-zero exit code")
        logging.warning("Stderr %s" % sanitise_proc.stderr.decode())


def main(args):
    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)

    # First sanitise fastq files
    if args.sanitiser:
        sanitise_fastq_files(args.fastq_path)

    # Read in dataset
    dataset = get_all_files(sequencing_summary_dir=os.path.abspath(args.sequencing_summary_path),
                            fastq_dir=os.path.abspath(args.fastq_path),
                            fast5_dir=os.path.abspath(args.fast5_path))
    logging.info("Obtained %d folders" % dataset.shape[0])

    # Add md5sum outputs
    output_md5sum_fast5 = os.path.join(os.path.abspath(args.fast5_path), "checksum.fast5.md5")
    output_md5sum_fastq = os.path.join(os.path.abspath(args.fastq_path), "checksum.fastq.md5")
    dataset['md5_fast5'] = output_md5sum_fast5
    dataset['md5_fastq'] = output_md5sum_fastq

    # Output the yaml file
    output_yaml(args.output_yaml_file, dataset)


if __name__ == "__main__":
    main()
