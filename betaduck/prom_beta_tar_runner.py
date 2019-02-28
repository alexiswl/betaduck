#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from datetime import datetime
import subprocess
import shutil
import time
import gzip

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

"""
Usage: Given a list of fast5, fastq and md5 outputs, gzip up fastq and fast5, write md5 to md5sum.
Rename paths to shorten run-id
"""


def get_args():
    parser = argparse.ArgumentParser(description="Compress nanopore data, fast5 to fast5.gz, fastq to fastq.gz")
    parser.add_argument('--fastq-input-pass-path',
                        help="Path to input fastq  pass file", required=True)
    parser.add_argument('--fastq-output-pass-path',
                        help="Path to output fastq pass file", required=True)
    parser.add_argument("--fast5-input-pass-path",
                        help="Path to input fast5 pass file", required=True)
    parser.add_argument("--fast5-output-pass-path",
                        help="Path to output fast5 pass file", required=True)
    parser.add_argument("--md5-pass-fast5",
                        help="File to append to md5sum for pass fast5 files", required=True)
    parser.add_argument("--md5-pass-fastq",
                        help="File to append to md5sum for pass fastq files", required=True)
    parser.add_argument('--fastq-input-fail-path',
                        help="Path to input fastq  fail file", required=True)
    parser.add_argument('--fastq-output-fail-path',
                        help="Path to output fastq fail file", required=True)
    parser.add_argument("--fast5-input-fail-path",
                        help="Path to input fast5 fail file", required=True)
    parser.add_argument("--fast5-output-fail-path",
                        help="Path to output fast5 fail file", required=True)
    parser.add_argument("--md5-fail-fast5",
                        help="File to append to md5sum for fail fast5 files", required=True)
    parser.add_argument("--md5-fail-fastq",
                        help="File to append to md5sum for fail fastq files", required=True)
    parser.add_argument("--inplace", action='store_true', default=False, help='Remove uncompressed files')
    parser.add_argument("--overwrite", action='store_true', default=False,
                        help="Overwrite output file rather than append to it")
    parser.add_argument("--dry-run", action='store_true', default=False,
                        help="Don't actually zip up anything, just output the logs")
    args = parser.parse_args()
    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)
    return args


def zip_and_move_file(file_path, file_output_path, overwrite=False, inplace=False, dry_run=False):
    if not dry_run:
        if os.path.isfile(file_output_path) and not overwrite:
            logging.info("Fastq file %s already exists in destination and overwrite not set. "
                         "Skipping" % file_output_path)

        # Zip file to .tmp file and then move to .gz
        tmp_output_path = file_output_path + ".tmp"

        # Zip and move the summary file
        with open(file_path, 'rb') as f_in, gzip.open(tmp_output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Now move to final dest, wait for filesystem to catch up first
        time.sleep(1)
        shutil.move(tmp_output_path, file_output_path)

        if inplace:
            # Wait for file system to catch up then remove
            time.sleep(1)
            os.remove(file_path)
    else:
        logging.info("Would have gzipped and moved %s into %s" % (file_path, file_output_path))


def get_md5sum(output_path):
    logging.info("Obtaining the md5sum for %s" % output_path)

    # Grab the md5sum of the file. Use the relative path
    output_file = os.path.basename(os.path.normpath(output_path))
    work_dir = os.path.dirname(os.path.normpath(output_path))
    md5_command = ['md5sum', output_file]
    md5_proc = subprocess.run(md5_command, capture_output=True, cwd=work_dir)
    md5_output = md5_proc.stdout.decode().splitlines()[0]
    logging.info("Obtained %s as md5 for %s" % (md5_output, output_path))

    # Return the md5 of the file for writing to a checksum file
    return md5_output


def write_md5sum(md5_sum, output_md5_file):
    logging.info("Writing md5sum '%s' to %s" % (md5_sum, output_md5_file))
    with open(output_md5_file, 'a') as md5_h:
        md5_h.write(md5_sum + "\n")


def main():
    # Get args
    args = get_args()

    # Log arguments
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)

    # Gzip up fast5
    zip_and_move_file(args.fast5_input_pass_path, args.fast5_output_pass_path,
                      overwrite=args.overwrite, inplace=args.inplace, dry_run=args.dry_run)
    zip_and_move_file(args.fast5_input_fail_path, args.fast5_output_fail_path,
                      overwrite=args.overwrite, inplace=args.inplace, dry_run=args.dry_run)

    # Move fastq folder
    zip_and_move_file(args.fastq_input_pass_path, args.fastq_output_pass_path,
                      overwrite=args.overwrite, inplace=args.inplace, dry_run=args.dry_run)
    zip_and_move_file(args.fastq_input_fail_path, args.fastq_output_fail_path,
                      overwrite=args.overwrite, inplace=args.inplace, dry_run=args.dry_run)

    # Get md5 for fastq and fast5
    if not args.dry_run:
        # Get md5
        md5sum_fast5_pass = get_md5sum(args.fast5_output_pass_path)
        md5sum_fast5_fail = get_md5sum(args.fast5_output_fail_path)
        md5sum_fastq_pass = get_md5sum(args.fastq_output_pass_path)
        md5sum_fastq_fail = get_md5sum(args.fastq_output_fail_path)

        # Append md5 to file
        write_md5sum(md5sum_fast5_pass, args.md5_pass_fast5)
        write_md5sum(md5sum_fast5_fail, args.md5_fail_fast5)
        write_md5sum(md5sum_fastq_pass, args.md5_pass_fastq)
        write_md5sum(md5sum_fastq_fail, args.md5_fail_fastq)


if __name__ == "__main__":
    main()
