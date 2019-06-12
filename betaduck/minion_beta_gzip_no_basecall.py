#!/usr/bin/env python3

import argparse
import os
import psutil
import logging
import subprocess
import sys
import concurrent.futures

# Set logging
logging.basicConfig(level=logging.INFO)

def get_args():
    # Tar command
    tar_parser = argparse.ArgumentParser(description="Mini Script to wrap up MinION data")
    tar_parser.add_argument("--run-dir", required=True,
                            help="Path to run directory")
    tar_parser.add_argument("--keep", default=False, action='store_true',
                            help="Keep the original files")
    tar_parser.add_argument("--dry-run", default=False, action='store_true',
                            help="Just log the commands, do not run them")
    tar_parser.add_argument("--overwrite", default=False, action='store_true',
                            help="Overwrite files if they already exist")
    tar_parser.add_argument("--threads", default=1, type=int,
                            help="Number of folders to zip up simultaneously")

    args = tar_parser.parse_args()

    return args


def set_args(args):
    """

    :param args:
    :return:
    """

    # Make sure run directory is accurate
    setattr(args, "run_dir", os.path.normpath(args.run_dir))

    return args


def is_open(fpath):
    """
    Uses the psutil function to determine if a file is open.
    Useful for determining if the fast5 file has been completely written to
    :param filename:
    :return:
    """
    for proc in psutil.process_iter():
        try:
            for item in proc.open_files():
                if fpath == item.path:
                    return True
        except Exception:
            pass
    return False


def check_overwrite(fpath, overwrite):
    """
    Check that the gzipped suffix doesn't already exist
    :param fpath:
    :param overwrite:
    :return:
    """
    if overwrite:
        return False
    else:
        if os.path.isfile(fpath + ".gz"):
            return True


def get_fast5_files(fast5_directory, overwrite=True):
    """
    :param fast5_directory: /var/lib/MinKNOW/data/RUN_ID/SAMPLE_ID/RUN_DIR/fast5
    :return:
    """
    fast5_files = [fast5_file
                   for fast5_file in os.listdir(fast5_directory)
                   if os.path.isfile(os.path.join(fast5_directory, fast5_file))
                   and fast5_file.endswith(".fast5")
                   and not is_open(fast5_file)
                   and not check_overwrite(os.path.join(fast5_directory, fast5_file), overwrite)]

    return fast5_files


def compress_fast5_file(fast5_file, fast5_directory, debug=False):
    """
    Simpler command, compress fastq file using gzip command
    :param fastq_file:
    :return:
    """

    # Generate command
    gzip_command = ['gzip', '--best', fast5_file]

    logging.info("Gzipping file %s with best compression" % fast5_file)
    logging.info("Command is %s" % gzip_command)

    if debug:
        return
    # Run through subprocess
    gzip_proc = subprocess.run(gzip_command, cwd=fast5_directory, capture_output=True)

    if not gzip_proc.returncode == 0:
        logging.error("Gzip command '%s' returned non-zero exit code" % ' '.join(gzip_command))
        sys.exit(1)
    else:
        logging.info("Gzip of %s was successful" % fast5_file)


def compress_all_fast5_files(input_fast5_files, fast5_directory, threads=1, debug=False):
    # Run simultaneously
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        iterator = {executor.submit(compress_fast5_file, fast5_file, fast5_directory, debug=debug)
                    for fast5_file in input_fast5_files}
        for item in concurrent.futures.as_completed(iterator):
            #output = iterator[item]
            success = item.result()


def main():
    # Args
    args = get_args()
    args = set_args(args)
    fast5_directory = os.path.join(args.run_dir, 'fast5')

    # Get files
    fast5_files = get_fast5_files(fast5_directory)

    # Zip files (through concurrent futures tool).
    compress_all_fast5_files(fast5_files, fast5_directory, threads=args.threads, debug=args.dry_run)


if __name__=="__main__":
    main()
