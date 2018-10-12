#!/usr/bin/env python3

import os
import sys
import argparse
import betaduck.version


# Return version number
def get_version():
    # Version should correspond to MinKNOW version
    return betaduck.__version__


# Run the function we select through argparse
def run_function(args):
    # Which function of the three did we choose?
    if args.command == "config":
        import betaduck.prom_beta_tar_gen_config as command_to_run
    if args.command == "tar":
        import betaduck.prom_beta_tar_wrapper as command_to_run
    if args.command == "plot":
        import betaduck.prom_beta_plotter_wrapper as command_to_run

    # Now run it!
    command_to_run.main(args)


# Define main script:
def main():
    # Create betaduck parser
    parser = argparse.ArgumentParser(prog='betaduck', description="betaduck package")
    parser.add_argument("--version", help="Get version of betaduck",
                        action="version",
                        version=betaduck.__version__)

    subparsers = parser.add_subparsers(help="Callable betaduck functions", dest="command")

    # Generate config arguments
    config_parser = subparsers.add_parser('config',
                                          help="Generate a config file that will be used to organise folder")
    config_parser.add_argument('--sequencing_summary_path',
                               help="Path to sequencing_summary_path", required=True)
    config_parser.add_argument('--fastq_path',
                               help="Path to fastq files", required=True)
    config_parser.add_argument("--fast5_path",
                               help="Path to folder of fast5 files", required=True)
    config_parser.add_argument("--output_yaml_file",
                               help="Yaml file to create", required=True)
    config_parser.add_argument("--sanitiser", action='store_true', default=False, 
                               help="run fastq sanitiser before generating config")
    config_parser.add_argument("--active", action='store_true', default=False,
                               help="Don't tar up the last folder as data may still be writing to there")
    config_parser.set_defaults(func=run_function)

    # Tar command
    tar_parser = subparsers.add_parser('tar',
                                       help="Use the config file to now tar up the directories")
    tar_parser.add_argument("--config", required=True,
                            help="Path to config file")
    tar_parser.add_argument("--keep", default=False, action='store_true',
                            help="Keep the original files")
    tar_parser.add_argument("--dry_run", default=False, action='store_true',
                            help="Just log the commands, do not run them")
    tar_parser.add_argument("--overwrite", default=False, action='store_true',
                            help="Overwrite files if they already exist")
    tar_parser.add_argument("--threads", default=1, type=int,
                            help="Number of types ")
    tar_parser.set_defaults(func=run_function)

    # Plotter
    plotter_parser = subparsers.add_parser('plot',
                                           help="Plot the run(s)")
    plotter_parser.add_argument("--summary_dir", type=str, required=True,
                                help="Contains the txt files (comma separated for multiple locations)")
    plotter_parser.add_argument("--fastq_dir", type=str, required=True,
                                help="Where are the fastq files (comma separated for multiple locations)")
    plotter_parser.add_argument("--plots_dir", type=str, required=True,
                                help="Where do the plots go")
    plotter_parser.add_argument("--name", type=str, required=True,
                                help="Titles for plots")
    plotter_parser.add_argument("--threads", type=str, required=True,
                                help="Read the dataframes in parallel")
    plotter_parser.set_defaults(func=run_function)

    args = parser.parse_args()

    # Print help if just 'betaduck' is typed in.
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    main()
