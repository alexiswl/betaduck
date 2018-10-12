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
    config_parser.set_defaults(func=run_function)

    # Tar command
    tar_parser = subparsers.add_parser('tar',
                                       help="Use the config file to now tar up the directories")
    tar_parser.add_argument("--config", required=True,
                            help="Path to config file")
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
    plotter_parser.set_defaults(func=run_function)

    args = parser.parse_args()

    # Print help if just 'betaduck' is typed in.
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    main()
