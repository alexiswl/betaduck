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
    elif args.command == "tidy":
        import betaduck.prom_beta_tar_wrapper as command_to_run
    elif args.command == "plot":
        import betaduck.prom_beta_plotter_wrapper as command_to_run
    if args.command == "align":
        import betaduck.prom_beta_align_wrapper as command_to_run
    else:
        return

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
    config_parser.add_argument('--run-dir',
                               help="Path to run", required=True)
    config_parser.add_argument("--output-yaml-file",
                               help="Yaml file to create", required=True)
    config_parser.add_argument("--sanitiser", action='store_true', default=False, 
                               help="run fastq sanitiser before generating config")
    config_parser.add_argument("--active", action='store_true', default=False,
                               help="Don't tar up the last folder as data may still be writing to there")
    config_parser.set_defaults(func=run_function)

    # Tar command
    tar_parser = subparsers.add_parser('tidy',
                                       help="Use the config file to now tidy and tar up the directories")
    tar_parser.add_argument("--config", required=True,
                            help="Path to config file")
    tar_parser.add_argument("--keep", default=False, action='store_true',
                            help="Keep the original files")
    tar_parser.add_argument("--dry-run", default=False, action='store_true',
                            help="Just log the commands, do not run them")
    tar_parser.add_argument("--overwrite", default=False, action='store_true',
                            help="Overwrite files if they already exist")
    tar_parser.add_argument("--threads", default=1, type=int,
                            help="Number of folders to zip up simultaneously")
    tar_parser.set_defaults(func=run_function)

    # Plotter
    plotter_parser = subparsers.add_parser('plot',
                                           help="Plot the run(s)")
    plotter_parser.add_argument("--summary-dir", type=str, required=True,
                                help="Contains the txt files (comma separated for multiple locations)")
    plotter_parser.add_argument("--fastq-dir", type=str, required=True,
                                help="Where are the fastq files (comma separated for multiple locations)")
    plotter_parser.add_argument("--plots-dir", type=str, required=True,
                                help="Where do the plots go")
    plotter_parser.add_argument("--name", type=str, required=True,
                                help="Titles for plots")
    plotter_parser.add_argument("--threads", type=int, default=1,
                                help="Read the dataframes in parallel")
    plotter_parser.set_defaults(func=run_function)

    # Args for alignment
    align_parser = subparsers.add_parser('align',
                                         help="Align reads to a reference genome")
    align_parser.add_argument("--fastq_dir", type=str, required=True,
                              help="Path to gzipped fastq files")
    align_parser.add_argument("--output_dir", type=str, required=True,
                              help="Folder to place bam files")
    align_parser.add_argument("--genome_dir", type=str, required=True,
                              help="Path to directory of genomes")
    align_parser.add_argument("--genome", type=str, required=True,
                              help="Name of genome in genome directory")
    align_parser.add_argument("--filter_lambda", dest="w_lambda",
                              default=False, action='store_true')
    align_parser.add_argument("--cs", type=str, default='long',
			      choices=["long", "short", "none"],
			      help="CS tag type")
    align_parser.add_argument("--md", default=False, action='store_true',
			      help="Add MD tag to BAM file")
    align_parser.add_argument("--threads", type=int, default=1)
    align_parser.set_defaults(func=run_function)

    args = parser.parse_args()

    # Print help if just 'betaduck' is typed in.
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    main()
