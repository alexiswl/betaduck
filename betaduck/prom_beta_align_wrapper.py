#!/usr/bin/env python3

import os
import subprocess
import logging
import sys
import re
import concurrent.futures
import shutil
from Bio import SeqIO
import wub.util.misc


class Genome:
    def __init__(self, genome_dir, genome, w_lambda=False):
        self.path = os.path.join(genome_dir, genome)
        self.name = genome
        self.host_genome_path = os.path.join(genome_dir, genome, "genome.fa")
        self.host_minimap2_index = os.path.join(genome_dir, genome, "genome.mmi")
        self.w_lambda = False

        # Change a few settings if we want to have lambda filtering in the analysis.
        if w_lambda:
            self.w_lambda = True
            self.host_w_lambda_genome_path = os.path.join(genome_dir, genome, "genome.w_lambda.fa")
            self.host_w_lambda_minimap2_index = os.path.join(genome_dir, genome, "genome.w_lambda.mmi")
            self.lambda_path = os.path.join(genome_dir, "lambda")
            self.lambda_genome_path = os.path.join(genome_dir, "lambda", "genome.fa")
            self.lambda_genome_bed = os.path.join(genome_dir, "lambda", "genome.bed")
            self.lambda_samtools_index = os.path.join(genome_dir, "lambda", "genome.fa.fai")


class Sample:
    def __init__(self, input_dir, output_dir, genome):
        self.fastq_files = collect_fastqs(input_dir)
        self.genome = genome
        self.w_lambda = genome.w_lambda
        self.alignment_objects = [SubFolder(fastq_file, output_dir, genome, w_lambda=self.w_lambda)
                                  for fastq_file in self.fastq_files]
        _, self.flowcell, self.rnumber = re.sub(".fastq.gz", "",
                                                os.path.basename(os.path.normpath(self.fastq_files[0]))).split("_", 3)
        self.sample_prefix = '_'.join([genome.name, self.flowcell, self.rnumber])
        self.unaligned_merged_bam_file = os.path.join(output_dir, 'merged',
                                                      '_'.join(['unaligned', self.flowcell, self.rnumber])
                                                      + ".sorted.merged.bam")
        self.unaligned_merged_pickle_file = os.path.join(output_dir, 'merged',
                                                         '_'.join(['unaligned', self.flowcell, self.rnumber])
                                                         + ".merged.pickle")
        self.unaligned_files = [align.unaligned for align in self.alignment_objects]
        self.unaligned_pickles = [align.unaligned_pickle for align in self.alignment_objects]

        if self.w_lambda:
            self.host_merged_bam_file = os.path.join(output_dir, 'merged',
                                                     '_'.join([genome.name+'.lambda-filt', self.flowcell, self.rnumber])
                                                     + ".sorted.merged.bam")
            self.lambda_merged_bam_file = os.path.join(output_dir, 'merged',
                                                       '_'.join(['lambda', self.flowcell, self.rnumber])
                                                       + ".sorted.merged.bam")
            self.host_merged_pickle_file = os.path.join(output_dir, 'merged',
                                                        '_'.join([genome.name+'.lambda-filt', self.flowcell, self.rnumber])
                                                        + ".merged.pickle")
            self.lambda_merged_pickle_file = os.path.join(output_dir, 'merged',
                                                          '_'.join(['lambda', self.flowcell, self.rnumber])
                                                          + ".merged.pickle")
            self.host_alignment_files = [align.host_aligned for align in self.alignment_objects]
            self.lambda_alignment_files = [align.lambda_aligned for align in self.alignment_objects]
            self.host_alignment_pickles = [align.host_pickle for align in self.alignment_objects]
            self.lambda_alignment_pickles = [align.lambda_pickle for align in self.alignment_objects]

            self.merged_bam_files = [self.host_merged_bam_file, self.unaligned_merged_bam_file,
                                     self.lambda_merged_bam_file]
            self.merged_reference_names = [self.genome.name, self.genome.name, "lambda"]
            self.merged_references = [genome.host_genome_path, genome.host_genome_path, genome.lambda_genome_path]
            self.merged_pickles = [self.host_merged_pickle_file, self.unaligned_merged_pickle_file,
                                   self.lambda_merged_pickle_file]
        else:
            self.host_merged_bam_file = os.path.join(output_dir, 'merged',
                                                     '_'.join([genome.name, self.flowcell, self.rnumber])
                                                     + ".sorted.merged.bam")
            self.host_merged_pickle_file = os.path.join(output_dir, 'merged',
                                                        '_'.join([genome.name, self.flowcell, self.rnumber])
                                                        + ".merged.pickle")
            self.host_alignment_files = [align.host_aligned for align in self.alignment_objects]
            self.merged_bam_files = [self.host_merged_bam_file, self.unaligned_merged_bam_file]
            self.merged_reference_names = [self.genome.name, self.genome.name]
            self.merged_references = [genome.host_genome_path, genome.host_genome_path]
            self.merged_pickles = [self.host_merged_pickle_file, self.unaligned_merged_pickle_file]


class SubFolder:
    def __init__(self, fastq_file, output_dir, genome, w_lambda=False):
        self.fastq_path = fastq_file
        self.prefix = re.sub(".fastq.gz$", "", os.path.basename(os.path.normpath(fastq_file)))
        self.w_lambda = w_lambda
        self.genome = genome
        if w_lambda:
            self.index = genome.host_w_lambda_minimap2_index
            self.host_aligned = os.path.join(output_dir, genome.name, self.prefix+".lambda-filt.sorted.bam")
            self.lambda_aligned = os.path.join(output_dir, 'lambda', self.prefix+".lambda.sorted.bam")
            self.host_pickle = os.path.join(output_dir, 'wub', self.prefix+".lambda-filt.pickle")
            self.lambda_pickle = os.path.join(output_dir, 'wub', self.prefix+".lambda.pickle")
            self.host_report = os.path.join(output_dir, 'wub', self.prefix+".lambda-filt.pdf")
            self.lambda_report = os.path.join(output_dir, 'wub', self.prefix+".lambda.pdf")
        else:
            self.host_pickle = os.path.join(output_dir, 'wub', self.prefix+".pickle")
            self.host_report = os.path.join(output_dir, 'wub', self.prefix+".pdf")
            self.index = genome.host_minimap2_index
            self.host_aligned = os.path.join(output_dir, genome.name, self.prefix + ".sorted.bam")
            self.lambda_aligned = None

        self.unaligned = os.path.join(output_dir, "unaligned", self.prefix+".unaligned.bam")
        self.unaligned_pickle = os.path.join(output_dir, 'wub', self.prefix+".unaligned.pickle")
        self.unaligned_report = os.path.join(output_dir, 'wub', self.prefix+".unaligned.pdf")


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M')

logger = logging.getLogger()

lambda_genome = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/GCF_000840245.1_ViralProj14204_genomic.fna.gz"

# Run minimap2 through a subprocess command.
# Import threads to determine how many to run at once.
# Generates data to directory 'alignment'
# Generates log documenting
# First need to generate genome index
# Then need to run through minimap2
# Create lambda index

# Example minimap2 index
# minimap2 -x map-ont -d /path/to/index /path/to/genome
# Example minimap2 -a /path/to/index /path/to/reads


def generate_index(index_file, genome_file):
    # Create index command
    minimap2_index_command = ["minimap2", "-x", "map-ont",
                              "-d", index_file, genome_file]
    # Run command
    # noinspection PyArgumentList
    minimap2_index_proc = subprocess.run(minimap2_index_command, capture_output=True)

    # Check output of command
    if not minimap2_index_proc.returncode == 0:
        logger.warning("Index command returned non-zero exit code")
        logger.warning("Stderr: %s" % minimap2_index_proc.stderr.decode())


def generate_alignment(subfolder, md, cs):
    # Create alignment command
    minimap2_alignment_command = ["minimap2", "-a"]
    # Add in parameters
    if md:
        minimap2_alignment_command.append("--MD")
    if cs == "long":
        minimap2_alignment_command.append("--cs=long")
    elif cs == "short":
        minimap2_alignment_command.append("--cs=short")
    # Add index and alignment output
    minimap2_alignment_command.extend([subfolder.index, subfolder.fastq_path])
    logging.info("Running command %s minimap2_alignment_command" % minimap2_alignment_command)
    logging.info("Writing to %s" % subfolder.host_aligned)
    logging.info("Completed command %s" % minimap2_alignment_command)
    # Write to alignment_file
    # Write to sorted bam file and sort bam flie
    dev_null = open(os.devnull, 'w')
    if not os.path.isdir(os.path.dirname(subfolder.unaligned)):
        os.mkdir(os.path.dirname(subfolder.unaligned))

    p1 = subprocess.Popen(minimap2_alignment_command,
                          stdout=subprocess.PIPE, stderr=dev_null)
    p2 = subprocess.Popen(["samtools", "view", '-bu', '-F4', '-U', subfolder.unaligned],
                          stdin=p1.stdout,
                          stdout=subprocess.PIPE, stderr=dev_null)
    p1.stdout.close()
    if not subfolder.w_lambda:
        p3 = subprocess.Popen(["samtools", "sort", "-o", subfolder.host_aligned],
                              stdin=p2.stdout,
                              stdout=dev_null, stderr=dev_null)
        p2.stdout.close()
        out, err = p3.communicate()
        if not p3.returncode == 0:
            logging.info("BAM Alignment returned non-zero exit code for alignment file %s" % subfolder.host_aligned)
        # Generate index for alignment file
        index_command = ["samtools", "index", subfolder.host_aligned]
        # noinspection PyArgumentList
        index_proc = subprocess.run(index_command, capture_output=True)
        if not index_proc.returncode == 0:
            logging.warning("Indexing of '%s' returned non-zero exit code" % ' '.join(index_command))
            logging.warning("Stderr: %s" % index_proc.stderr.decode())
    else:
        # If w_lambda, pipe sort to stdout
        p3 = subprocess.Popen(["samtools", "sort"],
                              stdin=p2.stdout,
                              stdout=subprocess.PIPE, stderr=dev_null)
        p2.stdout.close()
        p4 = subprocess.Popen(["samtools", "view", "-b", "-L", subfolder.genome.lambda_genome_bed,
                               "-U", subfolder.host_aligned, "-o", subfolder.lambda_aligned],
                              stdin=p3.stdout,
                              stdout=dev_null, stderr=dev_null)
        p3.stdout.close()
        out, err = p4.communicate()
        if not p4.returncode == 0:
            logging.warning("BAM Alignment returned non-zero exit code for alignment file %s" % subfolder.host_aligned)
        # Generate index for alignment file(s)
        logging.info("Generating indexes for %s and %s" % (subfolder.lambda_aligned, subfolder.host_aligned))
        for alignment_file in [subfolder.lambda_aligned, subfolder.host_aligned]:
            index_command = ["samtools", "index", alignment_file]
            # noinspection PyArgumentList
            index_proc = subprocess.run(index_command, capture_output=True)
            if not index_proc.returncode == 0:
                logging.warning("Indexing of '%s' returned non-zero exit code" % ' '.join(index_command))
                logging.warning("Stderr: %s" % index_proc.stderr.decode())


# Generate pickle for wub
def generate_pickle(subfolder):
    # Host
    bam_alignment_qc_command = ["bam_alignment_qc.py", "-x",
                                "-f", subfolder.genome.host_genome_path,
                                '-p', subfolder.host_pickle,
                                "-r", subfolder.host_report,
                                '-t', subfolder.genome.name,
                                "-Q", subfolder.host_aligned]
    logging.info("Performing bam alignment qc: %s" % ' '.join(bam_alignment_qc_command))
    # noinspection PyArgumentList
    bam_alignment_qc_proc = subprocess.run(bam_alignment_qc_command, capture_output=True)
    if not bam_alignment_qc_proc.returncode == 0:
        logging.warning("Bam QC returned non-zero exit code: %s" % bam_alignment_qc_command)
        logging.warning("Stderr: %s" % bam_alignment_qc_proc.stderr.decode())

    # Lambda
    if subfolder.w_lambda:
        bam_alignment_qc_command = ["bam_alignment_qc.py", "-x",
                                    "-f", subfolder.genome.lambda_genome_path,
                                    '-p', subfolder.lambda_pickle,
                                    "-r", subfolder.lambda_report,
                                    '-t', "lambda",
                                    "-Q", subfolder.lambda_aligned]
        logging.info("Performing bam alignment qc: %s" % ' '.join(bam_alignment_qc_command))
        # noinspection PyArgumentList
        bam_alignment_qc_proc = subprocess.run(bam_alignment_qc_command, capture_output=True)
        if not bam_alignment_qc_proc.returncode == 0:
            logging.warning("Bam QC returned non-zero exit code: %s" % bam_alignment_qc_command)
            logging.warning("Stderr: %s" % bam_alignment_qc_proc.stderr.decode())

    # Unaligned
    bam_alignment_qc_command = ["bam_alignment_qc.py", "-x",
                                "-f", subfolder.genome.host_genome_path,
                                '-p', subfolder.unaligned_pickle,
                                "-r", subfolder.unaligned_report,
                                '-t', subfolder.genome.name,
                                "-Q", subfolder.unaligned]
    logging.info("Performing bam alignment qc: %s" % ' '.join(bam_alignment_qc_command))
    # noinspection PyArgumentList
    #bam_alignment_qc_proc = subprocess.run(bam_alignment_qc_command, capture_output=True)
    #if not bam_alignment_qc_proc.returncode == 0:
    #    logging.warning("Bam QC returned non-zero exit code: %s" % bam_alignment_qc_command)
    #    logging.warning("Stderr: %s" % bam_alignment_qc_proc.stderr.decode())


# Check arguments
def check_args(args):
    # Make sure that the output_directory exists
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)


def collect_fastqs(fastq_dir):
    fastq_files = sorted([os.path.join(fastq_dir, fastq_file)
                          for fastq_file in os.listdir(fastq_dir)
                          if fastq_file.endswith(".fastq.gz")])

    if len(fastq_files) == 0:
        logging.warning("Could not find any fastq files")
        sys.exit()

    return fastq_files


def get_alignment_files(output_dir, fastq_files):
    alignment_files = [os.path.join(output_dir,
                                    re.sub(".fastq.gz", ".bam", os.path.basename(fastq_file)))
                       for fastq_file in fastq_files]
    return alignment_files


def check_index(genome):
    # If the index file exists.
    if genome.w_lambda and not genome.host_w_lambda_minimap2_index:
        logging.info("Index file does not exist. Creating index")
        generate_index(genome.host_w_lambda_minimap2_index, genome.host_w_lambda_genome_path)
    elif not genome.w_lambda and not genome.host_genome_path:
        logging.info("Index file does not exist. Creating index")
        generate_index(genome.host_minimap2_index, genome.host_genome_path)


def get_index(genome):
    # Make sure the genome_dir exists.
    if not os.path.isdir(genome.path):
        logging.error("No such directory %s" % genome.path)
        logging.error("Please refer to manual on generating genome structure")
        sys.exit()

    # Make sure host genome exists
    if not os.path.isdir(genome.path):
        logging.error("No such directory %s" % genome.path)
        logging.error("Please refer to manual on generating genome structure")
        sys.exit()
    if not os.path.isfile(os.path.join(genome.host_genome_path)):
        logging.error("No such file %s" % os.path.join(genome.genome_host_path))
        logging.error("Please refer to manual on generating genome structure")
        sys.exit()

    # Make sure lambda genome exists
    if genome.w_lambda and not os.path.isdir(genome.lambda_path):
        logging.info("Lambda directory does not exist but lambda specified. Creating directory")
        os.mkdir(genome.lambda_path)
    if genome.w_lambda and not os.path.isfile(genome.lambda_genome_path):
        logging.info("Cannot find lambda genome but lambda filtering was specified")
        logging.info("Downloading lambda genome")
        subprocess.run(['curl', '-O', lambda_genome], cwd=genome.lambda_path)
        logging.info("Unzipping lambda genome")
        subprocess.run(['gzip', '--decompress', os.path.join(genome.lambda_path, os.path.basename(lambda_genome))])
        shutil.move(os.path.join(genome.lambda_path, re.sub(".gz", "", os.path.basename(lambda_genome))),
                    genome.lambda_genome_path)
    # Run fasta index (used as a bed file for later)
    if genome.w_lambda and not os.path.isfile(genome.lambda_samtools_index):
        subprocess.run(["samtools", "faidx", genome.lambda_genome_path])
    # Create a genome file
    if genome.w_lambda and not os.path.isfile(genome.lambda_genome_bed):
        with open(genome.lambda_genome_path, 'r') as in_f, \
                open(genome.lambda_genome_bed, 'w') as out_f:
            for record in SeqIO.parse(in_f, 'fasta'):
                out_f.write('{}\t0\t{}\n'.format(record.id, len(record)))
    # Combine lambda and genome file
    if genome.w_lambda and not os.path.isfile(genome.host_w_lambda_genome_path):
        combine_lambda_genome(genome)

    # Create index if it doesn't exist
    check_index(genome)


def combine_lambda_genome(genome):
    # Open up lambda genome
    logging.info("Appending lambda genome to current genome")

    with open(genome.host_w_lambda_genome_path, 'w') as combined_out:
        with open(genome.host_genome_path, 'r') as genome_in:
            for fasta in SeqIO.parse(genome_in, 'fasta'):
                SeqIO.write(fasta, combined_out, "fasta")
        with open(genome.lambda_genome_path, 'r') as genome_in:
            fasta = SeqIO.read(genome_in, 'fasta')
            SeqIO.write(fasta, combined_out, "fasta")
    logging.info("Appending complete, combo file is in %s" % genome.host_w_lambda_genome_path)


def merge_bams(sample):
    # Get list of files to merge into and files to be merged
    merged_bam_file_list = [sample.host_merged_bam_file, sample.unaligned_merged_bam_file]
    alignment_file_list = [sample.host_alignment_files, sample.unaligned_files]

    if sample.w_lambda:
        merged_bam_file_list.append(sample.lambda_merged_bam_file)
        alignment_file_list.append(sample.lambda_alignment_files)

    # Run merge command for each in list of merged and files to be merged
    for merged_file, alignment_files in zip(merged_bam_file_list, alignment_file_list):
        logging.info("Merging %d bams into %s" % (len(alignment_files), merged_file))
        samtools_merge_command = ["samtools", "merge", "-f", merged_file]
        samtools_merge_command.extend(alignment_files)
        # noinspection PyArgumentList
        samtools_merge_proc = subprocess.run(samtools_merge_command, capture_output=True)
        # Merge processes
        if not samtools_merge_proc.returncode == 0:
            logging.warning("Error, merging did not go smoothly")
            logging.warning("Command: %s" % ' '.join(samtools_merge_command))
            logging.warning("Stderr: %s" % samtools_merge_proc.stderr.decode())
        # Index merged bam files
        samtools_index_command = ['samtools', 'index', merged_file]
        # noinspection PyArgumentList
        samtools_index_proc = subprocess.run(samtools_index_command, capture_output=True)
        # Merge processes
        if not samtools_index_proc.returncode == 0:
            logging.warning("Error, indexing of merged bam did not go smoothly")
            logging.warning("Command: %s" % ' '.join(samtools_index_command))
            logging.warning("Stderr: %s" % samtools_index_proc.stderr.decode())


def merge_pickles(sample):
    logging.info("Merging pickles")
    # Get list of files to merge into and files to be merged
    merged_pickle_file_list = [sample.host_merged_pickle_file]
    alignment_pickle_list = [sample.host_alignment_pickles]

    if sample.w_lambda:
        merged_pickle_file_list.append(sample.lambda_merged_pickle_file)
        alignment_pickle_list.append(sample.lambda_alignment_pickles)

    for merged_file, alignment_files in zip(merged_pickle_file_list, alignment_pickle_list):
        # Import pickles
        loaded_pickles = [wub.util.misc.pickle_load(pickle) for pickle in alignment_files]
        # Merge pickles
        merged_pickle = wub.util.misc.merge_pickles(loaded_pickles)
        # Dump merged pickle
        wub.util.misc.pickle_dump(obj=merged_pickle, fname=merged_file)


def run_wubber_multiqc(sample, qc_dir):
    # Run multiqc using all the generated pickle files in the directory.
    """
    usage: bam_multi_qc [-h] [-r report_pdf] [-x]
                    [input_pickles [input_pickles ...]]
    """
    bam_multi_qc_command = ["bam_multi_qc.py",
                            '-r', os.path.join(qc_dir, sample.sample_prefix + ".multiqc.pdf")]
    bam_multi_qc_command.extend(sample.merged_pickles)
    # Write bam
    logging.info("Performing multiqc on bam alignments %s" % ' '.join(bam_multi_qc_command))
    # noinspection PyArgumentList
    bam_multi_qc_proc = subprocess.run(bam_multi_qc_command, capture_output=True)
    if not bam_multi_qc_proc.returncode == 0:
        logging.warning("Bam multiqc returned non-zero exit code")
        logging.warning("Stderr: %s" % bam_multi_qc_proc.stderr.decode())


def main(args):
    # Make sure files and directories  exist
    check_args(args)

    # Get / create index
    genome = Genome(os.path.normpath(args.genome_dir), os.path.normpath(args.genome), args.w_lambda)
    get_index(genome)

    # Get fastq and alignment files
    sample = Sample(os.path.normpath(args.fastq_dir), os.path.normpath(args.output_dir), genome)

    # Reduce thread count unless already 1.
    threads = 1 if args.threads == 1 else args.threads - 1
    logging.info("Given we need to take of the parent script, "
                 "running %d jobs in parallel" % threads)
    # Create host dir
    if not os.path.isdir(os.path.join(args.output_dir, args.genome)):
        os.mkdir(os.path.join(args.output_dir, args.genome))

    # Create lambda_dir
    if args.w_lambda and not os.path.isdir(os.path.join(args.output_dir, "lambda")):
        os.mkdir(os.path.join(args.output_dir, "lambda"))

    # Run commands in parallel
    # Run in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        iterator = {executor.submit(generate_alignment, subfolder, args.md, args.cs):
                    subfolder for subfolder in sample.alignment_objects}
        for item in concurrent.futures.as_completed(iterator):
            parameter_input = iterator[item]
            success = item.result()

    # Merge bams
    merge_bams(sample)

    # Create QC dir
    qc_dir = os.path.join(args.output_dir, "wub")
    if not os.path.isdir(qc_dir):
        os.mkdir(qc_dir)

    # Run wubqc in parallel
    # Run in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        iterator = {executor.submit(generate_pickle, subfolder):
                        subfolder for subfolder in sample.alignment_objects}
        for item in concurrent.futures.as_completed(iterator):
            parameter_input = iterator[item]
            success = item.result()

    # Merge pickles
    merge_pickles(sample)

    # Run multiqc
    run_wubber_multiqc(sample, qc_dir)


if __name__ == "__main__":
    main()