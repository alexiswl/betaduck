#!/usr/bin/env python3

import os
import subprocess
import logging
import sys
import re
import concurrent.futures
import shutil
from Bio import SeqIO

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
    minimap2_index_proc = subprocess.run(minimap2_index_command, capture_output=True)

    # Check output of command
    if not minimap2_index_proc.returncode == 0:
        logger.warning("Index command returned non-zero exit code")
        logger.warning("Stderr: %s" % minimap2_index_proc.stderr.decode())


def generate_alignment(genome_dir, genome, index, input_fastq, alignment_file, md, cs, w_lambda=False):
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
    minimap2_alignment_command.extend([index, input_fastq])
    logging.info("Running command %s minimap2_alignment_command" % minimap2_alignment_command)
    logging.info("Writing to %s" % alignment_file)
    logging.info("Completed command %s" % minimap2_alignment_command)
    # Write to alignment_file
    # Write to sorted bam file and sort bam flie
    dev_null = open(os.devnull, 'w')
    if not os.path.isdir(os.path.join(os.path.dirname(alignment_file), "unaligned")):
        os.mkdir(os.path.join(os.path.dirname(alignment_file), "unaligned"))
    alignment_file_no_alignment = os.path.join(os.path.dirname(alignment_file), "unaligned",
                                               re.sub(".bam$", ".unaligned.bam", os.path.basename(alignment_file)))
    p1 = subprocess.Popen(minimap2_alignment_command,
                          stdout=subprocess.PIPE, stderr=dev_null)
    p2 = subprocess.Popen(["samtools", "view", '-bu', '-F4', '-U', alignment_file_no_alignment],
                          stdin=p1.stdout,
                          stdout=subprocess.PIPE, stderr=dev_null)
    p1.stdout.close()
    if not w_lambda:
        alignment_file = os.path.join(os.path.dirname(alignment_file), genome,
                                      re.sub(".bam", ".sorted.bam", os.path.basename(alignment_file)))
        p3 = subprocess.Popen(["samtools", "sort", "-o", alignment_file],
                              stdin=p2.stdout,
                              stdout=dev_null, stderr=dev_null)
        p2.stdout.close()
        out, err = p3.communicate()
        if not p3.returncode == 0:
            logging.info("BAM Alignment returned non-zero exit code for alignment file %s" % alignment_file)
        # Generate index for alignment file
        index_command = ["samtools", "index", alignment_file]
        index_proc = subprocess.run(index_command, capture_output=True)
        if not index_proc.returncode == 0:
            logging.warning("Indexing of '%s' returned non-zero exit code" % ' '.join(index_command))
            logging.warning("Stderr: %s" % index_proc.stderr.decode())
    else:
        # If w_lambda, pipe sort to stdout
        lambda_alignment_file = os.path.join(os.path.dirname(alignment_file), "lambda",
                                             re.sub(".bam$", ".lambda.sorted.bam", os.path.basename(alignment_file)))
        alignment_file_no_lambda = os.path.join(os.path.dirname(alignment_file), genome,
                                                re.sub(".bam$", ".lambda-filt.sorted.bam", os.path.basename(alignment_file)))
        lambda_bed = os.path.join(genome_dir, "lambda", "genome.bed")
        p3 = subprocess.Popen(["samtools", "sort"],
                              stdin=p2.stdout,
                              stdout=subprocess.PIPE, stderr=dev_null)
        p2.stdout.close()
        p4 = subprocess.Popen(["samtools", "view", "-b", "-L", lambda_bed,
                               "-U", alignment_file_no_lambda, "-o", lambda_alignment_file],
                              stdin=p3.stdout,
                              stdout=dev_null, stderr=dev_null)
        p3.stdout.close()
        out, err = p4.communicate()
        if not p4.returncode == 0:
            logging.warning("BAM Alignment returned non-zero exit code for alignment file %s" % alignment_file)
        # Generate index for alignment file(s)
        logging.info("Generating indexes for %s and %s" % (lambda_alignment_file, alignment_file_no_lambda))
        for alignment_file in [lambda_alignment_file, alignment_file_no_lambda]:
            index_command = ["samtools", "index", alignment_file]
            index_proc = subprocess.run(index_command, capture_output=True)
            if not index_proc.returncode == 0:
                logging.warning("Indexing of '%s' returned non-zero exit code" % ' '.join(index_command))
                logging.warning("Stderr: %s" % index_proc.stderr.decode())


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


def check_index(index_file, genome_file):
    # If the index file exists.
    if not os.path.isfile(index_file):
        logging.info("Index file does not exist. Creating index")
        generate_index(index_file, genome_file)


def get_index(genome_dir, genome, w_lambda):
    # Index is just the reference genome path with a different suffix
    if not w_lambda:
        index_file = os.path.join(genome_dir, genome, "genome.mmi")
        genome_file = os.path.join(genome_dir, genome, "genome.fa")
    else:
        index_file = os.path.join(genome_dir, genome, "genome.w_lambda.mmi")
        genome_file = os.path.join(genome_dir, genome, "genome.w_lambda.fa")

    # Make sure the genome_dir exists.
    if not os.path.isdir(genome_dir):
        logging.error("No such directory %s" % genome_dir)
        logging.error("Please refer to manual on generating genome structure")
        sys.exit()

    # Make sure host genome exists
    if not os.path.isdir(os.path.join(genome_dir, genome)):
        logging.error("No such directory %s" % os.path.join(genome_dir, genome))
        logging.error("Please refer to manual on generating genome structure")
        sys.exit()
    if not os.path.isfile(os.path.join(genome_dir, genome, "genome.fa")):
        logging.error("No such file %s" % os.path.join(genome_dir, genome, "genome.fa"))
        logging.error("Please refer to manual on generating genome structure")
        sys.exit()

    # Make sure lambda genome exists
    if w_lambda and not os.path.isdir(os.path.join(genome_dir, "lambda")):
        logging.info("Lambda directory does not exist but lambda specified. Creating directory")
        os.mkdir(os.path.join(genome_dir, "lambda"))
    if w_lambda and not os.path.isfile(os.path.join(genome_dir, "lambda", "genome.fa")):
        logging.info("Cannot find lambda genome but lambda filtering was specified")
        logging.info("Downloading lambda genome")
        subprocess.run(['curl', '-O', lambda_genome], cwd=os.path.join(genome_dir, "lambda"))
        logging.info("Unzipping lambda genome")
        subprocess.run(['gzip', '--decompress', os.path.join(genome_dir, "lambda", os.path.basename(lambda_genome))])
        shutil.move(os.path.join(genome_dir, "lambda", re.sub(".gz", "", os.path.basename(lambda_genome))),
                    os.path.join(genome_dir, "lambda", "genome.fa"))
    # Run fasta index (used as a bed file for later)
    if w_lambda and not os.path.isfile(os.path.join(genome_dir, "lambda", "genome.fa.fai")):
        subprocess.run(["samtools", "faidx", os.path.join(genome_dir, "lambda", "genome.fa")])
    # Create a genome file
    if w_lambda and not os.path.isfile(os.path.join(genome_dir, "lambda", "genome.bed")):
        with open(os.path.join(genome_dir, "lambda", "genome.fa")) as in_f, \
                open(os.path.join(genome_dir, "lambda", "genome.bed"),'w') as out_f:
            for record in SeqIO.parse(in_f, 'fasta'):
                out_f.write('{}\t0\t{}\n'.format(record.id, len(record)))
    # Combine lambda and genome file
    if w_lambda and not os.path.isfile(genome_file):
        combine_lambda_genome(genome_dir, genome)

    # Create index if it doesn't exist
    check_index(index_file, genome_file)

    return index_file


def combine_lambda_genome(genome_dir, genome):
    # Open up lambda genome
    logging.info("Appending lambda genome to current genome")
    genome_file = os.path.join(genome_dir, genome, "genome.fa")
    combined_file = os.path.join(genome_dir, genome, "genome.w_lambda.fa")
    lambda_genome_file = os.path.join(genome_dir, "lambda", "genome.fa")
    with open(combined_file, 'w') as combined_out:
        with open(genome_file, 'r') as genome_in:
            for fasta in SeqIO.parse(genome_in, 'fasta'):
                SeqIO.write(fasta, combined_out, "fasta")
        with open(lambda_genome_file, 'r') as genome_in:
            fasta = SeqIO.read(genome_in, 'fasta')
            SeqIO.write(fasta, combined_out, "fasta")
    logging.info("Appending complete, combo file is in %s" % combined_file)


def merge_bams(input_dir, genome, output_dir, output_prefix, w_lambda=False):
    if w_lambda:
        unaligned_output_prefix = re.sub('.lambda-filt.sorted.merged.bam$', '.sorted.merged.bam', output_prefix)
        unaligned_output_prefix = re.sub('^%s_' % genome, "unaligned_", unaligned_output_prefix)
    else:
        unaligned_output_prefix = re.sub('.sorted.merged.bam$', '.merged.bam', output_prefix)
        unaligned_output_prefix = re.sub('^%s_' % genome, "unaligned_", unaligned_output_prefix)

    sorted_bam_files = sorted([os.path.join(input_dir, genome, bam)
                               for bam in os.listdir(os.path.join(input_dir, genome))
                               if bam.endswith(".sorted.bam")])
    samtools_merge_command = ["samtools", "merge", "-f", os.path.join(output_dir, output_prefix)]
    samtools_merge_command.extend(sorted_bam_files)

    # Run merge command
    logging.info("Merging %s bams: %s" % (genome, ' '.join(samtools_merge_command)))
    samtools_merge_proc = subprocess.run(samtools_merge_command, capture_output=True)
    if not samtools_merge_proc.returncode == 0:
        logging.warning("Error, merging did not go smoothly")
        logging.warning("Stderr: %s" % samtools_merge_proc.stderr.decode())

    # Merge lambdas
    if w_lambda:
        lambda_output_prefix = re.sub('.lambda-filt.sorted.merged.bam$', '.sorted.merged.bam', output_prefix)
        lambda_output_prefix = re.sub('^%s_' % genome, "lambda_", lambda_output_prefix)
        sorted_bam_files = sorted([os.path.join(input_dir, "lambda", bam)
                                   for bam in os.listdir(os.path.join(input_dir, "lambda"))
                                   if bam.endswith(".sorted.bam")])
        samtools_merge_command = ["samtools", "merge", "-f", os.path.join(output_dir, lambda_output_prefix)]
        samtools_merge_command.extend(sorted_bam_files)
        # Run merge command
        logging.info("Merging lambda bams: %s" % ' '.join(samtools_merge_command))
        samtools_merge_proc = subprocess.run(samtools_merge_command, capture_output=True)
        if not samtools_merge_proc.returncode == 0:
            logging.warning("Error, merging did not go smoothly")
            logging.warning("Stderr: %s" % samtools_merge_proc.stderr.decode())

    # Merge the unaligned bam files
    sorted_bam_files = sorted([os.path.join(input_dir, "unaligned", bam)
                               for bam in os.listdir(os.path.join(input_dir, genome))
                               if bam.endswith(".bam")])
    samtools_merge_command = ["samtools", "merge", "-f", os.path.join(output_dir, unaligned_output_prefix)]
    samtools_merge_command.extend(sorted_bam_files)

    # Run merge command
    logging.info("Merging unaligned bams: %s" % ' '.join(samtools_merge_command))
    samtools_merge_proc = subprocess.run(samtools_merge_command, capture_output=True)
    if not samtools_merge_proc.returncode == 0:
        logging.warning("Error, merging did not go smoothly")
        logging.warning("Stderr: %s" % samtools_merge_proc.stderr.decode())


def run_wubber(input_dir, genome_dir, genome, output_dir, qc_prefix, w_lambda=False):
    # Get reference
    reference = os.path.join(genome_dir, genome, "genome.fa")

    # Run the wubber
    """
    bam_alignment_qc [-h] -f reference [-c region] [-n context_sizes] [-x]
                        [-t bam_tag] [-q aqual] [-i qual_ints] [-r report_pdf]
                        [-p results_pickle] [-Q]
                        bam
    """
    bam_files = [bam_file
                 for bam_file in os.listdir(input_dir)
                 if bam_file.endswith(".merged.bam")]
    for bam_file in bam_files:
        if bam_file.startswith("lambda") and w_lambda:
            bam_reference = os.path.join(genome_dir, "lambda", "genome.fa")
        else:
            bam_reference = reference
        bam_prefix = '_'.join([bam_file.rsplit("_", 3)[0], qc_prefix])
        bam_alignment_qc_command = ["bam_alignment_qc",
                                    "-f", bam_reference,
                                    "-Q", '-p', os.path.join(output_dir, '.'.join([bam_prefix, "pickle"])),
                                    "-r", os.path.join(output_dir, '.'.join([bam_prefix, "report", "pdf"])),
                                    bam_file]
        logging.info("Performing bam alignment qc: %s" % ' '.join(bam_alignment_qc_command))
        bam_alignment_qc_proc = subprocess.run(bam_alignment_qc_command, capture_output=True)
        if not bam_alignment_qc_proc.returncode == 0:
            logging.warning("Bam QC returned non-zero exit code")
            logging.warning("Stderr: %s" % bam_alignment_qc_proc.stderr.decode())

    # Run multiqc using all the pickle files
    pickles = [os.path.join(output_dir, pickle)
               for pickle in os.listdir(output_dir)
               if pickle.endswith(".pickle")]
    """
    usage: bam_multi_qc [-h] [-r report_pdf] [-x]
                    [input_pickles [input_pickles ...]]
    """
    bam_multi_qc_command = ["bam_multi_qc",
                            '-r', os.path.join(output_dir, qc_prefix + ".multiqc.pdf")]
    bam_multi_qc_command.extend(pickles)
    # Write bam
    logging.info("Performing multiqc on bam alignments %s" % ' '.join(bam_multi_qc_command))
    bam_multi_qc_proc = subprocess.run(bam_multi_qc_command, capture_output=True)
    if not bam_multi_qc_proc.returncode == 0:
        logging.warning("Bam multiqc returned non-zero exit code")
        logging.warning("Stderr: %s" % bam_multi_qc_proc.stderr.deocde())


def main(args):
    # Make sure files and directories  exist
    check_args(args)
    # Get fastq files
    fastq_files = collect_fastqs(args.fastq_dir)
    # Get alignmnet files
    alignment_files = get_alignment_files(args.output_dir, fastq_files)
    # Get / create index
    index = get_index(args.genome_dir, args.genome, args.w_lambda)
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
        iterator = {executor.submit(generate_alignment, args.genome_dir, args.genome, index, fastq_file, alignment_file, args.md, args.cs, args.w_lambda):
                    (fastq_file, alignment_file) for fastq_file, alignment_file in zip(fastq_files, alignment_files)}
        for item in concurrent.futures.as_completed(iterator):
            parameter_input = iterator[item]
            success = item.result()

    # Merge the bam files
    flowcell, rnumber = alignment_files[0].split("_", 2)[:2]
    if args.w_lambda:
        alignment_file_prefix = '_'.join([args.genome, flowcell, rnumber]) + ".lambda-filt.sorted.merged.bam"
    else:
        alignment_file_prefix = '_'.join([args.genome, flowcell, rnumber]) + ".sorted.merged.bam"
    if not os.path.isdir(os.path.join(args.output_dir, "merged")):
        os.mkdir(os.path.join(args.output_dir, "merged"))
    # Merge bams
    merge_bams(args.output_dir, args.genome, os.path.join(args.output_dir, "merged"),
               alignment_file_prefix, w_lambda=args.w_lambda)
    # Create QC dir
    if not os.path.isdir(os.path.join(args.output_dir, "wub")):
        os.mkdir(os.path.join(args.output_dir, "wub"))
    # Run wubqc
    qc_prefix = '_'.join([args.genome, flowcell, rnumber])
    run_wubber(args.output_dir, args.genome_dir, args.genome,
               os.path.join(args.output_dir, "wub"), qc_prefix,
               w_lambda=args.w_lambda)


if __name__ == "__main__":
    main()