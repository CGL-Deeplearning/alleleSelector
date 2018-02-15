import argparse
import math
import time
import os
import sys
import multiprocessing
import numpy
import copy

from modules.FilterCandidateFinder import CandidateFinder
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler
from modules.VcfHandler import VCFFileProcessor
from modules.FilterCandidateLabeler import CandidateLabeler
from modules.TextColor import TextColor
from modules.FileManager import FileManager
from modules.BedHandler import BedHandler
from modules.TsvHandler import TsvHandler

"""
alleleSelector finds possible variant sites in given bam file.
This script is responsible of creating candidate sites with true genotypes for neural network training.

It requires three parameters:
- bam_file_path: path to a bam file
- reference_file_path: path to a reference file
- vcf_file_path: path to a VCF file for true genotype labeling

Creates:
- Bed files containing candidate sites and their true genotypes for training.


Also, the terms "window" and "region" are NOT interchangeable.
Region: A genomic region of interest where we want to find possible variant candidate
Window: A window in genomic region where there can be multiple alleles

A region can have multiple windows and each window belongs to a region.

Example usage:

local:
python3 train_data_generator_filter.py --ref /Users/saureous/data/chr1.fa --bam /Users/saureous/data/chr1.sorted.bam --vcf /Users/saureous/data/NA12878_S1.genome.vcf.gz --bed /Users/saureous/data/ConfidentRegions.hg19.chr1.bed --chromosome_name chr1 --edit_threshold 2 --max_threads 1 >out.txt

jarvis:
python3 train_data_generator_filter.py --ref /data/users/ryan/chr1.fa --bam /data/users/ryan/chr1.sorted.bam --vcf /data/users/ryan/NA12878_S1.genome.vcf.gz --bed /data/users/ryan/ConfidentRegions.hg19.chr1.bed --chromosome_name chr1 --edit_threshold 4 --max_threads 80

"""

DEBUG_PRINT_CANDIDATES = False
DEBUG_TIME_PROFILE = False
DEBUG_SAVE_VCF = False


class View:
    """
    Works as a main class and handles user interaction with different modules.
    """
    def __init__(self, chromosome_name, bam_file_path, reference_file_path, output_file_path, vcf_file_path, bed_path, MIN_MISMATCH_PERCENT_THRESHOLD):
        # --- initialize handlers ---
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)
        self.output_dir = output_file_path
        self.vcf_handler = VCFFileProcessor(file_path=vcf_file_path)
        self.tsv_handler = TsvHandler(bed_path)

        # --- initialize parameters ---
        self.chromosome_name = chromosome_name
        self.MIN_MISMATCH_PERCENT_THRESHOLD = MIN_MISMATCH_PERCENT_THRESHOLD

    def write_training_set(self, frequencies, start, stop):
        """
        Create a npz training set of all labeled candidate sites found in the region
        :param frequencies: vectors of the form [f1, f2, f3, ..., fn, L] where f = freq and L = label 0/1
        :param start: start coord of region
        :param stop: stop coord of region
        :return:
        """
        filename = "allele_frequencies_" + str(start) + "_" + str(stop) + ".npz"
        numpy.savez_compressed(os.path.join(self.output_dir, filename), a=frequencies)

    def get_labeled_candidate_sites(self, selected_candidate_list, start_pos, end_pos, filter_hom_ref=True):
        """
        Takes a dictionary of allele data and compares with a VCF to determine which candidate alleles are supported.
        :param selected_candidate_list: List of all selected candidates with their alleles
        :param filter_hom_ref: whether to ignore hom_ref VCF records during candidate validation
        :param start_pos: start position of the region
        :param end_pos: end position of the region
        :return: labeled_sites: the parsed candidate list with the following structure for each entry:

        [chromosome_name, start, stop, is_insert, ref_seq, alt1, alt2, gt1, gt2]
        """
        # get dictionary of variant records for full region
        self.vcf_handler.populate_dictionary(contig=self.chromosome_name,
                                             start_pos=start_pos,
                                             end_pos=end_pos,
                                             hom_filter=filter_hom_ref)

        # get separate positional variant dictionaries for IN, DEL, and SNP
        positional_variants = self.vcf_handler.get_variant_dictionary()

        if DEBUG_SAVE_VCF:
            vcf_filename = "vcf_bed_" + str(start_pos) + "_" + str(end_pos)
            self.vcf_handler.save_positional_vcf_as_bed(self.chromosome_name, os.path.join(self.output_dir,vcf_filename))

        intervals = self.tsv_handler.get_subset_of_bed_intervals(start=start_pos,
                                                                 stop=end_pos,
                                                                 universal_offset=-1,
                                                                 start_offset=1)

        allele_labeler = CandidateLabeler(fasta_handler=self.fasta_handler)

        labeled_sites = allele_labeler.get_labeled_candidates(chromosome_name=self.chromosome_name,
                                                              positional_vcf=positional_variants,
                                                              candidate_sites=selected_candidate_list,
                                                              confident_intervals=intervals)
        print(labeled_sites.shape)

        return labeled_sites

    def parse_region(self, start_position, end_position):
        """
        Iterate through all the reads that fall in a region, find candidates, label candidates and output a bed file.
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :return:
        """
        reads = self.bam_handler.get_reads(chromosome_name=self.chromosome_name,
                                           start=start_position,
                                           stop=end_position)

        candidate_finder = CandidateFinder(reads=reads,
                                           fasta_handler=self.fasta_handler,
                                           chromosome_name=self.chromosome_name,
                                           region_start_position=start_position,
                                           region_end_position=end_position,
                                           MIN_MISMATCH_PERCENT_THRESHOLD=self.MIN_MISMATCH_PERCENT_THRESHOLD)

        # go through each read and find candidate positions and alleles
        selected_candidates = candidate_finder.parse_reads_and_select_candidates(reads=reads)

        if DEBUG_PRINT_CANDIDATES:
            for candidate in selected_candidates:
                print(candidate)

        labeled_sites = self.get_labeled_candidate_sites(selected_candidates, start_position, end_position, True)

        # print(labeled_sites)
        # bed_file = BedHandler.list_to_bed(labeled_sites)
        self.write_training_set(labeled_sites, start_position, end_position)

        # a = numpy.load(os.path.join(self.output_dir, "allele_frequencies.npz"))

        # print(numpy.array_equal(a['a'],labeled_sites))
        # print(a)

    def test(self):
        """
        Run a test
        :return:
        """
        start_time = time.time()
        # self.parse_region(start_position=121400000, end_position=121600000)

        self.parse_region(start_position=100000, end_position=200000)
        end_time = time.time()
        print("TOTAL TIME ELAPSED: ", end_time-start_time)


def parallel_run(chr_name, bam_path, ref_path, output_dir, vcf_path, bed_path, start_position, end_position, MIN_MISMATCH_PERCENT_THRESHOLD):
    """
    Run this method in parallel
    :param chr_name: Chromosome name
    :param bam_path: Bam file path
    :param ref_path: Ref file path
    :param output_dir: Output directory
    :param vcf_path: VCF file path
    :param start_position: Start position
    :param end_position: End position
    :return:
    """

    # create a view object
    view_ob = View(chromosome_name=chr_name,
                   bam_file_path=bam_path,
                   reference_file_path=ref_path,
                   output_file_path=output_dir,
                   vcf_file_path=vcf_path,
                   bed_path=bed_path,
                   MIN_MISMATCH_PERCENT_THRESHOLD=MIN_MISMATCH_PERCENT_THRESHOLD)

    # return the results
    view_ob.parse_region(start_position, end_position)


def chromosome_level_parallelization(chr_name, bam_path, ref_path, vcf_path, bed_path, output_dir, max_threads, MIN_MISMATCH_PERCENT_THRESHOLD):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param chr_name: Chromosome name
    :param bam_path: Bam file path name
    :param ref_path: Ref file path name
    :param vcf_path: VCF file path name
    :param bed_path: confident bed file path name
    :param output_dir: Output directory
    :param max_threads: Maximum number of threads
    :return: A list of results returned by the processes
    """
    # entire length of chromosome
    fasta_handler = FastaHandler(ref_path)
    whole_length = fasta_handler.get_chr_sequence_length(chr_name)

    # 0.2MB segments at once
    each_segment_length = 200000

    # chunk the chromosome into 1000 pieces
    chunks = int(math.ceil(whole_length / each_segment_length))

    for i in range(chunks):
        # parse window of the segment. Use a 1000 overlap for corner cases.
        start_position = i * each_segment_length
        end_position = min((i + 1) * each_segment_length, whole_length)

        args = (chr_name, bam_path, ref_path, output_dir, vcf_path, bed_path, start_position, end_position, MIN_MISMATCH_PERCENT_THRESHOLD)

        p = multiprocessing.Process(target=parallel_run, args=args)
        p.start()

        while True:
            if len(multiprocessing.active_children()) < max_threads:
                break


def create_output_dir_for_chromosome(output_dir, chr_name):
    """
    Create an internal directory inside the output directory to dump choromosomal bed files
    :param output_dir: Path to output directory
    :param chr_name: chromosome name
    :return: New directory path
    """
    path_to_dir = output_dir + chr_name + "/"
    if not os.path.exists(path_to_dir):
        os.mkdir(path_to_dir)

    return path_to_dir


# def genome_level_parallelization(bam_file, ref_file, vcf_file, output_dir, max_threads):
#     """
#     This method calls chromosome_level_parallelization for each chromosome.
#     :param bam_file: BAM file path
#     :param ref_file: Reference file path
#     :param vcf_file: VCF file path
#     :param output_dir: Output directory
#     :param max_threads: Maximum number of threads to create in chromosome level
#     :return: Saves a bed file
#     """
#     chr_list = ["chr1", "chr2", "chr3", "chr4", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
#                 "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
#     program_start_time = time.time()
#
#     # chr_list = ["chr1"]
#
#     # each chormosome in list
#     for chr in chr_list:
#         sys.stderr.write(TextColor.BLUE + "STARTING " + str(chr) + " PROCESSES" + "\n")
#         start_time = time.time()
#
#         # create dump directory inside output directory
#         output_dir = create_output_dir_for_chromosome(output_dir, chr)
#
#         # do a chromosome level parallelization
#         chromosome_level_parallelization(chr, bam_file, ref_file, vcf_file, output_dir, max_threads)
#
#         end_time = time.time()
#         sys.stderr.write(TextColor.PURPLE + "FINISHED " + str(chr) + " PROCESSES" + "\n")
#         sys.stderr.write(TextColor.CYAN + "TIME ELAPSED: " + str(end_time - start_time) + "\n")
#
#     # wait for the last process to end before file processing
#     while True:
#         if len(multiprocessing.active_children()) == 0:
#             break
#
#     for chr in chr_list:
#         # here we dumped all the bed files
#         path_to_dir = output_dir
#         concatenated_file_name = output_dir + chr + "_labeled.bed"
#         filemanager_object = FileManager()
#         # get all bed file paths from the directory
#         file_paths = filemanager_object.get_file_paths_from_directory(path_to_dir)
#         # dump all bed files into one
#         filemanager_object.concatenate_files(file_paths, concatenated_file_name)
#         # delete all temporary files
#         filemanager_object.delete_files(file_paths)
#
#     program_end_time = time.time()
#     sys.stderr.write(TextColor.RED + "PROCESSED FINISHED SUCCESSFULLY" + "\n")
#     sys.stderr.write(TextColor.CYAN + "TOTAL TIME FOR GENERATING ALL RESULTS: " + str(program_end_time-program_start_time) + "\n")


def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # create an internal directory so we don't overwrite previous runs
    timestr = time.strftime("%m%d%Y_%H%M%S")
    internal_directory = "run_" + timestr + "/"
    output_dir = output_dir + internal_directory

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return output_dir


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file containing reads of interest."
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF file path."
    )
    parser.add_argument(
        "--bed",
        type=str,
        required=True,
        help="BED file path, containing confident regions of the VCF."
    )
    parser.add_argument(
        "--chromosome_name",
        type=str,
        help="Desired chromosome number E.g.: 3"
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    parser.add_argument(
        "--test",
        type=bool,
        default=False,
        help="If true then a dry test is run."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/",
        help="Path to output directory."
    )
    parser.add_argument(
        "--edit_threshold",
        type=int,
        default=1,
        help="Percent frequency threshold for a candidate site to be considered."
    )

    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = handle_output_directory(FLAGS.output_dir)

    # view = View(chromosome_name=FLAGS.chromosome_name,
    #             bam_file_path=FLAGS.bam,
    #             reference_file_path=FLAGS.ref,
    #             output_file_path=FLAGS.output_dir,
    #             vcf_file_path=FLAGS.vcf,
    #             MIN_MISMATCH_PERCENT_THRESHOLD=FLAGS.edit_threshold)

    # if FLAGS.test is True:
    #     view = View(chromosome_name=FLAGS.chromosome_name,
    #                 bam_file_path=FLAGS.bam,
    #                 reference_file_path=FLAGS.ref,
    #                 output_file_path=FLAGS.output_dir,
    #                 vcf_file_path=FLAGS.vcf,
    #                 MIN_MISMATCH_PERCENT_THRESHOLD=FLAGS.edit_threshold)
    #     view.test()

    if FLAGS.chromosome_name is not None:
        chromosome_level_parallelization(chr_name=FLAGS.chromosome_name,
                                         bam_path=FLAGS.bam,
                                         ref_path=FLAGS.ref,
                                         vcf_path=FLAGS.vcf,
                                         bed_path=FLAGS.bed,
                                         output_dir=FLAGS.output_dir,
                                         max_threads=FLAGS.max_threads,
                                         MIN_MISMATCH_PERCENT_THRESHOLD=FLAGS.edit_threshold)
    # else:
    #     genome_level_parallelization(FLAGS.bam, FLAGS.ref, FLAGS.vcf, FLAGS.output_dir, FLAGS.max_threads)
