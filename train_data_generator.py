import argparse
import math
import time
import os
import sys
import multiprocessing

from modules.CandidateFinder import CandidateFinder
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler
from modules.AlleleFinder import AlleleFinder
from modules.VcfHandler import VCFFileProcessor
from modules.CandidateLabeler import CandidateLabeler
from modules.BedHandler import BedHandler
from modules.TextColor import TextColor

"""
alignmentPolish finds possible variant sites in given bam file.

It requires three parameters:
- bam_file_path: path to a bam file
- reference_file_path: path to a reference file

Creates:
- CandidateFinder object that contains windows of possible variants.


Also, the terms "window" and "region" are NOT interchangeable.
Region: A genomic region of interest where we want to find possible variant candidate
Window: A window in genomic region where there can be multiple alleles

A region can have multiple windows and each window belongs to a region.

 Example Usage:
 python3 main.py --bam [path_to_bam] --ref [path_to_reference_fasta_file] --chromosome_name chr3 --max_threads [max_number_of_threads] --test [True/False] --output_dir [path_to_output] 
"""

DEBUG_PRINT_WINDOWS = False
DEBUG_PRINT_CANDIDATES = False
DEBUG_TIME_PROFILE = False


class AllCandidatesInRegion:
    """
    Creates a list of candidates in a region.
    """
    def __init__(self, chromosome_name, start_position, end_position):
        """
        Initialize object
        :param chromosome_name: Name of the chromosome
        :param start_position: Region start
        :param end_position: Region end
        """
        self.chromosome_name = chromosome_name
        self.start_position = start_position
        self.end_position = end_position
        self.all_candidates = []

    def add_candidate_to_list(self, alignment_candidates_tuple):
        """
        Add a candidate to the list
        :param alignment_candidates_object: Candidate object to add
        :return:
        """
        self.all_candidates.append(alignment_candidates_tuple)


class View:
    """
    Works as a main class and handles user interaction with different modules.
    """
    def __init__(self, chromosome_name, bam_file_path, reference_file_path, output_file_path, vcf_file_path):
        # --- initialize handlers ---
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)
        self.output_dir = output_file_path
        self.vcf_handler = VCFFileProcessor(file_path=vcf_file_path)

        # --- initialize parameters ---
        self.chromosome_name = chromosome_name

    def get_labeled_candidate_sites(self, AllCandidatesInRegion_object, filter_hom_ref=False):
        """
        Takes a dictionary of allele data and compares with a VCF to determine which candidate alleles are supported.
        :param candidate_dictionary: dictionary with the list of allele sites under the top-level entry "all_candidates"
        :param filter_hom_ref: whether to ignore hom_ref VCF records during candidate validation
        :return: labeled_sites: the parsed candidate list with the following structure for each entry:

        [chromosome_name, start, stop, is_insert, ref_seq, alt1, alt2, gt1, gt2]
        """
        candidate_sites = AllCandidatesInRegion_object.all_candidates

        # find start and stop position of region that covers all candidates

        chromosome_name = AllCandidatesInRegion_object.chromosome_name
        start = AllCandidatesInRegion_object.start_position
        stop = AllCandidatesInRegion_object.end_position

        # get dictionary of variant records for full region
        self.vcf_handler.populate_dictionary(contig=chromosome_name,
                                             start_pos=start,
                                             end_pos=stop,
                                             hom_filter=filter_hom_ref)

        # get separate positional variant dictionaries for IN, DEL, and SNP
        positional_variants = self.vcf_handler.get_variant_dictionary()

        allele_selector = CandidateLabeler(fasta_handler=self.fasta_handler)

        labeled_sites = allele_selector.get_labeled_candidates(chromosome_name=chromosome_name,
                                                               positional_vcf=positional_variants,
                                                               candidate_sites=candidate_sites)

        return labeled_sites

    def parse_region(self, start_position, end_position):
        """
        Find possible candidate windows.
        - All candidate lists
        """
        reads = self.bam_handler.get_reads(chromosome_name=self.chromosome_name,
                                           start=start_position,
                                           stop=end_position)

        candidate_finder = CandidateFinder(reads=reads,
                                           fasta_handler=self.fasta_handler,
                                           chromosome_name=self.chromosome_name,
                                           region_start_position=start_position,
                                           region_end_position=end_position)
        # parse reads to find candidate positions
        candidate_finder.parse_reads(reads=reads)
        # merge candidate positions
        candidate_finder.merge_positions()
        # print the windows we got
        if DEBUG_PRINT_WINDOWS:
            candidate_finder.print_windows()

        candidate_windows = candidate_finder.get_candidate_windows()
        all_candidate_lists = AllCandidatesInRegion(self.chromosome_name, start_position, end_position)

        # for each window find list of possible alleles
        for chr_name, window_start, window_end in candidate_windows:
            # get the reference sequence
            reference_sequence = self.fasta_handler.get_sequence(chr_name, window_start, window_end+1)
            # get all pileup columns in that window
            pileup_columns = self.bam_handler.get_pileupcolumns_aligned_to_a_region(chr_name, window_start, window_end+1)

            allele_finder = AlleleFinder(chr_name, window_start, window_end, pileup_columns, reference_sequence)
            # generate base dictionaries
            allele_finder.generate_base_dictionaries()
            # generate candidate allele list
            in_alleles, snp_alleles, del_alleles = allele_finder.generate_candidate_allele_list()

            if DEBUG_PRINT_CANDIDATES:
                print(chr_name, window_start, window_end, "INs: ", in_alleles, "SNPs: ", snp_alleles, 'DELs', del_alleles)

            # add alleles to candidate
            all_candidate_lists.add_candidate_to_list((window_start, window_end, in_alleles, snp_alleles, del_alleles))

        labeled_sites = self.get_labeled_candidate_sites(all_candidate_lists, True)

        return labeled_sites

    def test(self):
        """
        Run a test
        :return:
        """
        self.parse_region(start_position=100000, end_position=200000)


def write_bed(chromosome_name, bedTools_object, output_dir):
    """
    Create a bed output of all candidates found in the region
    :param start: Candidate region start
    :param end: Candidate region end
    :param all_candidate_lists: Candidate list to be saved
    :return:
    """
    if not os.path.exists(output_dir + "bed_output/"):
        os.mkdir(output_dir + "bed_output/")
    bedTools_object.saveas(output_dir + "bed_output/" + "Labeled_sites" + '_' + chromosome_name + ".bed")


def parallel_run(args):
    """
    This method is run in parallel
    :param args: Tuple containing arguments for parallelization
    :return: Result got from that region
    """
    chr_name, bam_file, ref_file, output_dir, vcf_file, start_pos, end_pos = args

    # create a view object
    view_ob = View(chromosome_name=chr_name,
                   bam_file_path=bam_file,
                   reference_file_path=ref_file,
                   output_file_path=output_dir,
                   vcf_file_path=vcf_file)

    # return the results
    return view_ob.parse_region(start_pos, end_pos)


def chromosome_level_parallelization(chr_name, bam_file, ref_file, vcf_file, output_dir, max_threads):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param chr_name: Chromosome name
    :param bam_file: Bam file
    :param ref_file: Ref file
    :param vcf_file: VCF file
    :param json_out: If output as json or not
    :param output_dir: Output directory
    :param max_threads: Maximum number of threads
    :return: A list of results returned by the processes
    """
    # entire length of chromosome
    fasta_handler = FastaHandler(ref_file)
    whole_length = fasta_handler.get_chr_sequence_length(chr_name)
    whole_length = 8000000 # for testing purpose

    # expected length of each segment
    each_segment_length = int(math.ceil(whole_length / max_threads))
    args = list()

    for i in range(max_threads):
        # parse window of the segment. Use a 1000 overlap for corner cases.
        start_position = i * each_segment_length
        end_position = min((i + 1) * each_segment_length + 1000, whole_length)
        args.append((chr_name, bam_file, ref_file, output_dir, vcf_file, start_position, end_position))

    # create a pool of workers
    pool = multiprocessing.Pool(processes=max_threads)

    # run and get results of those threads
    results = pool.map(parallel_run, args)

    # wait for all the processes to finish
    pool.close()
    pool.join()

    # return results
    return results


def genome_level_parallelization(bam_file, ref_file, vcf_file, output_dir, max_threads):
    """
    This method calls chromosome_level_parallelization for each chromosome.
    :param chr_name: Character name
    :param bam_file: BAM file path
    :param ref_file: Reference file path
    :param vcf_file: VCF file path
    :param json_out: JSON out parameter
    :param output_dir: Output directory
    :param max_threads: Maximum number of threads to create in chromosome level
    :return: Saves a bed file
    """
    # chr_list = ["chr1", "chr2", "chr3", "chr4", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    #             "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    program_start_time = time.time()

    chr_list = ["chr3"]

    labeled_sites = []
    # each chormosome in list
    for chr in chr_list:
        sys.stderr.write(TextColor.BLUE + "STARTING " + str(chr) + " PROCESSES" + "\n")
        start_time = time.time()

        # do a chromosome level parallelization
        results = chromosome_level_parallelization(chr, bam_file, ref_file, vcf_file, output_dir, max_threads)

        # add all results to labeled sites
        for result in results:
            labeled_sites.extend(result)

        end_time = time.time()
        sys.stderr.write(TextColor.PURPLE + "FINISHED " + str(chr) + " PROCESSES" + "\n")
        sys.stderr.write(TextColor.CYAN + "TIME ELAPSED: " + str(end_time-start_time) + "\n")

    sys.stderr.write(TextColor.BLUE + "CONVERTING TO BED FILE" + "\n")

    bed_file = BedHandler.list_to_bed(labeled_sites)
    write_bed("Whole_genome", bed_file, output_dir)

    program_end_time = time.time()
    sys.stderr.write(TextColor.RED + "PROCESSED FINISHED SUCESSFULLY" + "\n")
    sys.stderr.write(TextColor.CYAN + "TOTAL TIME FOR GENERATING ALL RESULTS: " + str(program_end_time-program_start_time) + "\n")



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

    FLAGS, unparsed = parser.parse_known_args()
    # process the output directory
    if FLAGS.output_dir[-1] != '/':
        FLAGS.output_dir += '/'
    if not os.path.exists(FLAGS.output_dir):
        os.mkdir(FLAGS.output_dir)

    view = View(chromosome_name=FLAGS.chromosome_name,
                bam_file_path=FLAGS.bam,
                reference_file_path=FLAGS.ref,
                output_file_path=FLAGS.output_dir,
                vcf_file_path=FLAGS.vcf)

    if FLAGS.test is True:
        view = View(chromosome_name=FLAGS.chromosome_name,
                    bam_file_path=FLAGS.bam,
                    reference_file_path=FLAGS.ref,
                    output_file_path=FLAGS.output_dir,
                    vcf_file_path=FLAGS.vcf)
        view.test()
    elif FLAGS.chromosome_name is not None:
        chromosome_level_parallelization(FLAGS.chromosome_name, FLAGS.bam, FLAGS.ref,
                                         FLAGS.vcf, FLAGS.output_dir, FLAGS.max_threads)
    else:
        genome_level_parallelization(FLAGS.bam, FLAGS.ref, FLAGS.vcf, FLAGS.output_dir, FLAGS.max_threads)
