import argparse
import math
import time
import json
import os

from modules.CandidateFinder import CandidateFinder
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler
from modules.AlleleFinder import AlleleFinder
from multiprocessing import Process
from modules.VcfHandler import VCFFileProcessor
from modules.CandidateLabeler import CandidateLabeler
from modules.BedHandler import BedHandler

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
 python3 main.py --bam [path_to_bam] --ref [path_to_reference_fasta_file] --chromosome_name chr3 --max_threads [max_number_of_threads] --test [True/False] --json [True/False] --output_dir [path_to_JSON_output] 
"""

DEBUG_PRINT_WINDOWS = False
DEBUG_PRINT_CANDIDATES = False


class ComplexEncoder(json.JSONEncoder):
    """
    JSON encoder for class attributes
    """
    def default(self, obj):
        if hasattr(obj, 'reprJSON'):
            return obj.reprJSON()
        else:
            return json.JSONEncoder.default(self, obj)


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

    def reprJSON(self):
        """
        Report all attributes of this object as a dictionary that can be saved as a JSON
        :return: A dictionary with key value to be saved in json format
        """
        return dict(chromosome_name=self.chromosome_name, start_position=self.start_position,
                    end_position=self.end_position, all_candidates=self.all_candidates)


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

    def write_json(self, start, end, all_candidate_lists):
        """
        Create a json output of all candidates found in the region
        :param start: Candidate region start
        :param end: Candidate region end
        :param all_candidate_lists: Candidate list to be saved
        :return:
        """
        if not os.path.exists(self.output_dir + "json_output/"):
            os.mkdir(self.output_dir + "json_output/")
        json_file = open(self.output_dir + "json_output/" + "Candidates" + '_' + self.chromosome_name + '_'
                         + str(start) + '_' + str(end) + ".json", 'w')
        json_file.write(json.dumps(all_candidate_lists.reprJSON(), cls=ComplexEncoder, indent=4, sort_keys=True))

    def write_bed(self, start, end, bedTools_object):
        """
        Create a json output of all candidates found in the region
        :param start: Candidate region start
        :param end: Candidate region end
        :param all_candidate_lists: Candidate list to be saved
        :return:
        """
        if not os.path.exists(self.output_dir + "bed_output/"):
            os.mkdir(self.output_dir + "bed_output/")

        bedTools_object.saveas(self.output_dir + "bed_output/" + "Labeled_sites" + '_' +
                               self.chromosome_name + '_'+ str(start) + '_' + str(end) + ".bed")

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

    def parse_region(self, start_position, end_position, json_out):
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

        if json_out:
            self.write_json(start_position, end_position, all_candidate_lists)

        labeled_sites = self.get_labeled_candidate_sites(all_candidate_lists, True)

        bed_file = BedHandler.list_to_bed(labeled_sites)
        self.write_bed(start_position, end_position, bed_file)

    def test(self, json_out):
        """
        Run a test
        :param json_out:
        :return:
        """
        start_time = time.time()
        self.parse_region(start_position=100000, end_position=200000, json_out=json_out)
        end_time = time.time()
        print('TOTAL TIME: ', end_time-start_time)


def do_parallel(chr_name, bam_file, ref_file, vcf_file, json_out, output_dir, max_threads):
    """
    Split chromosome in different ranges for parallel processing
    :param chr_name: Chromosome name
    :param bam_file: Bam file path
    :param ref_file: Reference file path
    :param output_dir: Directory for saving output
    :param json_out: JSON out flag
    :param max_threads: Maximum number of threads
    :return:
    """
    # entire length of chromosome
    fasta_handler = FastaHandler(ref_file)
    whole_length = fasta_handler.get_chr_sequence_length(chr_name)

    # expected length of each segment
    each_segment_length = int(math.ceil(whole_length / max_threads))

    for i in range(max_threads):
        # parse window of the segment. Use a 1000 overlap for corner cases.
        view = View(chromosome_name=chr_name,
                    bam_file_path=bam_file,
                    reference_file_path=ref_file,
                    output_file_path=output_dir,
                    vcf_file_path=vcf_file)
        start_position = i*each_segment_length
        end_position = (i+1) * each_segment_length + 1000
        p = Process(target=view.parse_region, args=(start_position, end_position, json_out))
        p.start()


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
        default="3",
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
        "--json",
        type=bool,
        default=False,
        help="If true then output will be in a json file in json folder."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/",
        help="If true then output will be in a json file in json folder."
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
        view.test(FLAGS.json)
    else:
        do_parallel(FLAGS.chromosome_name, FLAGS.bam, FLAGS.ref, FLAGS.vcf, FLAGS.json,
                    FLAGS.output_dir, FLAGS.max_threads)
