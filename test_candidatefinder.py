import argparse
import math
import time
import json
import os
import collections
from multiprocessing import Process
from collections import defaultdict

from modules.CandidateFinder import CandidateFinder
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler
from modules.AlleleFinder import AlleleFinder
from modules.VcfHandler import VCFFileProcessor
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
"""


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

    def add_candidate_to_list(self, alignment_candidates_object):
        """
        Add a candidate to the list
        :param alignment_candidates_object: Candidate object to add
        :return:
        """
        self.all_candidates.append(alignment_candidates_object)

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
    def __init__(self, chromosome_name, bam_file_path, reference_file_path):
        # --- initialize handlers ---
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)

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
        if not os.path.exists("json_output/"):
            os.mkdir("json/")
        smry = open("json/" + "Candidates" + '_' + self.chromosome_name + '_' + str(start) + '_' + str(end) + ".json", 'w')
        smry.write(json.dumps(all_candidate_lists.reprJSON(), cls=ComplexEncoder, indent=4, sort_keys=True))

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
        # candidate_finder.print_windows()

        candidate_windows = candidate_finder.get_candidate_windows()
        all_candidate_lists = AllCandidatesInRegion(self.chromosome_name, start_position, end_position)

        # for each window find list of possible alleles
        for chr_name, window_start, window_end in candidate_windows:
            reference_sequence = self.fasta_handler.get_sequence(chr_name, window_start, window_end+1)
            pileup_columns = self.bam_handler.get_pileupcolumns_aligned_to_a_region(chr_name, window_start, window_end+1)
            allele_finder = AlleleFinder(chr_name, window_start, window_end, pileup_columns, reference_sequence)
            allele_finder.generate_base_dictionaries()
            candidate_list = allele_finder.generate_candidate_allele_list()
            # candidate_list.print_all_candidates()
            all_candidate_lists.add_candidate_to_list(candidate_list)

        if json_out:
            self.write_json(start_position, end_position, all_candidate_lists)

        return all_candidate_lists

    def validate_candidates(self, all_candidate_lists, vcf_dict, summary_file, output_file):
        self.all_pos_in_vcf = defaultdict(int)
        for pos in vcf_dict.keys():
            for record in vcf_dict[pos]:
                self.all_pos_in_vcf[record.pos] = True

        for AlleleCandidateList_object in all_candidate_lists:
            start = AlleleCandidateList_object.window_start
            end = AlleleCandidateList_object.window_end
            ref_seq = AlleleCandidateList_object.reference_sequence

            CandidateInformation_objects = AlleleCandidateList_object.candidate_alleles
            all_alleles = []
            for candidate in CandidateInformation_objects:
                allele_sequence = candidate.allele_sequence
                all_alleles.append(allele_sequence)
            counter = collections.Counter(all_alleles)
            top_2_alleles = counter.most_common(2)
            if len(top_2_alleles) == 2:
                alt_1, alt_2 = top_2_alleles[0], top_2_alleles[1]
            else:
                alt_1 = top_2_alleles[0]
                alt_2 = ''

            output_file.write("Found allele: " + '\n')
            output_file.write(str(start) + " " + str(end) + str(ref_seq) + str(alt_1) + str(alt_2) + '\n')
            for pos in range(start, end+1):
                if pos in vcf_dict.keys():
                    for record in vcf_dict[pos]:
                        self.all_pos_in_vcf[record.pos] = False
                        output_file.write(str(record) + '\n')
            output_file.write("-----------" + '\n')

        total_positions = 0
        total_not_found = 0
        for pos in self.all_pos_in_vcf.keys():
            total_positions += 1
            if self.all_pos_in_vcf[pos] is True:
                total_not_found += 1
                summary_file.write("IN VCF BUT NOT FOUND: " + str(pos) + '\n')

        return total_positions, total_not_found

    def test(self, vcf_file_path, start_position, end_position):
        if not os.path.exists("tmp/"):
            os.mkdir("tmp/")
        smry = open("tmp/" + "summary" + '_' + self.chromosome_name + '_' + str(start_position) + '_' + str(
            end_position) + ".csv", 'w')
        output = open("tmp/" + "output" + '_' + self.chromosome_name + '_' + str(start_position) + '_' + str(
            end_position) + ".txt", 'w')

        start_time = time.time()
        AllCandidatesInRegion_object = self.parse_region(start_position=start_position, end_position=end_position, json_out=False)
        end_time = time.time()

        vcf_handler = VCFFileProcessor(vcf_file_path)
        vcf_handler.populate_dictionary(self.chromosome_name, start_position, end_position, hom_filter=True)

        # get the vcf dictionary of that region
        vcf_dict = vcf_handler.get_variant_dictionary()

        total_positions, total_not_found = self.validate_candidates(AllCandidatesInRegion_object.all_candidates, vcf_dict, smry, output)

        smry.write("SUMMARY: " + '\n')
        smry.write("TIME ELAPSED: " + str(end_time - start_time) + '\n')
        smry.write("TOTAL POSITIONS IN VCF: " + str(total_positions) + '\n')
        smry.write("TOTAL POSITIONS NOT PICKED: " + str(total_not_found) + '\n')


def do_parallel(chr_name, bam_file, ref_file, vcf_file, max_threads=5):
    """
    Split the chromosome into different ranges for parallel processing.
    :param max_threads: Maximum number of threads to create
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
                    reference_file_path=ref_file
                    )
        start_position = i*each_segment_length
        end_position = (i+1) * each_segment_length + 1000
        p = Process(target=view.test, args=(vcf_file, start_position, end_position))
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

    FLAGS, unparsed = parser.parse_known_args()

    view = View(chromosome_name=FLAGS.chromosome_name,
                bam_file_path=FLAGS.bam,
                reference_file_path=FLAGS.ref)

    if FLAGS.test is True:
        view = View(chromosome_name=FLAGS.chromosome_name,
                    bam_file_path=FLAGS.bam,
                    reference_file_path=FLAGS.ref)
        view.test(FLAGS.vcf, start_position=100000, end_position=200000)
    else:
        do_parallel(FLAGS.chromosome_name, FLAGS.bam, FLAGS.ref, FLAGS.vcf, FLAGS.max_threads)