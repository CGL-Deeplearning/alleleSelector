import argparse
import math
import os
import time
from modules.CandidateFinder import CandidateFinder
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler
from  modules.vcf_handler import VCFFileProcessor
from collections import defaultdict
from multiprocessing import Process

"""
alignmentPolish finds possible variant sites in given bam file.

It requires three parameters:
- bam_file_path: path to a bam file
- reference_file_path: path to a reference file

Creates:
- CandidateFinder object that contains windows of possible variants.
"""

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

        self.test_positions = defaultdict(int)

    def parse_window(self, start_position, end_position):
        """
        Find possible candidate windows.
        """
        reads = self.bam_handler.get_reads(chromosome_name=self.chromosome_name,
                                           start=start_position,
                                           stop=end_position)

        candidate_finder = CandidateFinder(reads=reads,
                                           fasta_handler=self.fasta_handler,
                                           chromosome_name=self.chromosome_name,
                                           window_start_position=start_position,
                                           window_end_position=end_position)

        # parse reads to find candidate positions
        candidate_finder.parse_reads(reads=reads)
        # print(sorted(candidate_finder.candidate_positions))
        # merge candidate positions to windows
        candidate_finder.merge_positions()
        # print(candidate_finder.merged_windows)
        # print the windows we got
        return candidate_finder.merged_windows

    def test(self, vcf_file_path, start_position, end_position):
        if not os.path.exists("tmp/"):
            os.mkdir("tmp/")
        smry = open("tmp/" + "summary" + '_' + self.chromosome_name + '_' + str(start_position) + '_' + str(end_position) + ".csv", 'w')

        start_time = time.time()
        merged_windows = self.parse_window(start_position=start_position, end_position=end_position)
        end_time = time.time()

        total_positions = 0
        for windows in merged_windows:
            start_pos = windows[1]
            end_pos = windows[2]
            for pos in range(start_pos, end_pos+1):
                self.test_positions[pos] = 1
                total_positions += 1

        vcf_handler = VCFFileProcessor(vcf_file_path)
        vcf_handler.populate_dictionary(self.chromosome_name, start_position, end_position, hom_filter=True)

        # get the vcf dictionary of that region
        vcf_dict = vcf_handler.get_variant_dictionary()
        not_found_count = 0
        total_positions_in_vcf = 0
        for pos in vcf_dict.keys():
            for rec in vcf_dict[pos]:
                rec_len = len(rec.ref)
                rec_start = pos
                rec_end = rec_start + rec_len
                for i in range(rec_start, rec_end):
                    total_positions_in_vcf += 1
                    if self.test_positions[i] != 1:
                        smry.write(str(i) + '\n')
                        not_found_count += 1

        smry.write("TIME ELAPSED: " + str(end_time - start_time) + '\n')
        smry.write("TOTAL POSITIONS: " + str(total_positions)  +'\n')
        smry.write("TOTAL POSITIONS IN VCF: " + str(total_positions_in_vcf) + '\n')
        smry.write("TOTAL POSITIONS NOT PICKED: " + str(not_found_count) + '\n')


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
                    reference_file_path=ref_file)
        p = Process(target=view.test, args=(vcf_file, i*each_segment_length, (i+1)*each_segment_length + 1000))
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
        help="VCF file for testing."
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

    FLAGS, unparsed = parser.parse_known_args()

    if FLAGS.test is True:
        view = View(chromosome_name=FLAGS.chromosome_name,
                    bam_file_path=FLAGS.bam,
                    reference_file_path=FLAGS.ref)
        view.test(FLAGS.vcf, start_position=100000, end_position=400000)
    else:
        do_parallel(FLAGS.chromosome_name, FLAGS.bam, FLAGS.ref, FLAGS.vcf, FLAGS.max_threads)

# usage example:
# python3 main.py --bam /Users/saureous/data/chr3_200k.bam --ref /Users/saureous/data/chr3.fa --chromosome_name chr3 --window_size 1000