import modules.BamHandler as bam_processor
import argparse
import sys
import pysam
from modules.CandidateFinder import CandidateFinder
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler

class View:
    def __init__(self, chromosome_name, bam_file_path, reference_file_path, window_size): #, output_dir, window_size, window_cutoff, coverage_cutoff, map_quality_cutoff, vcf_quality_cutoff, coverage_threshold, threads):
        # --- initialize handlers ---
        self.candidate_finder = CandidateFinder()
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)

        # --- initialize parameters ---
        self.chromosome_name = chromosome_name
        self.window_size = window_size
        self.current_position = 0

    def iterate_windows(self):
        '''
        Given a user defined window size, iterate through windows in the BAM and collect candidate alleles
        '''

        reference_sequence = self.fasta_handler.get_sequence(chromosome_name=self.chromosome_name,
                                                             start=self.current_position,
                                                             stop=self.current_position+self.window_size)

        reads = self.bam_handler.get_reads(chromosome_name=self.chromosome_name,
                                           start=self.current_position,
                                           stop=self.current_position+self.window_size)

        print(reference_sequence)
        print(reads)

        # self.candidates = self.candidate_finder.parse_reads()

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
        "--chromosome_name",
        type=str,
        default="3",
        help="Desired chromosome number E.g.: 3"
    )
    parser.add_argument(
        "--window_size",
        type=int,
        default=1000,
        help="Window size of query region."
    )
    # parser.add_argument(
    #     "--output_dir",
    #     type=str,
    #     default="output/",
    #     help="Name of output directory"
    # )
    # parser.add_argument(
    #     "--parallel",
    #     type=bool,
    #     default=False,
    #     help="Option to run threaded pileup generator"
    # )
    # parser.add_argument(
    #     "--max_threads",
    #     type=int,
    #     default=5,
    #     help="Number of maximum threads for this region."
    # )

    FLAGS, unparsed = parser.parse_known_args()

    view = View(chromosome_name = FLAGS.chromosome_name,
                bam_file_path = FLAGS.bam,
                reference_file_path = FLAGS.ref,
                window_size = FLAGS.window_size)
                # output_dir = FLAGS.output_dir,
                # threads = FLAGS.max_threads)

    view.iterate_windows()

# usage example:
# python3 main.py --bam /Users/saureous/data/chr3_200k.bam --ref /Users/saureous/data/chr3.fa --chromosome_name chr3 --window_size 1000