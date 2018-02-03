from VcfHandler import VCFFileProcessor
from BedHandler import BedHandler
from FileManager import FileManager
import pybedtools
import argparse
import multiprocessing
import os.path
import sys

'''
EXAMPLE USAGE:

python3 modules/FilterTest.py --vcf /users/saureous/data/NA12878_S1.genome.vcf.gz --allele_bed /users/saureous/data/bed_alleles/Label_chr3_120000_140010.bed --confident_bed /users/saureous/data/ConfidentRegions.bed.gz

'''

# VCF record indexes
REF,ALT,GT = 0,1,2

# BED record indexes
CHROMOSOME_NAME,START,STOP,REF_BED,ALT_BED,GT_BED = 0,1,2,3,4,5

class View:
    """
    I/O manager for FilterTest
    """
    def __init__(self, vcf_file_path, confident_bed_file_path, allele_bed_file_path):
        self.confident_bed_file_path = confident_bed_file_path
        self.allele_bed_file_path = allele_bed_file_path
        self.vcf_handler = VCFFileProcessor(vcf_file_path)
        self.bed_handler_confident = BedHandler(confident_bed_file_path)
        self.bed_handler_allele = BedHandler(allele_bed_file_path)
        self.tester = FilterTest()

    def get_allele_validation_stats(self, deallocate=False):
        confident_alleles = self.bed_handler_allele.intersect(self.bed_handler_confident)

        chromosome_name, bed_start, bed_stop = self.get_region_from_file_path()

        if deallocate:
            del self.bed_handler_allele
            del self.bed_handler_confident

        length = len(confident_alleles)

        # if there are records contained in the intersected BED
        if length > 0:
            chromosome_name = confident_alleles[0].chrom
            # allele_start = confident_alleles[0].start
            # allele_stop = confident_alleles[length-1].stop

            print("BED region:                ", chromosome_name, bed_start, bed_stop)
            print("Items before intersection: ", len(self.bed_handler_allele))
            print("Items after intersection:  ", len(confident_alleles))

            # get dictionary of variant records for full region
            self.vcf_handler.populate_dictionary(contig=chromosome_name,
                                                 start_pos=bed_start,
                                                 end_pos=bed_stop,
                                                 hom_filter=True)

            # get separate positional variant dictionaries for IN, DEL, and SNP
            positional_variants = self.vcf_handler.get_variant_dictionary()

            self.tester.validate_alleles(positional_variants, confident_alleles)

            confident_alleles.save("confident.bed")

        else:
            self.print_empty_region(chromosome_name, bed_start, bed_stop, len(confident_alleles))

    def get_region_from_file_path(self):
        # Label_chr1_97600000_97800000.bed

        filename = os.path.basename(self.allele_bed_file_path).split('.')[0]
        filename_words = filename.split('_')
        chromosome_name = filename_words[1]
        start = int(filename_words[2])
        stop = int(filename_words[3])

        return chromosome_name, start, stop

    def print_empty_region(self,chromosome_name,bed_start,bed_stop,confident_length):
        print("BED region:                ", chromosome_name, bed_start, bed_stop)
        print("Items before intersection: ", len(self.bed_handler_allele))
        print("Items after intersection:  ", confident_length)
        print("False Negative: ", None)
        print("False Positive: ", None)
        print("True Positive:  ", None)


class FilterTest:
    """
    Performs summary statistics on BED files
    """
    def validate_alleles(self, variants, alleles):
        """
        Find the TP,FP,FN for a set of alleles given the VCF
        :param variants:
        :param alleles:
        :return:
        """
        validated_vcf_positions = {key: 0 for key in variants.keys()}
        n_false_negative = 0
        n_false_positive = 0
        n_true_positive = 0

        for entry in alleles:
            genotype = int(entry[GT_BED])
            allele_start = int(entry[START])

            # test if the site has any true variants (not Hom)
            if self._is_supported(genotype):
                if allele_start in validated_vcf_positions:
                    validated_vcf_positions[allele_start] = 1
                    n_true_positive += 1
            else:
                n_false_positive += 1

        n_false_negative = self.count_false_negatives(validated_vcf_positions=validated_vcf_positions, variants=variants)
        self.print_outcome(n_false_negative, n_false_positive, n_true_positive)

    def _is_supported(self,genotype):
        supported = False

        if genotype > 0:
            supported = True

        return supported

    @staticmethod
    def count_false_negatives(validated_vcf_positions, variants):
        n_false_negative = 0

        for pos in validated_vcf_positions:
            if validated_vcf_positions[pos] == 0:
                n_false_negative += 1
                print("\nWARNING: Unsupported VCF position: %d", pos)
                print("\tRecord: ", variants[pos])

        return n_false_negative

    @staticmethod
    def print_outcome(n_false_negative, n_false_positive, n_true_positive):
        print("False Negative: ", n_false_negative)
        print("False Positive: ", n_false_positive)
        print("True Positive:  ", n_true_positive)


def parallel_run(vcf_file_path, confident_bed_file_path, allele_bed_file_path):
    """
    Run this method in parallel
    :param vcf_file_path: Chromosome name
    :param confident_bed_file_path: Bam file path
    :param allele_bed_file_path: Ref file path
    :return:
    """

    # create a view object
    view = View(vcf_file_path=vcf_file_path,
                confident_bed_file_path=confident_bed_file_path,
                allele_bed_file_path=allele_bed_file_path)

    # return the results
    view.get_allele_validation_stats()


def chromosome_level_parallelization(vcf_file_path, confident_bed_file_path, allele_bed_directory_path, max_threads):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param vcf_file_path: Chromosome name
    :param confident_bed_file_path: Bam file path
    :param allele_bed_directory_path: Ref file path
    :return:
    """
    # entire set of files
    file_manager = FileManager()
    allele_bed_file_paths = file_manager.get_file_paths_from_directory(allele_bed_directory_path)

    for file_path in allele_bed_file_paths:
        # parse window of the segment. Use a 1000 overlap for corner cases.
        args = (vcf_file_path, confident_bed_file_path, file_path)

        p = multiprocessing.Process(target=parallel_run, args=args)
        p.start()

        while True:
            if len(multiprocessing.active_children()) < max_threads:
                break


def chromosome_level_iteration(vcf_file_path, confident_bed_file_path, allele_bed_directory_path):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param vcf_file_path: Chromosome name
    :param confident_bed_file_path: Bam file path
    :param allele_bed_directory_path: Ref file path
    :return:
    """
    # entire set of files
    file_manager = FileManager()
    allele_bed_file_paths = file_manager.get_file_paths_from_directory(allele_bed_directory_path)

    for file_path in allele_bed_file_paths:
        # parse window of the segment. Use a 1000 overlap for corner cases.
        view = View(vcf_file_path=vcf_file_path, confident_bed_file_path=confident_bed_file_path, allele_bed_file_path=file_path)
        view.get_allele_validation_stats()


if __name__ == '__main__':
    """
    Processes arguments and performs tasks to generate the pileup.
    """
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF file containing all known variant records."
    )
    parser.add_argument(
        "--confident_bed",
        type=str,
        required=True,
        help="BED file containing ranges of confident variant calls for the VCF."
    )
    parser.add_argument(
        "--allele_dir",
        type=str,
        required=True,
        help="Directory containing set of chunked BED files to be parsed."
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=6,
        required=False,
        help="Number of threads to use in parallel."
    )
    parser.add_argument(
        "--test",
        type=bool,
        default=False,
        required=False,
        help="Run without parallelization"
    )

    FLAGS, unparsed = parser.parse_known_args()

    if FLAGS.test:
        chromosome_level_iteration(vcf_file_path=FLAGS.vcf,
                                   confident_bed_file_path=FLAGS.confident_bed,
                                   allele_bed_directory_path=FLAGS.allele_dir)
    else:
        chromosome_level_parallelization(vcf_file_path=FLAGS.vcf,
                                         confident_bed_file_path=FLAGS.confident_bed,
                                         allele_bed_directory_path=FLAGS.allele_dir,
                                         max_threads=FLAGS.max_threads)


