from VcfHandler import VCFFileProcessor
from BedHandler import BedHandler
from FileManager import FileManager
from timeit import default_timer as timer
import argparse
import multiprocessing
import os.path
from os import remove
from os.path import join
import csv
import sys
from multiprocessing import Pool

'''
EXAMPLE USAGE:

python3 modules/FilterTest.py --vcf /users/saureous/data/NA12878_S1.genome.vcf.gz --allele_bed /users/saureous/data/bed_alleles/Label_chr3_120000_140010.bed --confident_bed /users/saureous/data/ConfidentRegions.bed.gz

'''

# VCF record indexes
REF,ALT,GT = 0,1,2

# BED record indexes
CHROMOSOME_NAME,START,STOP,REF_BED,ALT_BED,GT_BED = 0,1,2,3,4,5

# validation stats
TP,FN,FP,UNVALIDATED_ALLELES = 0,1,2,3

# unvalidated alleles
POS,VARIANT = 0,1

class View:
    """
    I/O manager for FilterTest
    """
    def __init__(self, vcf_file_path, confident_bed_file_path, allele_bed_file_path, output_directory_path):
        self.confident_bed_file_path = confident_bed_file_path
        self.allele_bed_file_path = allele_bed_file_path
        self.vcf_handler = VCFFileProcessor(vcf_file_path)
        self.bed_handler_confident = BedHandler(confident_bed_file_path)
        self.bed_handler_allele = BedHandler(allele_bed_file_path)
        self.bed_handler_confident_alleles = None
        self.bed_handler_confident_vcf = None
        self.bed_handler_vcf = None
        self.tester = FilterTest()

        # print(self.allele_bed_file_path)
        self.chromosome_name, self.bed_start, self.bed_stop = self.get_region_from_file_path()

        self.output_filename = self.chromosome_name + "_" + str(self.bed_start) + "_" + str(self.bed_stop) + ".txt"
        self.output_bed_filename = "confident_" + self.output_filename[:-4] + ".bed"
        self.output_vcf_filename = "vcf_" + self.output_filename[:-4] + ".bed"
        self.output_confident_vcf_filename = "confident_vcf_" + self.output_filename[:-4] + ".bed"

        self.output_directory_path = output_directory_path
        if not os.path.exists(output_directory_path):
            os.mkdir(output_directory_path)

        self.output_path = join(self.output_directory_path, self.output_filename)
        self.output_bed_path = join(self.output_directory_path, self.output_bed_filename)
        self.output_vcf_bed_path = join(self.output_directory_path, self.output_vcf_filename)
        self.output_confident_vcf_bed_path = join(self.output_directory_path, self.output_confident_vcf_filename)

        self.confident_alleles_file = None
        self.output_file = None
        self.length = None
        self.unintersected_length = None

    def intersect_alleles_with_confident_bed(self, deallocate=True):
        """
        Use BedHandler to intersect the labeled candidates with the confident VCF regions. Save intersected BED.
        :param deallocate:
        :return:
        """

        self.bed_handler_confident_alleles = self.bed_handler_allele.intersect(self.bed_handler_confident)
        self.bed_handler_confident_alleles.save(self.output_bed_path)

        self.length = len(self.bed_handler_confident_alleles)
        self.unintersected_length = len(self.bed_handler_allele)

    def intersect_vcf_with_confident_bed(self):
        self.vcf_handler.save_positional_vcf_as_bed(chromosome_name=self.chromosome_name,
                                                    output_path_name=self.output_vcf_bed_path)

        self.bed_handler_vcf = BedHandler(self.output_vcf_bed_path)

        self.bed_handler_confident_vcf = self.bed_handler_vcf.intersect(self.bed_handler_confident)
        self.bed_handler_confident_vcf.save(self.output_confident_vcf_bed_path)

        variants = self.vcf_handler.read_positional_vcf_from_bed(self.output_confident_vcf_bed_path)

        remove(self.output_vcf_bed_path)
        remove(self.output_confident_vcf_bed_path)

        del self.bed_handler_allele
        del self.bed_handler_confident
        del self.bed_handler_confident_alleles
        del self.bed_handler_confident_vcf
        del self.bed_handler_vcf

        return variants

    def parse_region(self):
        """
        Parse the input data from end to end, return the output stats
        :return:
        """
        self.intersect_alleles_with_confident_bed()

        bed_file = open(self.output_bed_path, 'r')
        confident_alleles = csv.reader(bed_file, delimiter='\t')

        stats = [0, 0, 0, []]

        # if there are records contained in the intersected BED
        if self.length > 0:

            # get dictionary of variant records for full region
            self.vcf_handler.populate_dictionary(contig=self.chromosome_name,
                                                 start_pos=self.bed_start,
                                                 end_pos=self.bed_stop,
                                                 hom_filter=True)

            positional_variants = self.intersect_vcf_with_confident_bed()

            stats = self.tester.validate_alleles(positional_variants, confident_alleles)

        bed_file.close()
        remove(self.output_bed_path)

        print("COMPLETED: ", self.chromosome_name, self.bed_start, self.bed_stop)

        return [self.chromosome_name, self.bed_start, self.bed_stop] + list(stats)

    def get_region_from_file_path(self):
        """
        Use the systematic file name of the input BED allele file to get genomic coordinates
        :return: chromosome_name, start, stop
        """
        # Label_chr1_97600000_97800000.bed
        filename = os.path.basename(self.allele_bed_file_path).split('.')[0]
        # print(filename)
        filename_words = filename.split('_')
        # print(filename_words)
        chromosome_name = filename_words[1]
        start = int(filename_words[2])
        stop = int(filename_words[3])

        return chromosome_name, start, stop


class FilterTest:
    """
    Performs summary statistics on BED files
    """
    def validate_alleles(self, variants, alleles):
        """
        Find the TP,FP,FN for a set of alleles given the VCF. Store any False Negatives for debugging/printing
        :param variants: positional vcf
        :param alleles: bed_handler object for the alleles of interest
        :return: n_true_positive, n_false_negative, n_false_positive, unvalidated_positions (a tuple of (pos,variant))
        """

        validated_vcf_positions = {key: 0 for key in variants.keys()}
        n_false_negative = 0
        n_false_positive = 0
        n_true_positive = 0

        for e, entry in enumerate(alleles):
            # print(entry)

            genotype = int(entry[GT_BED])
            allele_start = int(entry[START])

            # test if the site has any true variants (not Hom)
            if self._is_supported(genotype):
                if allele_start in validated_vcf_positions:
                    validated_vcf_positions[allele_start] = 1
                    n_true_positive += 1
            else:
                n_false_positive += 1

        result = self.count_false_negatives(validated_vcf_positions=validated_vcf_positions,
                                            variants=variants)

        n_false_negative, unvalidated_positions = result

        return n_true_positive, n_false_negative, n_false_positive, unvalidated_positions

    @staticmethod
    def _is_supported(genotype):
        """
        Test whether a candidate allele has a corresponding VCF variant. i.e.: is the allele labeled as het/hom_alt ?
        :param genotype: the label/genotype code for an allele
        :return: boolean T/F stating whether the position is validated
        """
        supported = False

        if genotype > 0:
            supported = True

        return supported

    @staticmethod
    def count_false_negatives(validated_vcf_positions, variants):
        """
        Iterate a dictionary of VCF positions that have been marked as supported by a candidate and count any unmarked.
        :param validated_vcf_positions:
        :param variants:
        :return:
        """
        unvalidated_positions = list()

        n_false_negative = 0

        for p,pos in enumerate(validated_vcf_positions):
            # print(p)
            if validated_vcf_positions[pos] == 0:
                n_false_negative += 1
                unvalidated_positions.append([pos, variants[pos]])

        return n_false_negative, unvalidated_positions


def print_unvalidated_positions(output_file, unvalidated_alleles):
    for allele in unvalidated_alleles:
        output_file.write("\nWARNING: Unsupported VCF position: %d\n"%allele[POS])
        output_file.write("\tRecord: %s\n"%str(allele[VARIANT]))
    output_file.write("\n")


def print_region_info(output_file, chromosome_name, bed_start, bed_stop, n_false_negative, n_false_positive, n_true_positive):
    output_file.write("BED region:     %s %d %d\n"%(chromosome_name, bed_start, bed_stop))
    output_file.write("False Negative: %d\n"%n_false_negative)
    output_file.write("False Positive: %d\n"%n_false_positive)
    output_file.write("True Positive:  %d\n"%n_true_positive)


def print_region_tsv_line(tsv_writer, chromosome_name, bed_start, bed_stop, n_false_negative, n_false_positive, n_true_positive):
    tsv_writer.writerow([chromosome_name,
                         bed_start,
                         bed_stop,
                         n_false_negative,
                         n_false_positive,
                         n_true_positive])


def parallel_run(vcf_file_path, confident_bed_file_path, allele_bed_file_path, output_dir_path, return_dict):
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
                allele_bed_file_path=allele_bed_file_path,
                output_directory_path=output_dir_path)

    data = view.parse_region()
    key = data[1]               # unique key identifying chromosome coordinate
    return_dict[key] = data


def chromosome_level_parallelization(vcf_file_path, confident_bed_file_path, allele_bed_directory_path,
                                     max_threads, output_dir_path):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param vcf_file_path: Chromosome name
    :param confident_bed_file_path: Bam file path
    :param allele_bed_directory_path: Ref file path
    :return:
    """
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []

    # entire set of files
    file_manager = FileManager()
    allele_bed_file_paths = file_manager.get_file_paths_from_directory(allele_bed_directory_path)

    for file_path in allele_bed_file_paths:
        if file_path.endswith(".bed"):
            # parse window of the segment. Use a 1000 overlap for corner cases.
            args = (vcf_file_path, confident_bed_file_path, file_path, output_dir_path, return_dict)

            p = multiprocessing.Process(target=parallel_run, args=args)
            jobs.append(p)
            p.start()

            while True:
                if len(multiprocessing.active_children()) < max_threads:
                    break

    for proc in jobs:
        proc.join()

    write_results_dictionary_to_file(output_dir_path=output_dir_path, results=return_dict)


def write_results_dictionary_to_file(output_dir_path, results):
    output_filename = join(output_dir_path, "chromosome_output.txt")
    output_file = open(output_filename, 'w')
    tsv_writer = csv.writer(output_file, delimiter='\t')

    total_true_positive = 0
    total_false_negative = 0
    total_false_positive = 0

    all_unvalidated_alleles = list()

    for key in sorted(results.keys()):
        chromosome_name, \
        bed_start, \
        bed_stop, \
        n_true_positive, \
        n_false_negative, \
        n_false_positive, \
        unvalidated_alleles = results[key]

        total_true_positive += n_true_positive
        total_false_negative += n_false_negative
        total_false_positive += n_false_positive

        print_region_tsv_line(tsv_writer=tsv_writer,
                              chromosome_name=chromosome_name,
                              bed_start=bed_start,
                              bed_stop=bed_stop,
                              n_false_negative=n_false_negative,
                              n_false_positive=n_false_positive,
                              n_true_positive=n_true_positive)

        all_unvalidated_alleles.append(unvalidated_alleles)

    for allele_list in all_unvalidated_alleles:
        print_unvalidated_positions(output_file=output_file,
                                    unvalidated_alleles=allele_list)

    output_file.close()


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
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/validator_output/",
        # required=False,
        help="Where to save results"
    )

    FLAGS, unparsed = parser.parse_known_args()

    chromosome_level_parallelization(vcf_file_path=FLAGS.vcf,
                                     confident_bed_file_path=FLAGS.confident_bed,
                                     allele_bed_directory_path=FLAGS.allele_dir,
                                     max_threads=FLAGS.max_threads,
                                     output_dir_path=FLAGS.output_dir)


