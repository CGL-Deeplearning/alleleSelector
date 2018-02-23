import argparse
from modules.IntervalTree import IntervalTree
from modules.TsvHandler import TsvHandler
import csv


def intersect(bed_file_path_query, bed_file_path_reference):
    """
    Takes regions from file a and subsets them based on whether they are entirely contained in regions from file b
    :param bed_file_path_query: path to the BED file that will be subsetted (output contains these regions)
    :param bed_file_path_reference: path to the BED file that contains regions to subset with
    :return:
    """
    tsv_handler_reference = TsvHandler(tsv_file_path=bed_file_path_reference)
    intervals_chromosomal_reference = tsv_handler_reference.get_bed_intervals_by_chromosome(start_offset=1,
                                                                                            universal_offset=-1)

    interval_trees_chromosomal = build_chromosomal_interval_trees(intervals_chromosomal=intervals_chromosomal_reference)

    subset_intervals(bed_file_path_query=bed_file_path_query, interval_trees_chromosomal=interval_trees_chromosomal)


def build_chromosomal_interval_trees(intervals_chromosomal):
    """
    Produce a diictionary of intervals trees, with one tree per chromosome
    :param intervals_chromosomal: dictionary of intervals per chromosome name
    :return: trees_chromosomal
    """
    trees_chromosomal = dict()

    for chromosome_name in intervals_chromosomal:
        intervals = intervals_chromosomal[chromosome_name]
        tree = IntervalTree(intervals)

        trees_chromosomal[chromosome_name] = tree

    return trees_chromosomal


def subset_intervals(bed_file_path_query, interval_trees_chromosomal):
    """
    Test a raw BED file against a dictionary of interval trees based on chromosome, output list of lists (BED format)
    :param bed_file_path_query: path to BED file to subset
    :param interval_trees_chromosomal: dictionary of interval trees per chromosome name
    :return: intervals_subset_chromosomal: a list of lists with first 3 items of each list of the format: chr start stop
    """
    tsv_file = open(bed_file_path_query, 'r')
    reader = csv.reader(tsv_file, delimiter='\t')

    for line in reader:
        chromosome_name, start, stop = line[0:3]
        in_confident = "0"
        # if there is an interval tree for this chromosome, test whether the current query interval is a subset of tree
        if chromosome_name in interval_trees_chromosomal:
            tree = interval_trees_chromosomal[chromosome_name]
            interval = [int(start), int(stop)]

            # if interval is a subset, add it to output
            if tree.contains_interval_subset(interval):
                in_confident = "1"
        print("\t".join(line), end='')
        print("\t"+in_confident)


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--candidate_bed",
        type=str,
        required=True,
        help="Candidate bed file."
    )
    parser.add_argument(
        "--confident_bed",
        type=str,
        required=True,
        help="Bed file containing confident windows."
    )
    FLAGS, unparsed = parser.parse_known_args()
    intersect(bed_file_path_query=FLAGS.candidate_bed, bed_file_path_reference=FLAGS.confident_bed)