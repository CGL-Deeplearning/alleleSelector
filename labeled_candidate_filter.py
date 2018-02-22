import argparse
import csv


def filter_bed(bed_file_path):
    """
    Test a raw BED file against a dictionary of interval trees based on chromosome, output list of lists (BED format)
    :param bed_file_path_query: path to BED file to subset
    :param interval_trees_chromosomal: dictionary of interval trees per chromosome name
    :return: intervals_subset_chromosomal: a list of lists with first 3 items of each list of the format: chr start stop
    """
    tsv_file = open(bed_file_path, 'r')
    reader = csv.reader(tsv_file, delimiter='\t')

    for line in reader:
        rec_qual = int(line[6])
        rec_filter = line[7]
        in_confident = int(line[8])
        if rec_qual >=60 and rec_filter == 'PASS' and in_confident == 1:
            print('\t'.join(line))


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

    FLAGS, unparsed = parser.parse_known_args()
    filter_bed(bed_file_path=FLAGS.candidate_bed)