import argparse
import csv

def filter_or_not(alt1, filter1, alt2, filter2, in_confident):
    if alt2 == '.' and (filter1 == 'PASS' or filter1 == '.') and in_confident == 1:
        return 1
    elif (filter1 == 'PASS' or filter1 == '.') and (filter2 == 'PASS' or filter2 == '.') and in_confident == 1:
        return 1

def filter_bed(bed_file_path):
    """
    Test a raw BED file against a dictionary of interval trees based on chromosome, output list of lists (BED format)
    :param bed_file_path_query: path to BED file to subset
    :param interval_trees_chromosomal: dictionary of interval trees per chromosome name
    :return: intervals_subset_chromosomal: a list of lists with first 3 items of each list of the format: chr start stop
    """
    tsv_file = open(bed_file_path, 'r')
    reader = csv.reader(tsv_file, delimiter='\t')
    # chr3	111205	111205	C	T	.	SUB	1	TruthSensitivityTranche99.00to99.90	0	.	1
    for line in reader:
        alt1 = line[3]
        alt2 = line[4]
        filter1 = line[8]
        filter2 = line[10]
        in_confident = int(line[11])

        if filter_or_not(alt1, filter1, alt2, filter2, in_confident):
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