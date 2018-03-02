import os
import sys

SNP = 0
IN = 1


class BedHandler:
    """
    Bed file interface that converts the bed file to a dictionary.
    """
    def __init__(self, bed_file_path):
        """
        Initialize BedHandler object
        :param bed_file_path: Path to a bed file.
        """
        if not os.path.isfile(bed_file_path):
            sys.stderr.write("INVALID BED FILE PATH")
            exit(1)
        self.file_path = bed_file_path
        self.all_bed_records = list()

        # dictionaries to convert the bed file to dictionaries
        self.bed_dictionary = {}
        self.alt_dictionary = {}

        # genotype counts
        self.total_hom = 0
        self.total_het = 0
        self.total_hom_alt = 0

        # convert the bed file to a dictionary
        self.convert_bed_to_dictionary()

    def initialize_dictionary(self, chr, start_pos, end_pos):
        """
        Initialize the dictionaries.
        :param chr: Chromosome name
        :param start_pos: Start position
        :param end_pos: End position
        :return:
        """
        if chr not in self.bed_dictionary.keys():
            self.bed_dictionary[chr] = {}

        for pos in range(start_pos, end_pos):
            if pos not in self.bed_dictionary[chr].keys():
                self.bed_dictionary[chr][pos] = {}
                self.bed_dictionary[chr][pos][SNP] = list()
                self.bed_dictionary[chr][pos][IN] = list()
            if pos not in self.alt_dictionary.keys():
                self.alt_dictionary[pos] = list()

    def _update_genotype_count(self, genotype):
        if genotype == 0:
            self.total_hom += 1
        elif genotype == 1:
            self.total_het += 1
        elif genotype == 2:
            self.total_hom_alt += 1

    def add_record_to_dictionary(self, record):
        """
        Add a record to the dictionary.
        :param record: A bed file record
        :return:
        """
        chr, pos_start, pos_end, ref, alt, genotype = record[0:6]

        self._update_genotype_count(int(genotype))

        pos_start = int(pos_start)
        pos_end = int(pos_end) + 1
        genotype = int(genotype)

        self.initialize_dictionary(chr, pos_start, pos_end)
        for pos in range(pos_start, pos_end):
            # if the alt is already in the dictionary don't add it again.
            if alt in self.alt_dictionary[pos]:
                continue
            if len(ref) == len(alt):
                self.bed_dictionary[chr][pos][SNP].append((ref, alt, genotype))
            else:
                self.bed_dictionary[chr][pos][IN].append((ref, alt, genotype))

            self.alt_dictionary[pos].append(alt)

    def convert_bed_to_dictionary(self):
        """
        Convert a bed file to a dictionary.
        :return:
        """
        with open(self.file_path) as f:
            self.all_bed_records = f.readlines()

        for record in self.all_bed_records:
            if record.rstrip() is None:
                continue
            self.add_record_to_dictionary(tuple(record.rstrip().split('\t')))
