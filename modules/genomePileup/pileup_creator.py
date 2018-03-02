"""
Implemented by: Kishwar Shafin
Date: 02/01/2018
"""

import pysam
import sys
import numpy as np
from PIL import Image
from scipy import misc

"""
This script creates pileup images given a vcf record, bam alignment file and reference fasta file.

imageChannels: Handles how many channels to create for each base and their structure

"""

MAX_COLOR_VALUE = 254.0
BASE_QUALITY_CAP = 40.0
# this is used to define qualities for '*' where map quality is undefined
BASE_QUALITY_MIN = 10.0
MAP_QUALITY_CAP = 60.0
MAP_QUALITY_FILTER = 10.0
MATCH_CIGAR_CODE = 0
INSERT_CIGAR_CODE = 1
DELETE_CIGAR_CODE = 2

class imageChannels:
    """
    Handles how many channels to create for each base and their way of construction.
    """

    def __init__(self, pileup_attributes, ref_base, is_supporting):
        """
        Initialize a base with it's attributes
        :param pileup_attributes: Attributes of a pileup base
        :param ref_base: Reference base corresponding to that pileup base
        """
        self.pileup_base = pileup_attributes[0]
        self.base_qual = pileup_attributes[1]
        self.map_qual = pileup_attributes[2]
        self.cigar_code = pileup_attributes[3]
        self.is_rev = pileup_attributes[4]
        self.ref_base = ref_base
        self.is_match = True if self.ref_base == self.pileup_base else False
        self.is_supporting = is_supporting

    @staticmethod
    def get_base_color(base):
        """
        Get color based on a base.
        - Uses different band of the same channel.
        :param base:
        :return:
        """
        if base == 'A':
            return 250.0
        if base == 'C':
            return 100.0
        if base == 'G':
            return 180.0
        if base == 'T':
            return 30.0
        if base == '*' or 'N':
            return 5.0

    @staticmethod
    def get_base_quality_color(base_quality):
        """
        Get a color spectrum given base quality
        :param base_quality: value of base quality
        :return:
        """
        c_q = min(base_quality, BASE_QUALITY_CAP)
        color = MAX_COLOR_VALUE * c_q / BASE_QUALITY_CAP
        return color

    @staticmethod
    def get_map_quality_color(map_quality):
        """
        Get a color spectrum given mapping quality
        :param map_quality: value of mapping quality
        :return:
        """
        c_q = min(map_quality, MAP_QUALITY_CAP)
        color = MAX_COLOR_VALUE * c_q / MAP_QUALITY_CAP
        return color

    @staticmethod
    def get_strand_color(is_rev):
        """
        Get color for forward and reverse reads
        :param is_rev: True if read is reversed
        :return:
        """
        if is_rev is True:
            return 240
        else:
            return 70

    @staticmethod
    def get_match_ref_color(is_match):
        """
        Get color for base matching to reference
        :param is_match: If true, base matches to reference
        :return:
        """
        if is_match is True:
            return MAX_COLOR_VALUE * 0.2
        else:
            return MAX_COLOR_VALUE * 1.0

    @staticmethod
    def get_alt_support_color(is_in_support):
        """
        Get support color
        :param is_in_support: Boolean value of support
        :return:
        """
        if is_in_support is True:
            return MAX_COLOR_VALUE * 1.0
        else:
            return MAX_COLOR_VALUE * 0.6

    @staticmethod
    def get_cigar_color(cigar_code):
        """
        ***NOT USED YET***
        :param is_in_support:
        :return:
        """
        if cigar_code == 0:
            return MAX_COLOR_VALUE
        if cigar_code == 1:
            return MAX_COLOR_VALUE * 0.6
        if cigar_code == 2:
            return MAX_COLOR_VALUE * 0.3

    @staticmethod
    def get_empty_channels():
        """
        Get empty channel values
        :return:
        """
        return [0, 0, 0, 0, 0, 0, 0]

    def get_channels(self):
        """
        Get a bases's channel construction
        :return: [color spectrum of channels based on base attributes]
        """
        base_color = self.get_base_color(self.pileup_base)
        base_quality_color = imageChannels.get_base_quality_color(self.base_qual)
        map_quality_color = imageChannels.get_map_quality_color(self.map_qual)
        strand_color = imageChannels.get_strand_color(self.is_rev)
        match_color = imageChannels.get_match_ref_color(self.is_match)
        support_color = imageChannels.get_alt_support_color(self.is_supporting)
        cigar_color = imageChannels.get_cigar_color(self.cigar_code)
        return [base_color, base_quality_color, map_quality_color, strand_color, match_color, support_color, cigar_color]

    @staticmethod
    def get_channels_for_ref(base):
        """
        Get a reference bases's channel construction
        :param base: Reference base
        :return: [color spectrum of channels based on some default values]
        """
        base_color = imageChannels.get_base_color(base)
        base_quality_color = imageChannels.get_base_quality_color(60)
        map_quality_color = imageChannels.get_map_quality_color(60)
        strand_color = imageChannels.get_strand_color(is_rev=False)
        match_color = imageChannels.get_match_ref_color(is_match=True)
        support_color = imageChannels.get_alt_support_color(is_in_support=True)
        cigar_color = imageChannels.get_cigar_color(MATCH_CIGAR_CODE)
        return [base_color, base_quality_color, map_quality_color, strand_color, match_color, support_color, cigar_color]

    # RGB image creator
    # ---ONLY USED FOR TESTING--- #
    @staticmethod
    def get_empty_rgb_channels():
        return [0, 0, 0, 255]

    @staticmethod
    def get_color_for_base_rgb(ref, base):
        if ref == base and ref != '*':
            return 255, 255, 255
        elif base == 'A':
            return 255, 0, 0
        elif base == 'C':
            return 255, 255, 0
        elif base == 'T':
            return 0, 0, 255
        elif base == 'G':
            return 0, 255, 0
        else:
            return 255, 0, 255

    def get_channels_only_rgb(self):
        base_color = self.get_base_color(self.pileup_base)
        base_quality_color = imageChannels.get_base_quality_color(self.base_qual)
        map_quality_color = imageChannels.get_map_quality_color(self.map_qual)
        strand_color = imageChannels.get_strand_color(self.is_rev)
        match_color = imageChannels.get_match_ref_color(self.is_match)
        support_color = imageChannels.get_alt_support_color(self.is_supporting)
        r, g, b = self.get_color_for_base_rgb(self.ref_base, self.pileup_base)

        return [r, g, b, support_color]

    @staticmethod
    def get_channels_for_ref_only_rgb(base):
        base_color = imageChannels.get_base_color(base)
        base_quality_color = imageChannels.get_base_quality_color(60)
        map_quality_color = imageChannels.get_map_quality_color(60)
        strand_color = imageChannels.get_strand_color(is_rev=False)
        get_match_color = imageChannels.get_match_ref_color(is_match=True)
        r, g, b = imageChannels.get_color_for_base_rgb('', base)
        support_color = imageChannels.get_alt_support_color(is_in_support=True)

        return [r, g, b, support_color]


class PileupProcessor:
    """
    Processes a pileup around a positoin
    """
    def __init__(self, ref_object, pileupcolumns, contig, pos, genotype, alt):
        """
        Initialize PileupProcessor object with required dictionaries
        :param ref_object: pysam FastaFile object that contains the reference
        :param pileupcolumns: pysam AlignmentFIle.pileup object that contains reads of a position
        :param contig: Contig (ex chr3)
        :param pos: Position in contig
        :param genotype: Genotype
        :param alt: Alternate allele
        """
        self.ref_object = ref_object
        self.pileupcolumns = pileupcolumns
        self.contig = contig
        self.pos = pos
        self.genotype = genotype
        self.alt = alt
        # [genomic_position] = [max_insert_length]
        self.insert_length_dictionary = {} # used
        # [read_id] = {{genomic_position}->base}
        self.read_dictionary = {} # used
        # [read_id] = {{genomic_position}->insert_bases}
        self.read_insert_dictionary = {} #used
        # List of Read ids in a genomic position
        self.reads_aligned_to_pos = {}
        # genomic_position_1, genomic_position_2...
        self.position_list = [] # used
        self.leftmost_genomic_position = -1
        self.rightmost_genomic_position = -1
        self.genomic_position_projection = {}
        self.reference_base_projection = {}
        self.ref_sequence = ''
        self.process_pileup()
        self.project_genomic_positions()

    def project_genomic_positions(self):
        """
        Generate reference sequence with inserts based on two dictionaries
        :return:
        """
        if self.leftmost_genomic_position < 0:
            self.leftmost_genomic_position = 0
        if self.rightmost_genomic_position < 0:
            self.rightmost_genomic_position = 0

        # get the reference sequence
        ref_seq, error_val = self.ref_object.get_ref_of_region(self.contig,
                                                    ":"+str(self.leftmost_genomic_position+1)+ "-"
                                                    + str(self.rightmost_genomic_position+1))
        if error_val == 1:
            print("ERROR IN FETCHING REFERENCE: ", self.contig, self.pos, self.alt, self.genotype)

        ref_seq_with_insert = ''
        idx = 0
        for i in range(self.leftmost_genomic_position, self.rightmost_genomic_position+1):
            # projection of genomic position
            self.genomic_position_projection[i] = idx
            # get reference of that position
            self.reference_base_projection[i] = ref_seq[i-self.leftmost_genomic_position]
            ref_seq_with_insert += ref_seq[i-self.leftmost_genomic_position]
            idx += 1
            # if genomic position has insert
            if i in self.insert_length_dictionary:
                # append inserted characters to the reference
                ref_seq_with_insert += (self.insert_length_dictionary[i] * '*')
                idx += self.insert_length_dictionary[i]
        # set the reference sequence
        self.ref_sequence = ref_seq_with_insert

        # return index
        return idx

    def length_of_region(self):
        """
        Return the length of the sequence from left to rightmost genomic position.
        :return:
        """
        length = 0
        for i in range(self.leftmost_genomic_position, self.rightmost_genomic_position):
            length += 1
            if i in self.insert_length_dictionary:
                length += self.insert_length_dictionary[i]
        return length

    def initialize_dictionaries(self, genomic_position, read_id, is_insert):
        """
        Initialize all the dictionaries for a specific position
        :param genomic_position: Genomic position of interest
        :param read_id: Read id for which dictionaries should be initialized
        :param is_insert: If the position is an insert
        :return:
        """
        if self.leftmost_genomic_position < 0 or genomic_position < self.leftmost_genomic_position:
            self.leftmost_genomic_position = genomic_position
        if self.rightmost_genomic_position < 0 or genomic_position > self.rightmost_genomic_position:
            self.rightmost_genomic_position = genomic_position

        if genomic_position not in self.reads_aligned_to_pos:
            self.reads_aligned_to_pos[genomic_position] = []

        if read_id not in self.read_dictionary:
            self.read_dictionary[read_id] = {}
            self.read_dictionary[read_id][genomic_position] =''

        if is_insert:
            if genomic_position not in self.insert_length_dictionary:
                self.insert_length_dictionary[genomic_position] = 0
            if read_id not in self.read_insert_dictionary:
                self.read_insert_dictionary[read_id] = {}
                self.read_insert_dictionary[read_id][genomic_position] =''

    def save_info_of_a_position(self, genomic_position, read_id, base, base_qual, map_qual, is_rev, cigar_code, is_in):
        """
        Given the attributes of a base at a position
        :param genomic_position: Genomic position
        :param read_id: Read id of a read
        :param base: Base at the position
        :param base_qual: Base quality
        :param map_qual: Map quality
        :param is_rev: If read is reversed
        :param is_in:
        :return:
        """
        self.initialize_dictionaries(genomic_position, read_id, is_in)

        if is_in is False:
            self.read_dictionary[read_id][genomic_position] = (base, base_qual, map_qual, cigar_code, is_rev)
        else:
            self.read_insert_dictionary[read_id][genomic_position] = (base, base_qual, map_qual, cigar_code, is_rev)
            self.insert_length_dictionary[genomic_position] = max(self.insert_length_dictionary[genomic_position],
                                                                  len(base))

    @staticmethod
    def get_attributes_to_save_indel( pileupcolumn, pileupread):
        insert_start = pileupread.query_position + 1
        insert_end = insert_start + pileupread.indel

        return pileupcolumn.pos, \
               pileupread.alignment.query_name, \
               pileupread.alignment.query_sequence[insert_start:insert_end], \
               pileupread.alignment.query_qualities[insert_start:insert_end], \
               pileupread.alignment.mapping_quality, \
               pileupread.alignment.is_reverse, \
               INSERT_CIGAR_CODE # CIGAR OPERATION IS INSERT

    @staticmethod
    def get_attributes_to_save(pileupcolumn, pileupread):
        if pileupread.is_del:
            return pileupcolumn.pos, \
                   pileupread.alignment.query_name,\
                   '*', \
                   BASE_QUALITY_MIN,\
                   pileupread.alignment.mapping_quality, \
                   pileupread.alignment.is_reverse, \
                   DELETE_CIGAR_CODE # CIGAR OPERATION DELETE
        else:
            return pileupcolumn.pos, \
                   pileupread.alignment.query_name, \
                   pileupread.alignment.query_sequence[pileupread.query_position],  \
                   pileupread.alignment.query_qualities[pileupread.query_position], \
                   pileupread.alignment.mapping_quality, \
                   pileupread.alignment.is_reverse, \
                   MATCH_CIGAR_CODE  # CIGAR OPERATION MATCH

    @staticmethod
    def save_image_as_png(pileup_array, save_dir, file_name):
        pileupArray2d = pileup_array.reshape((pileup_array.shape[0], -1))
        misc.imsave(save_dir + file_name + ".png", pileupArray2d, format="PNG")

    def process_pileup(self):
        for pileupcolumn in self.pileupcolumns:
            self.position_list.append(pileupcolumn.pos)
            self.reads_aligned_to_pos[pileupcolumn.pos] = []
            for pileupread in pileupcolumn.pileups:
                self.reads_aligned_to_pos[pileupcolumn.pos].append(pileupread.alignment.query_name)

                if pileupread.indel > 0:
                    gen_pos, read_id, base, base_qual, map_qual, is_rev, cigar_code = \
                        self.get_attributes_to_save_indel(pileupcolumn, pileupread)
                    self.save_info_of_a_position(gen_pos, read_id, base, base_qual, map_qual, is_rev, cigar_code, is_in=True)

                gen_pos, read_id, base, base_qual, map_qual, is_rev, cigar_code = \
                    self.get_attributes_to_save(pileupcolumn, pileupread)
                self.save_info_of_a_position(gen_pos, read_id, base, base_qual, map_qual, is_rev, cigar_code, is_in=False)

    def create_text_pileup(self, query_pos):
        left_most_pos = -1
        for read_id in self.reads_aligned_to_pos[query_pos]:
            read_list = []
            aligned_positions = sorted(self.read_dictionary[read_id].keys())
            if left_most_pos < 0:
                left_most_pos = aligned_positions[0]
            left_most_pos = min(left_most_pos, aligned_positions[0])
            inserts_in_between = sum(self.insert_length_dictionary[val] if val in self.insert_length_dictionary.keys() else 0 for val in range(left_most_pos, aligned_positions[0]))
            padding = aligned_positions[0] - left_most_pos + inserts_in_between
            for pad in range(padding):
                read_list.append(' ')
            for pos in aligned_positions:
                read_list.append((self.read_dictionary[read_id][pos][0]))
                if pos in self.insert_length_dictionary.keys() and self.insert_length_dictionary[pos] > 0:
                    this_has_insert = read_id in self.read_insert_dictionary and pos in self.read_insert_dictionary[read_id]
                    inserted_bases = 0
                    if this_has_insert is True:
                        for base in self.read_insert_dictionary[read_id][pos][0]:
                            read_list.append(base)
                            inserted_bases += 1
                    for i in range(inserted_bases, self.insert_length_dictionary[pos]):
                        read_list.append('*')

            print(''.join(read_list))

    def check_for_support(self, read_id, ref, alt, poi):
        genomic_start_position = poi
        genomic_end_position = poi + len(ref)
        allele = ''
        for pos in range(genomic_start_position, genomic_end_position):
            if pos in self.read_dictionary[read_id]:
                allele += self.read_dictionary[read_id][pos][0]
            if len(alt) > 1 and read_id in self.read_insert_dictionary and pos in self.read_insert_dictionary[read_id]:
                allele += self.read_insert_dictionary[read_id][pos][0]
        allele = allele.replace('*', '')
        alt = alt.replace('*', '')
        if allele == alt:
            return True
        return False

    def get_row(self, read_id, poi, ref, alt):
        read_list = {}
        read_insert_list = {}
        is_supporting = self.check_for_support(read_id, ref, alt, poi)

        aligned_positions = sorted(self.read_dictionary[read_id].keys())
        for pos in aligned_positions:
            read_list[pos] = []
            read_list[pos].append(self.read_dictionary[read_id][pos])

            if pos in self.insert_length_dictionary.keys() and self.insert_length_dictionary[pos] > 0:
                read_insert_list[pos] = []
                inserted_bases = 0
                if read_id in self.read_insert_dictionary and pos in self.read_insert_dictionary[read_id]:
                    inserted_bases = len(self.read_insert_dictionary[read_id][pos][0])
                    read_insert_list[pos].append(self.read_insert_dictionary[read_id][pos])

                for i in range(inserted_bases, self.insert_length_dictionary[pos]):
                    read_attribute_tuple = ('*', [BASE_QUALITY_CAP], self.read_dictionary[read_id][pos][2],
                                            INSERT_CIGAR_CODE, self.read_dictionary[read_id][pos][4])
                    read_insert_list[pos].append(read_attribute_tuple)
        return read_list, read_insert_list, is_supporting

    def get_reference_row_rgb(self, image_width):
        image_row = [imageChannels.get_empty_rgb_channels() for i in range(image_width)]
        for i in range(0, min(len(self.ref_sequence), image_width)):
            image_row[i] = imageChannels.get_channels_for_ref_only_rgb(self.ref_sequence[i])
        return image_row

    def create_image_rgb(self, query_pos, image_height, image_width, ref_band, ref, alt):
        in_support = 0
        not_in_support = 0
        whole_image = []
        for i in range(ref_band):
            whole_image.append(self.get_reference_row_rgb(image_width))

        for read_id in self.reads_aligned_to_pos[query_pos]:
            row_list, row_insert_list, is_supporting = self.get_row(read_id, query_pos, ref, alt)
            in_support = in_support + 1 if is_supporting is True else in_support
            not_in_support = not_in_support + 1 if is_supporting is False else not_in_support

            image_row = [imageChannels.get_empty_rgb_channels() for i in range(image_width)]

            for position in sorted(row_list):
                imagechannels_object = imageChannels(row_list[position][0], self.reference_base_projection[position],
                                                     is_supporting)
                if self.genomic_position_projection[position] < image_width:
                    image_row[self.genomic_position_projection[position]] = imagechannels_object.get_channels_only_rgb()

                if position in row_insert_list.keys():
                    insert_ref = 0
                    for bases in row_insert_list[position]:
                        for base_idx in range(len(bases[0])):
                            insert_ref += 1
                            attribute_tuple = (bases[0][base_idx], bases[1][base_idx], bases[2], bases[3])
                            imagechannels_object = imageChannels(attribute_tuple, '*', is_supporting)
                            if self.genomic_position_projection[position] + insert_ref < image_width:
                                image_row[self.genomic_position_projection[position] + insert_ref] = \
                                    imagechannels_object.get_channels_only_rgb()

            whole_image.append(image_row)

        empty_rows_to_add = image_height - len(whole_image)
        for i in range(empty_rows_to_add):
            whole_image.append([imageChannels.get_empty_rgb_channels() for i in range(image_width)])

        image_array = np.array(whole_image).astype(np.uint8)
        img = Image.fromarray(image_array)
        return img, in_support, not_in_support

    # TEST FIVE CHANNELS
    def get_reference_row(self, image_width):
        image_row = [imageChannels.get_empty_channels() for i in range(image_width)]
        for i in range(0, min(len(self.ref_sequence), image_width)):
            image_row[i] = imageChannels.get_channels_for_ref(self.ref_sequence[i])
        return image_row

    @staticmethod
    def add_empty_rows(image, empty_rows_to_add, image_width):
        for i in range(empty_rows_to_add):
            image.append([imageChannels.get_empty_channels() for i in range(image_width)])
        return image

    def create_image(self, query_pos, image_height, image_width, ref_band, ref, alt):
        whole_image = []
        for i in range(ref_band):
            whole_image.append(self.get_reference_row(image_width))

        for read_id in self.reads_aligned_to_pos[query_pos]:
            row_list, row_insert_list, is_supporting = self.get_row(read_id, query_pos, ref, alt)

            image_row = [imageChannels.get_empty_channels() for i in range(image_width)]

            filter_row = False
            for position in sorted(row_list):
                if row_list[position][0][2] < MAP_QUALITY_FILTER:
                    filter_row = True
                    break

                imagechannels_object = imageChannels(row_list[position][0], self.reference_base_projection[position],
                                                     is_supporting)

                if self.genomic_position_projection[position] < image_width:
                    image_row[self.genomic_position_projection[position]] = imagechannels_object.get_channels()

                if position in row_insert_list.keys():
                    insert_ref = 0
                    for bases in row_insert_list[position]:
                        for base_idx in range(len(bases[0])):
                            insert_ref += 1
                            attribute_tuple = (bases[0][base_idx], bases[1][base_idx], bases[2], bases[3], bases[4])
                            imagechannels_object = imageChannels(attribute_tuple, '*', is_supporting)

                            if self.genomic_position_projection[position] + insert_ref < image_width:
                                image_row[self.genomic_position_projection[position] + insert_ref] = \
                                    imagechannels_object.get_channels()

            if filter_row is False and len(whole_image) < image_height:
                whole_image.append(image_row)

        whole_image = self.add_empty_rows(whole_image, image_height - len(whole_image), image_width)

        image_array = np.array(whole_image).astype(np.uint8)
        return image_array, image_array.shape

