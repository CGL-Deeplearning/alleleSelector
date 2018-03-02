from collections import defaultdict
from modules.Kaller.ImageChannels import ImageChannels
import numpy as np
from scipy import misc

"""
This script creates an image from a given bed record. 
"""

DEFAULT_MIN_MAP_QUALITY = 5
IMAGE_HEIGHT = 50
IMAGE_WIDTH = 300
IMAGE_BUFFER = 0
CIGAR_MATCH = 0
CIGAR_IN = 1
CIGAR_DEL = 2
MAX_COLOR_VALUE = 254.0
BASE_QUALITY_CAP = 40.0
MAP_QUALITY_CAP = 60.0
MAP_QUALITY_FILTER = 10.0
REF_BAND = 5


class ImageCreatorRGB:
    """
    Create image given a bed record.
    """
    def __init__(self, fasta_handler, chromosome_name, allele_start_position, allele_end_position):
        """
        Initialize image creator object
        :param fasta_handler: Reference file handler
        :param chromosome_name: Chromosome name
        :param allele_start_position: Start position of the allele in question
        :param allele_end_position: End position of the allele in question
        """
        self.chromosome_name = chromosome_name
        self.fasta_handler = fasta_handler

        # the base and the insert dictionary for finding alleles
        self.base_dictionary = {}
        self.insert_dictionary = {}
        self.reference_dictionary = {}

        # supplementary dictionaries and other values
        self.read_id_in_allele_position = list()
        self.longest_insert_in_position = {}
        self.leftmost_alignment_position = allele_start_position
        self.rightmost_alignment_position = allele_end_position
        self.read_rev_dict = {}
        self.read_mq_dict = {}
        # ref position to index projection to handle inserts
        self.ref_to_index_projection = {}

    @staticmethod
    def get_read_stop_position(read):
        """
        Returns the stop position of the reference to where the read stops aligning
        :param read: The read
        :return: stop position of the reference where the read last aligned
        """
        ref_alignment_stop = read.reference_end

        # only find the position if the reference end is fetched as none from pysam API
        if ref_alignment_stop is None:
            positions = read.get_reference_positions()

            # find last entry that isn't None
            i = len(positions) - 1
            ref_alignment_stop = positions[-1]
            while i > 0 and ref_alignment_stop is None:
                i -= 1
                ref_alignment_stop = positions[i]

        return ref_alignment_stop

    def _update_base_dictionary(self, pos, read_id, base, map_quality, base_qualities, direction, cigar_op):
        """
        In base dictionary add attributes to create the image.
        :param pos: Genomic position
        :param read_id: Read id
        :param base: Nucleotide base
        :param map_quality: Mapping quality
        :param base_qualities: Base quality
        :param direction: True if the read  is reverse
        :param cigar_op: CIGAR operation
        :return:
        """
        if read_id not in self.base_dictionary:
            self.base_dictionary[read_id] = {}
        if pos not in self.base_dictionary[read_id]:
            self.base_dictionary[read_id][pos] = []

        self.base_dictionary[read_id][pos] = (base, map_quality, base_qualities, direction, cigar_op)

    def _update_insert_dictionary(self, pos, read_id, base, map_quality, base_qualities, direction, cigar_op):
        """
        In insert dictionary add attributes to create the image.
        :param pos: Genomic position
        :param read_id: Read id
        :param base: Nucleotide base
        :param map_quality: Mapping quality
        :param base_qualities: Array containing base qualities
        :param direction: True if the read  is reverse
        :param cigar_op: CIGAR operation
        :return:
        """
        if read_id not in self.insert_dictionary:
            self.insert_dictionary[read_id] = {}
        if pos not in self.insert_dictionary[read_id]:
            self.insert_dictionary[read_id][pos] = []
        self.insert_dictionary[read_id][pos] = (base, map_quality, base_qualities, direction, cigar_op)

    def _process_match(self, pos, length, read_sequence, read_name, mapping_quality, base_qualities, direction):
        """
        Process a cigar match operation in a read
        :param pos: Starting position of the cigar operation
        :param length: Length of the operation
        :param read_sequence: Read sequence where this operation happens
        :param read_name: Read name
        :param mapping_quality: Mapping quality
        :param base_qualities: Array containing base qualities
        :param direction: True if the read  is reverse
        :return:
        """
        start = pos
        stop = start + length
        for i in range(start, stop):
            read_base = read_sequence[i-pos]
            base_quality = base_qualities[i-pos]
            self._update_base_dictionary(i, read_name, read_base, mapping_quality, base_quality, direction, CIGAR_MATCH)

    def _process_delete(self, pos, length, read_name, mapping_quality, base_qualities, direction):
        """
        Process a cigar delete operation in a read
        :param pos: Starting position of the cigar operation
        :param length: Length of the operation
        :param read_name: Read name
        :param mapping_quality: Mapping quality
        :param base_qualities: Array containing base qualities
        :param direction: True if the read  is reverse
        :return:
        """

        # actual delete position starts one after the anchor
        start = pos
        stop = start + length

        for i in range(start, stop):
            read_base = "*"
            base_quality = 0
            # update the base dictionary
            self._update_base_dictionary(i, read_name, read_base, mapping_quality, base_quality, direction, CIGAR_DEL)

    def _process_insert(self, pos, read_sequence, read_name, mapping_quality, base_qualities, direction):
        """
        Process a cigar delete operation in a read
        :param pos: Starting position of the cigar operation
        :param read_name: Read name
        :param mapping_quality: Mapping quality
        :param base_qualities: Array containing base qualities
        :param direction: True if the read  is reverse
        :return:
        """
        read_bases = read_sequence
        self._update_insert_dictionary(pos, read_name, read_bases, mapping_quality, base_qualities, direction, CIGAR_IN)

        if pos not in self.longest_insert_in_position:
            self.longest_insert_in_position[pos] = 0

        self.longest_insert_in_position[pos] = max(self.longest_insert_in_position[pos], len(read_bases))

    def _update_reference_dictionary(self, position, ref_base):
        """
        Update the reference dictionary
        :param position: Genomic position
        :param ref_base: Reference base at that position
        :return:
        """
        self.reference_dictionary[position] = ref_base

    def parse_cigar_tuple(self, cigar_code, length, alignment_position, read_sequence,
                          read_name, base_qualities, mapping_quality, direction):
        """
        Parse through a cigar operation to find possible candidate variant positions in the read
        :param cigar_code: Cigar operation code
        :param length: Length of the operation
        :param alignment_position: Alignment position corresponding to the reference
        :param read_sequence: Read sequence
        :param read_name: Read ID
        :param base_qualities: Array containing base quality of the read
        :param mapping_quality: Mapping quality of the read
        :param direction: If true then the read is reversed
        :return:

        cigar key map based on operation.
        details: http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        0: "MATCH",
        1: "INSERT",
        2: "DELETE",
        3: "REFSKIP",
        4: "SOFTCLIP",
        5: "HARDCLIP",
        6: "PAD"
        """
        # get what kind of code happened
        ref_index_increment = length
        read_index_increment = length

        # deal different kinds of operations
        if cigar_code == 0:
            # match
            self._process_match(pos=alignment_position,
                                length=length,
                                read_sequence=read_sequence,
                                read_name=read_name,
                                mapping_quality=mapping_quality,
                                base_qualities=base_qualities,
                                direction=direction
                                )
        elif cigar_code == 1:
            # insert
            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            self._process_insert(pos=alignment_position - 1,
                                 read_sequence=read_sequence,
                                 read_name=read_name,
                                 mapping_quality=mapping_quality,
                                 base_qualities=base_qualities,
                                 direction=direction
                                 )
            ref_index_increment = 0
        elif cigar_code == 2 or cigar_code == 3:
            # delete or ref_skip
            self._process_delete(pos=alignment_position,
                                 length=length,
                                 read_name=read_name,
                                 mapping_quality=mapping_quality,
                                 base_qualities=base_qualities,
                                 direction=direction
                                 )
            read_index_increment = 0
        elif cigar_code == 4:
            # soft clip
            ref_index_increment = 0
            # print("CIGAR CODE ERROR SC")
        elif cigar_code == 5:
            # hard clip
            ref_index_increment = 0
            read_index_increment = 0
            # print("CIGAR CODE ERROR HC")
        elif cigar_code == 6:
            # pad
            ref_index_increment = 0
            read_index_increment = 0
            # print("CIGAR CODE ERROR PAD")
        else:
            raise("INVALID CIGAR CODE: %s" % cigar_code)

        return ref_index_increment, read_index_increment

    def _update_image_bounderies(self, read_start, read_end):
        """
        Update the leftmost and rightmost alignment positions.
        :param read_start: Read alignment start
        :param read_end: Read alignment end
        :return:
        """
        self.leftmost_alignment_position = min(self.leftmost_alignment_position, read_start)
        self.rightmost_alignment_position = max(self.rightmost_alignment_position, read_end)

    def _process_read(self, read):
        """
        Process a read that aligns to the allele position
        :param read:
        :return:
        """
        ref_alignment_start = read.reference_start
        ref_alignment_stop = self.get_read_stop_position(read)
        self._update_image_bounderies(ref_alignment_start, ref_alignment_stop)

        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence

        base_qualities = read.query_qualities
        mapping_quality = read.mapping_quality
        direction = read.is_reverse
        self.read_mq_dict[read.query_name] = mapping_quality
        self.read_rev_dict[read.query_name] = direction

        read_index = 0
        ref_index = 0

        for cigar in cigar_tuples:
            cigar_code = cigar[0]
            length = cigar[1]
            # get the sequence segments that are effected by this operation
            read_sequence_segment = read_sequence[read_index:read_index+length]
            base_quality_segment = base_qualities[read_index:read_index+length]

            # send the cigar tuple to get attributes we got by this operation
            ref_index_increment, read_index_increment = \
                self.parse_cigar_tuple(cigar_code=cigar_code,
                                       length=length,
                                       alignment_position=ref_alignment_start+ref_index,
                                       read_sequence=read_sequence_segment,
                                       read_name=read.query_name,
                                       base_qualities=base_quality_segment,
                                       mapping_quality=mapping_quality,
                                       direction=direction)

            # increase the read index iterator
            read_index += read_index_increment
            ref_index += ref_index_increment

    def process_reads(self, reads):
        """
        Parse reads to aligned to a site to find variants
        :param reads: Set of reads aligned
        :return:
        """
        i = 0
        for read in reads:
            if i > IMAGE_HEIGHT-REF_BAND:
                break
            # check if the mapping quality of the read is above threshold
            if read.mapping_quality > DEFAULT_MIN_MAP_QUALITY:
                self.read_id_in_allele_position.append(read.query_name)
                self._process_read(read=read)
                i += 1

    def get_start_and_end_positions(self, position):
        """
        If leftmost and rightmost positions are out of window then find the left and right boundaries.
        :param position: Position where the allele resides
        :return:
        """
        distance_array = defaultdict(int)
        distance_array[self.leftmost_alignment_position] = 0

        for pos in range(self.leftmost_alignment_position, self.rightmost_alignment_position):
            if pos in self.longest_insert_in_position.keys():
                distance_array[pos + 1] += self.longest_insert_in_position[pos]

            distance_array[pos] = distance_array[pos-1] + 1

        if self.rightmost_alignment_position - self.leftmost_alignment_position + 1 <= IMAGE_WIDTH:
            return self.leftmost_alignment_position, self.rightmost_alignment_position

        left_side = right_side = int((IMAGE_WIDTH-IMAGE_BUFFER) / 2)

        left_val = max(0, distance_array[position] - left_side)
        right_val = min(len(distance_array.keys()), distance_array[position] + right_side)
        left_pos, right_pos = position, position

        for pos in sorted(distance_array.keys()):
            if distance_array[pos] < left_val:
                left_pos = pos
            if distance_array[pos] < right_val:
                right_pos = pos

        return left_pos, right_pos

    def get_reference_row(self, start_pos, end_pos):
        """
        Get the reference row.
        :param start_pos: Start position of the reference.
        :param end_pos: End position of the reference
        :return:
        """
        ref_row = [ImageChannels.get_empty_rgb_channels() for i in range(IMAGE_WIDTH)]
        for i in range(start_pos, end_pos):
            base = self.reference_dictionary[i]
            if self.ref_to_index_projection[i] < IMAGE_WIDTH:
                ref_row[self.ref_to_index_projection[i]] = ImageChannels.get_ref_channels_rgb(base)

            if i in self.longest_insert_in_position:
                for j in range(self.longest_insert_in_position[i]):
                    if self.ref_to_index_projection[i] + j + 1 < IMAGE_WIDTH:
                        ref_row[self.ref_to_index_projection[i]+j+1] = ImageChannels.get_ref_channels_rgb('*')

        return ref_row

    def _if_read_supports_alt(self, read_id, position, alt):
        """
        Check if read supports the alt allele in question
        :param read_id: Read id
        :param position: Position of the alt allele
        :param alt: The alt allele
        :return:
        """
        read_base = ''
        if read_id in self.base_dictionary and position in self.base_dictionary[read_id]:
            read_base += self.base_dictionary[read_id][position][0]
        if len(alt) > 1 and read_id in self.insert_dictionary and position in self.insert_dictionary[read_id]:
            read_base += self.insert_dictionary[read_id][position][0]

        if read_base == alt:
            return True

        return False

    def get_read_row(self, read_id, left_pos, right_pos, alts, alt_position):
        """
        Convert a read to an image row
        :param read_id: Read id
        :param left_pos: Leftmost position of the image
        :param right_pos: Rightmost position of the image
        :param alts: Alternate alleles
        :param alt_position: Alternate allele position
        :return:
        """
        image_row = [ImageChannels.get_empty_rgb_channels() for i in range(IMAGE_WIDTH)]
        is_supporting = False
        is_match = False

        for alt in alts:
            is_supporting = is_supporting or self._if_read_supports_alt(read_id, alt_position, alt)

        for pos in range(left_pos, right_pos):
            if read_id in self.base_dictionary and pos in self.base_dictionary[read_id]:
                base, map_q, base_q, is_rev, cigar_op = self.base_dictionary[read_id][pos]
                if base == self.reference_dictionary[pos]:
                    is_match = True

                attribute_tuple = (base, base_q, map_q, is_rev, is_match, is_supporting)
                # create channels for the base in that position
                channels = ImageChannels.get_channels_only_rgb(attribute_tuple, self.reference_dictionary[pos])
                if self.ref_to_index_projection[pos] < IMAGE_WIDTH:
                    image_row[self.ref_to_index_projection[pos]] = channels

            is_match = False
            if read_id in self.insert_dictionary and pos in self.insert_dictionary[read_id]:
                # if there is an insert
                bases, map_q, base_qs, is_rev, cigar_op = self.insert_dictionary[read_id][pos]
                row_index = self.ref_to_index_projection[pos] + 1

                # for each base of the insert
                for i, base in enumerate(bases):
                    attribute_tuple = (base, base_qs[i], map_q, is_rev, is_match, is_supporting)
                    channels = ImageChannels.get_channels_only_rgb(attribute_tuple, self.reference_dictionary[pos])
                    if row_index < IMAGE_WIDTH:
                        image_row[row_index] = channels
                    row_index += 1

                # if the insert is not the longest insert of that position
                if len(bases) < self.longest_insert_in_position[pos]:
                    for i in range(self.longest_insert_in_position[pos]-len(bases)):
                        attribute_tuple = ('*', 0, map_q, is_rev, is_match, is_supporting)
                        channels = ImageChannels.get_channels_only_rgb(attribute_tuple, '')
                        if row_index < IMAGE_WIDTH:
                            image_row[row_index] = channels
                        row_index += 1
            # if there is an insert at this position but not in the read
            elif pos in self.longest_insert_in_position:
                row_index = self.ref_to_index_projection[pos] + 1
                for i in range(self.longest_insert_in_position[pos]):
                    attribute_tuple = ('*', 0, self.read_mq_dict[read_id], self.read_rev_dict[read_id], False, is_supporting)
                    channels = ImageChannels.get_channels_only_rgb(attribute_tuple, '')
                    if row_index < IMAGE_WIDTH:
                        image_row[row_index] = channels
                    row_index += 1

        return image_row

    def generate_read_pileups(self, left_pos, right_pos, position, alts):
        """
        Generate rows for the reads that align to an allele position.
        :param left_pos: Leftmost position in the image
        :param right_pos: Rightmost position in the image
        :param position: Alternate allele position
        :param alts: Alternate alleles in question
        :return:
        """
        all_read_ids = self.read_id_in_allele_position
        image_rows = list()
        for read_id in all_read_ids:
            image_rows.append(self.get_read_row(read_id, left_pos, right_pos, alts, position))
            if len(image_rows) == IMAGE_HEIGHT - REF_BAND:
                break

        return image_rows

    def project_ref_positions(self, left_pos, right_pos):
        """
        Calculate the index where each reference position should go as inserts distort their positions.
        :param left_pos: Leftmost genomic position in the image
        :param right_pos: Rightmost genomic position in the image
        :return:
        """
        index = 0
        for pos in range(left_pos, right_pos):
            self.ref_to_index_projection[pos] = index
            if pos in self.longest_insert_in_position.keys():
                index += self.longest_insert_in_position[pos]
            index += 1

    def _update_ref_sequence(self, start_position, end_position):
        """
        Update the reference sequence
        :param start_position: Start position
        :param end_position: End position
        :return:
        """
        ref_sequence = self.fasta_handler.get_sequence(chromosome_name=self.chromosome_name,
                                                       start=start_position,
                                                       stop=end_position)

        for i, ref_base in enumerate(ref_sequence):
            self._update_reference_dictionary(start_position + i, ref_base)

    def generate_image(self, position, alts):
        """
        Generate an image given allele position
        :param position: Allele position
        :param alts: Alternate alleles
        :return:
        """
        left_pos, right_pos = self.get_start_and_end_positions(position)
        self._update_ref_sequence(left_pos, right_pos)
        self.project_ref_positions(left_pos, right_pos)
        img_data = list()
        for i in range(REF_BAND):
            img_data.append(self.get_reference_row(left_pos, right_pos))

        img_data = img_data + self.generate_read_pileups(left_pos, right_pos, position, [alts])

        while len(img_data) < IMAGE_HEIGHT:
            image_row = [ImageChannels.get_empty_rgb_channels() for i in range(IMAGE_WIDTH)]
            img_data.append(image_row)

        image_array = np.array(img_data).astype(np.uint8)

        return image_array
