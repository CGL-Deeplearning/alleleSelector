from collections import defaultdict

"""
CandidateFinder finds possible positions based on mismatches we see in reads.
It parses through each read of a site and finds possible candidate positions.
"""


class CandidateFinder:
    """
    Given reads that align to a site and a pointer to the reference fasta file handler,
    candidate finder finds possible variant candidates of that site.
    """
    def __init__(self, reads, fasta_handler, chromosome_name, window_ref_start_position):
        """
        Initialize a candidate finder object.
        :param reads: Reads that align to the site
        :param fasta_handler: Reference sequence handler
        :param chromosome_name: Chromosome name
        :param window_ref_start_position: Start site of the window
        """
        self.window_ref_start_position = window_ref_start_position
        self.chromosome_name = chromosome_name
        self.fasta_handler = fasta_handler
        self.reads = reads

        # cigar key map based on operation.
        # details: http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        self.cigar_key = {0: "MATCH",
                          1: "INSERT",
                          2: "DELETE",
                          3: "REFSKIP",
                          4: "SOFTCLIP",
                          5: "HARDCLIP",
                          6: "PAD"}

        # for which operation we should increase the reference position
        self.ref_advancement_by_cigar = {"MATCH": True,
                                         "INSERT": False,
                                         "DELETE": True,
                                         "REFSKIP": False,
                                         "SOFTCLIP": False,
                                         "HARDCLIP": False,
                                         "PAD": False}

        # for which cigar operations we should increase the read position
        self.read_advancement_by_cigar = {"MATCH": True,
                                          "INSERT": True,
                                          "DELETE": False,
                                          "REFSKIP": True,
                                          "SOFTCLIP": True,
                                          "HARDCLIP": False,
                                          "PAD": True}

        # the candidates found in reads
        self.candidates = defaultdict(list)

    def get_read_stop_position(self, read):
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

    def parse_reads(self, reads):
        """
        Parse reads to aligned to a site to find variants
        :param reads: Set of reads aligned
        :return:
        """

        for read in reads:
            # check if the mapping quality of the read is above threshold
            if read.mapping_quality > 0:
                self.find_read_candidates(read=read)

    def parse_match(self, alignment_position, read_sequence, ref_sequence, read_name):
        """
        Process a cigar operation that is a match
        :param alignment_position: Position where this match happened
        :param read_sequence: Read sequence
        :param ref_sequence: Reference sequence
        :param read_name: Read ID that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        for i in range(len(read_sequence)):
            if read_sequence[i] != ref_sequence[i]:
                self.candidates[alignment_position].append(read_name)

    def parse_delete(self, alignment_position, length, read_name):
        """
        Process a cigar operation that is a delete
        :param alignment_position: Alignment position
        :param length: Length of the delete
        :param read_name: Read Id that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        start = alignment_position-1
        stop = alignment_position + length + 1
        for i in range(start,stop):
            self.candidates[i].append(read_name)

    def parse_insert(self, alignment_position, length, read_name):
        """
        Process a cigar operation where there is an insert
        :param alignment_position: Position where the insert happpened
        :param length: Length of the insert
        :param read_name: Read Id that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        self.candidates[alignment_position].append(read_name)

    def find_read_candidates(self, read):
        """
        This method finds candidates given a read. We walk through the cigar string to find these candidates.
        :param read: Read from which we need to find the variant candidate positions.
        :return:

        Read candidates use a set data structure to find all positions in the read that has a possible variant.
        """
        read_candidates = set()
        ref_alignment_start = read.reference_start
        ref_alignment_stop = self.get_read_stop_position(read)
        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence
        ref_sequence = self.fasta_handler.get_sequence(chromosome_name=self.chromosome_name,
                                                       start=ref_alignment_start,
                                                       stop=ref_alignment_stop)

        corresponding_ref_sequence = ''
        # DEBUG
        # expanded_cigar used for debugging purpose, giving us the cigar indepent of the number
        # so 3M4I cigar in extended cigar will be MMMIIII.
        expanded_cigar = ''

        # read_index: index of read sequence
        # alignment_position: position that the current cigar tuple aligns to in reference
        read_index = 0
        ref_index = 0
        for cigar in cigar_tuples:
            cigar_code = cigar[0]
            length = cigar[1]
            # get the sequence segments that are effected by this operation
            ref_sequence_segment = ref_sequence[ref_index:ref_index+length]
            read_sequence_segment = read_sequence[read_index:read_index+length]

            # send the cigar tuple to get attributes we got by this operation
            ref_index_increment, read_index_increment, expanded_cigar_segment = \
                self.parse_cigar_tuple(cigar_code=cigar_code,
                                       length=length,
                                       alignment_position=ref_alignment_start+ref_index,
                                       ref_sequence=ref_sequence_segment,
                                       read_sequence=read_sequence_segment,
                                       read_name=read.query_name)

            # increase the read index iterator
            read_index += read_index_increment
            ref_index += ref_index_increment

            # DEBUG
            # for debugging create the extended cigar and reference sequence
            expanded_cigar += expanded_cigar_segment
            corresponding_ref_sequence += ref_sequence_segment

            # print("----------")
            # print(ref_sequence_segment)
            # print(expanded_cigar_segment)
            # print(read_sequence_segment)

        # DEBUG
        # print("--------------------------")
        print(ref_alignment_start)
        print(corresponding_ref_sequence)
        print(expanded_cigar)
        print(read_sequence)
        print(read_candidates)
        print("--------------------------")

    def parse_cigar_tuple(self, cigar_code, length, alignment_position, ref_sequence, read_sequence, read_name):
        """
        Parse through a cigar operation to find possible candidate variant positions in the read
        :param cigar_code: Cigar operation code
        :param length: Length of the operation
        :param alignment_position: Alignment position corresponding to the reference
        :param ref_sequence: Reference sequence
        :param read_sequence: Read sequence
        :param read_name: Read ID
        :return:
        """
        # get what kind of code happened
        cigar_operation = self.cigar_key[cigar_code]
        ref_index_increment = 0
        read_index_increment = 0
        expanded_cigar_string = ''
        # why is this here?
        candidates = set()

        # deal different kinds of operations
        if cigar_operation == "MATCH":
            self.parse_match(alignment_position=alignment_position,
                             read_sequence=read_sequence,
                             ref_sequence=ref_sequence,
                             read_name=read_name)

        elif cigar_operation == "INSERT":
            self.parse_insert(alignment_position=alignment_position, length=length, read_name=read_name)

        elif cigar_operation == "DELETE":
            self.parse_delete(alignment_position=alignment_position, length=length, read_name=read_name)

        elif cigar_operation == "REFSKIP":
            pass

        elif cigar_operation == "SOFTCLIP":
            pass

        elif cigar_operation == "PAD":
            pass

        else:
            raise("INVALID CIGAR CODE: %s"%cigar_code)

        if self.ref_advancement_by_cigar[cigar_operation]:
            ref_index_increment = length

        if self.read_advancement_by_cigar[cigar_operation]:
            read_index_increment = length

        expanded_cigar_string += cigar_operation[0]*length

        return ref_index_increment, read_index_increment, expanded_cigar_string


