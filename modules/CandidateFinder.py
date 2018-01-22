from collections import defaultdict

"""
CandidateFinder finds possible positions based on edits we see in reads.
It parses through each read of a site and finds possible candidate positions.

Dictionaries it updates:
- candidate_by_read:    records at each position what reads had mismatches  {int -> list}
- coverage:             records coverage of a position                      {int -> int}
- edit_count:           records number of mismatches in a position          {int -> int}

Other data structures:
- candidate_positions (set):    set of positions that of potential candidate variants
- merged_windows (list):        list of windows of potential candidate variants.

Merged_windows is constructed from candidate_positions. If two positions fall within
MERGE_WINDOW_DISTANCE we merge them in a single window.
"""
DEFAULT_MIN_MAP_QUALITY = 5
MERGE_WINDOW_DISTANCE = 1
MERGE_WINDOW_OFFSET = 1
MIN_MISMATCH_THRESHOLD = 5


class CandidateFinder:
    """
    Given reads that align to a site and a pointer to the reference fasta file handler,
    candidate finder finds possible variant candidates_by_read of that site.
    """
    def __init__(self, reads, fasta_handler, chromosome_name, region_start_position, region_end_position):
        """
        Initialize a candidate finder object.
        :param reads: Reads that align to the site
        :param fasta_handler: Reference sequence handler
        :param chromosome_name: Chromosome name
        :param region_start_position: Start position of the region
        :param region_end_position: End position of the region
        """
        self.region_start_position = region_start_position
        self.region_end_position = region_end_position
        self.chromosome_name = chromosome_name
        self.fasta_handler = fasta_handler
        self.reads = reads

        # the store which reads are creating candidates in that position
        self.candidates_by_read = defaultdict(list)
        self.coverage = defaultdict(int)
        self.edit_count = defaultdict(int)
        self.candidate_positions = set()
        # all merged windows
        self.merged_windows = []

    def print_positions(self):
        for pos in sorted(self.candidate_positions):
            print(pos, self.edit_count[pos], self.coverage[pos])

    def print_windows(self):
        print('Total windows: ', len(self.merged_windows))
        for window in self.merged_windows:
            print(window[0], window[1], window[2])

    def get_candidate_windows(self):
        return self.merged_windows

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

    def merge_positions(self):
        """
        Merge candidate positions we discovered in candidate positions set.
        If two positions are within a window they are merged in one window.

        This method updates merged_windows list where each entry is:
        (chromosome_name, start_position, end_position)
        :return:
        """
        start_pos = -1
        end_pos = -1
        for pos in sorted(self.candidate_positions):
            if start_pos == -1 and end_pos == -1:
                start_pos = pos
                end_pos = pos
            # if position is within window from previous position
            elif end_pos + MERGE_WINDOW_DISTANCE >= pos:
                end_pos = pos
            # else call it a window and start a new window
            else:
                self.merged_windows.append((self.chromosome_name, start_pos - MERGE_WINDOW_OFFSET, end_pos + MERGE_WINDOW_OFFSET))
                start_pos = pos
                end_pos = pos

        # if a window was left open inside loop
        if start_pos != -1:
            self.merged_windows.append((self.chromosome_name, start_pos - MERGE_WINDOW_OFFSET, end_pos + MERGE_WINDOW_OFFSET))

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
            self.coverage[alignment_position+i] += 1
            if read_sequence[i] != ref_sequence[i]:
                self.mismatch_count[alignment_position+i] += 1
                self.edit_count[alignment_position+i] += 1
                self.candidates_by_read[alignment_position+i].append(read_name)
                yield alignment_position + i

    def parse_delete(self, alignment_position, length, read_name):
        """
        Process a cigar operation that is a delete
        :param alignment_position: Alignment position
        :param length: Length of the delete
        :param read_name: Read Id that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        start = alignment_position
        stop = alignment_position + length
        for i in range(start, stop):
            self.mismatch_count[i] += 1
            self.coverage[i] += 1
            self.candidates_by_read[i].append(read_name)
            yield i

    def parse_insert(self, alignment_position, length, read_name):
        """
        Process a cigar operation where there is an insert
        :param alignment_position: Position where the insert happpened
        :param length: Length of the insert
        :param read_name: Read Id that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        self.edit_count[alignment_position] += 1
        self.candidates_by_read[alignment_position].append(read_name)
        self.mismatch_count[alignment_position] += 1
        yield alignment_position

    def parse_reads(self, reads):
        """
        Parse reads to aligned to a site to find variants
        :param reads: Set of reads aligned
        :return:
        """
        for read in reads:
            # check if the mapping quality of the read is above threshold
            if read.mapping_quality > DEFAULT_MIN_MAP_QUALITY:
                self.candidate_positions.update(self.find_read_candidates(read=read))

    def find_read_candidates(self, read):
        """
        This method finds candidates given a read. We walk through the cigar string to find these candidates.
        :param read: Read from which we need to find the variant candidate positions.
        :return:

        Read candidates use a set data structure to find all positions in the read that has a possible variant.
        """
        ref_alignment_start = read.reference_start
        ref_alignment_stop = self.get_read_stop_position(read)
        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence
        ref_sequence = self.fasta_handler.get_sequence(chromosome_name=self.chromosome_name,
                                                       start=ref_alignment_start,
                                                       stop=ref_alignment_stop)

        # read_index: index of read sequence
        # ref_index: index of reference sequence
        read_index = 0
        ref_index = 0
        for cigar in cigar_tuples:
            cigar_code = cigar[0]
            length = cigar[1]

            # get the sequence segments that are effected by this operation
            ref_sequence_segment = ref_sequence[ref_index:ref_index+length]
            read_sequence_segment = read_sequence[read_index:read_index+length]

            # send the cigar tuple to get attributes we got by this operation
            ref_index_increment, read_index_increment, candidate_positions = \
                self.parse_cigar_tuple(cigar_code=cigar_code,
                                       length=length,
                                       alignment_position=ref_alignment_start+ref_index,
                                       ref_sequence=ref_sequence_segment,
                                       read_sequence=read_sequence_segment,
                                       read_name=read.query_name)

            # increase the read index iterator
            read_index += read_index_increment
            ref_index += ref_index_increment

            for pos in candidate_positions:
                if self.region_start_position <= pos <= self.region_end_position:
                    if self.mismatch_count[pos] > MIN_MISMATCH_THRESHOLD:
                        yield pos

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

        # why is this here?
        candidate_positions = set()

        # deal different kinds of operations
        if cigar_code == 0:
            # match
            candidate_positions.update(self.parse_match(alignment_position=alignment_position,
                                                        read_sequence=read_sequence,
                                                        ref_sequence=ref_sequence,
                                                        read_name=read_name))
        elif cigar_code == 1:
            # insert
            candidate_positions.update(self.parse_insert(alignment_position=alignment_position,
                                                         length=length,
                                                         read_name=read_name))
            ref_index_increment = 0
        elif cigar_code == 2:
            # delete
            candidate_positions.update(self.parse_delete(alignment_position=alignment_position,
                                                         length=length,
                                                         read_name=read_name))
            read_index_increment = 0
        elif cigar_code == 3:
            # ref skip
            ref_index_increment = 0

        elif cigar_code == 4:
            # soft clip
            ref_index_increment = 0
        elif cigar_code == 5:
            # hard clip
            ref_index_increment = 0
            read_index_increment = 0
        elif cigar_code == 6:
            # pad
            ref_index_increment = 0
            read_index_increment = 0
        else:
            raise("INVALID CIGAR CODE: %s" % cigar_code)

        return ref_index_increment, read_index_increment, candidate_positions


