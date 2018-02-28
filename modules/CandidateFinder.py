from collections import defaultdict
import operator

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
MERGE_WINDOW_DISTANCE = 0
MERGE_WINDOW_OFFSET = 0
MIN_MISMATCH_PERCENT_THRESHOLD = 4
MIN_COVERAGE_THRESHOLD = 10
PLOIDY = 2


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
        self.mismatch_count = defaultdict(int)
        self.edit_count = defaultdict(int)
        self.candidate_positions = set()

        # the base and the insert dictionary for finding alleles
        self.snp_dictionary = {}
        self.insert_dictionary = {}
        self.reference_dictionary = {}

    def print_positions(self):
        for pos in sorted(self.candidate_positions):
            print(pos, self.edit_count[pos], self.coverage[pos])

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

    def _update_snp_dictionary(self, pos, allele):
        """
        Update the snp dictionary with allele frequency in a position
        :param pos: Position
        :param allele: Allele at that position
        :return:
        """
        if pos not in self.snp_dictionary:
            self.snp_dictionary[pos] = {}
        if allele not in self.snp_dictionary[pos]:
            self.snp_dictionary[pos][allele] = 0

        self.snp_dictionary[pos][allele] += 1

    def _update_insert_dictionary(self, pos, allele):
        """
        Update insert dictionary with allele frequency in a position
        :param pos: Position
        :param allele: Insert allele at that position
        :return:
        """
        if pos not in self.insert_dictionary:
            self.insert_dictionary[pos] = {}

        if allele not in self.insert_dictionary[pos]:
            self.insert_dictionary[pos][allele] = 0

        self.insert_dictionary[pos][allele] += 1

    def parse_match(self, alignment_position, length, read_sequence, ref_sequence, read_name):
        """
        Process a cigar operation that is a match
        :param alignment_position: Position where this match happened
        :param read_sequence: Read sequence
        :param ref_sequence: Reference sequence
        :param length: Length of the operation
        :param read_name: Read ID that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        start = alignment_position
        stop = start + length
        for i in range(start, stop):
            self.coverage[i] += 1
            if read_sequence[i-alignment_position] != ref_sequence[i-alignment_position]:
                # it's a true mismatch to the reference
                # update the base dictionary
                self._update_snp_dictionary(i, read_sequence[i-alignment_position])
                # increase mismatch count
                self.mismatch_count[i] += 1
                # increase the edit count
                self.edit_count[i] += 1
                # update the candidates by read list
                self.candidates_by_read[i].append(read_name)
                # yield the position to be added to the candidate positions
                yield i

    def parse_delete(self, alignment_position, length, read_name):
        """
        Process a cigar operation that is a delete
        :param alignment_position: Alignment position
        :param length: Length of the delete
        :param read_name: Read Id that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        # actual delete position starts one after the anchor
        start = alignment_position + 1
        stop = start + length

        for i in range(start, stop):
            # update the base dictionary
            self._update_snp_dictionary(i, "*")
            # increase the mismatch count
            self.mismatch_count[i] += 1
            # increase the coverage
            self.coverage[i] += 1
            # append the read name to candidate read list
            self.candidates_by_read[i].append(read_name)
            yield i

    def parse_insert(self, alignment_position, read_sequence, read_name):
        """
        Process a cigar operation where there is an insert
        :param alignment_position: Position where the insert happened
        :param read_sequence: The insert read sequence
        :param read_name: Read Id that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        self.edit_count[alignment_position] += 1
        self.candidates_by_read[alignment_position].append(read_name)
        self.mismatch_count[alignment_position] += 1
        self._update_insert_dictionary(alignment_position, self.reference_dictionary[alignment_position] + read_sequence)
        yield alignment_position

    def _select_alleles(self, position):
        """
        Given a genomic position, naively find the top 2 represented sequences and return two lists of alleles
        :param position: Position where to find the alleles
        :return: top_2_alleles: the most frequent two alleles
        """
        # keep only the top n most frequent alleles
        n = PLOIDY
        snp_alleles = list()
        in_alleles = list()

        if position in self.snp_dictionary:
            snp_alleles = self.snp_dictionary[position]
            snp_alleles = sorted(snp_alleles.items(), key=operator.itemgetter(1), reverse=True)

        if position in self.insert_dictionary:
            in_alleles = self.insert_dictionary[position]
            in_alleles = sorted(in_alleles.items(), key=operator.itemgetter(1), reverse=True)

        if len(snp_alleles) > n:
            snp_alleles = snp_alleles[:n]
        if len(in_alleles) > n:
            in_alleles = in_alleles[:n]

        return snp_alleles, in_alleles

    def parse_reads_and_select_candidates(self, reads):
        """
        Parse reads to aligned to a site to find variants
        :param reads: Set of reads aligned
        :return:
        """
        for read in reads:
            # check if the mapping quality of the read is above threshold
            if read.mapping_quality > DEFAULT_MIN_MAP_QUALITY:
                self.candidate_positions.update(self.find_read_candidates(read=read))

        filtered_candidate_positions = list()
        for pos in self.candidate_positions:
            if self.region_start_position <= pos <= self.region_end_position:
                percent_mismatch = int((self.mismatch_count[pos] * 100) / self.coverage[pos])
                if self.coverage[pos] > MIN_COVERAGE_THRESHOLD and \
                        percent_mismatch > MIN_MISMATCH_PERCENT_THRESHOLD:
                    filtered_candidate_positions.append(pos)

        self.candidate_positions = filtered_candidate_positions

        selected_allele_list = []
        for pos in self.candidate_positions:
            n_inserts, n_snps = self._select_alleles(position=pos)
            ref_base = self.reference_dictionary[pos]
            selected_allele_list.append((pos, ref_base, n_snps, n_inserts))

        return selected_allele_list

    def _update_reference_dictionary(self, position, ref_base):
        """
        Update the reference dictionary
        :param position: Genomic position
        :param ref_base: Reference base at that position
        :return:
        """
        self.reference_dictionary[position] = ref_base

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

        for i, ref_base in enumerate(ref_sequence):
            self._update_reference_dictionary(ref_alignment_start + i, ref_base)

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
                                                        length=length,
                                                        read_sequence=read_sequence,
                                                        ref_sequence=ref_sequence,
                                                        read_name=read_name))
        elif cigar_code == 1:
            # insert
            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            candidate_positions.update(self.parse_insert(alignment_position=alignment_position-1,
                                                         read_sequence=read_sequence,
                                                         read_name=read_name))
            ref_index_increment = 0
        elif cigar_code == 2 or cigar_code == 3:
            # delete or ref_skip
            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            candidate_positions.update(self.parse_delete(alignment_position=alignment_position-1,
                                                         length=length,
                                                         read_name=read_name))
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

        return ref_index_increment, read_index_increment, candidate_positions


