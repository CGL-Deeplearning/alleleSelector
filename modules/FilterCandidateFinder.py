from collections import defaultdict
import operator
import time
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
MIN_MISMATCH_THRESHOLD = 3
MIN_MISMATCH_PERCENT_THRESHOLD = 4
MIN_COVERAGE_THRESHOLD = 10
PLOIDY = 8
MATCH_ALLELE = 0
MISMATCH_ALLELE = 1
INSERT_ALLELE = 2
DELETE_ALLELE = 3


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
        self.delete_dictionary = {}
        self.insert_dictionary = {}
        self.positional_allele_dictionary = {}
        self.read_allele_dictionary = {}
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

    def _update_read_allele_dictionary(self, pos, allele, type):
        if pos not in self.read_allele_dictionary:
            self.read_allele_dictionary[pos] = {}
        if (allele, type) not in self.read_allele_dictionary[pos]:
            self.read_allele_dictionary[pos][(allele, type)] = 0

        self.read_allele_dictionary[pos][(allele, type)] += 1

    def _update_positional_allele_dictionary(self, pos, allele, type):
        if pos not in self.positional_allele_dictionary:
            self.positional_allele_dictionary[pos] = {}
        if type not in self.positional_allele_dictionary[pos]:
            self.positional_allele_dictionary[pos][type] = {}
        if allele not in self.positional_allele_dictionary[pos][type]:
            self.positional_allele_dictionary[pos][type][allele] = 0

        self.positional_allele_dictionary[pos][type][allele] += 1

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
            allele = read_sequence[i-alignment_position]
            ref = ref_sequence[i-alignment_position]

            if allele != ref:
                self._update_read_allele_dictionary(i, allele, MISMATCH_ALLELE)
                # increase mismatch count
                self.mismatch_count[i] += 1
                # increase the edit count
                self.edit_count[i] += 1
            # this slows things down a lot. Don't add reference allele to the dictionary if we don't use them
            # else:
                # self._update_read_allele_dictionary(i, allele, MATCH_ALLELE)

    def parse_delete(self, alignment_position, length, ref_sequence, read_name):
        """
        Process a cigar operation that is a delete
        :param alignment_position: Alignment position
        :param length: Length of the delete
        :param read_name: Read Id that we need to save
        :return:

        This method updates the candidates dictionary. Mostly by adding read IDs to the specific positions.
        """
        if alignment_position < self.region_start_position or alignment_position > self.region_end_position:
            return

        # actual delete position starts one after the anchor
        start = alignment_position + 1
        stop = start + length

        for i in range(start, stop):
            # increase the coverage
            self.coverage[i] += 1

        # the allele is the anchor + what's being deleted
        allele = self.reference_dictionary[alignment_position] + ref_sequence

        # record the delete where it first starts
        self._update_read_allele_dictionary(alignment_position + 1, allele, DELETE_ALLELE)

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
        self.mismatch_count[alignment_position] += 1

        # the allele is the anchor + what's being deleted
        allele = self.reference_dictionary[alignment_position] + read_sequence
        # print(alignment_position, allele)
        # record the insert where it first starts
        self._update_read_allele_dictionary(alignment_position + 1, allele, INSERT_ALLELE)

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
        self.read_allele_dictionary = {}
        ref_alignment_start = read.reference_start
        ref_alignment_stop = self.get_read_stop_position(read)
        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence.upper()
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
            ref_index_increment, read_index_increment = \
                self.parse_cigar_tuple(cigar_code=cigar_code,
                                       length=length,
                                       alignment_position=ref_alignment_start+ref_index,
                                       ref_sequence=ref_sequence_segment,
                                       read_sequence=read_sequence_segment,
                                       read_name=read.query_name)

            # increase the read index iterator
            read_index += read_index_increment
            ref_index += ref_index_increment

        # after collecting all alleles from reads, update the global dictionary
        for position in self.read_allele_dictionary.keys():
            if position < self.region_start_position or position > self.region_end_position:
                continue
            for record in self.read_allele_dictionary[position]:
                # there can be only one record per position in a read
                allele, allele_type = record
                if allele_type == MATCH_ALLELE or allele_type == MISMATCH_ALLELE:
                    # If next allele is indel then group it with the current one, don't make a separate one
                    if position + 1 <= ref_alignment_stop and position + 1 in self.read_allele_dictionary.keys():
                        next_allele, next_allele_type = list(self.read_allele_dictionary[position + 1].keys())[0]
                        if next_allele_type == INSERT_ALLELE or next_allele_type == DELETE_ALLELE:
                            continue
                    self._update_positional_allele_dictionary(position, allele, allele_type)
                else:
                    # it's an insert or delete, so, add to the previous position
                    self._update_positional_allele_dictionary(position-1, allele, allele_type)

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

        # deal different kinds of operations
        if cigar_code == 0:
            # match
            self.parse_match(alignment_position=alignment_position,
                             length=length,
                             read_sequence=read_sequence,
                             ref_sequence=ref_sequence,
                             read_name=read_name)
        elif cigar_code == 1:
            # insert
            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            self.parse_insert(alignment_position=alignment_position-1,
                              read_sequence=read_sequence,
                              read_name=read_name)
            ref_index_increment = 0
        elif cigar_code == 2 or cigar_code == 3:
            # delete or ref_skip
            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            self.parse_delete(alignment_position=alignment_position-1,
                              ref_sequence=ref_sequence,
                              length=length,
                              read_name=read_name)
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

    def _filter_alleles(self, position, allele_frequency_list):
        filtered_list = list()
        for allele, count in allele_frequency_list:
            frequency = self.coverage[position] * 100 / count
            if count > MIN_MISMATCH_THRESHOLD and frequency > MIN_MISMATCH_PERCENT_THRESHOLD:
                filtered_list.append(allele)
        return filtered_list

    def _get_substitution_record(self, pos, alt1, alt2, ref):
        return [self.chromosome_name, pos, pos, ref, alt1, alt2, 'SUB']

    def _get_insert_record(self, pos, alt1, alt2, ref):
        return [self.chromosome_name, pos, pos, ref, alt1, alt2, 'IN']

    def _get_delete_record(self, pos, alt1, alt2, ref):
        if len(alt1) < len(alt2):
            alt1, alt2 = alt2, alt1
        ref, alt1 = alt1, ref
        pos_end = pos + len(ref) - 1
        return [self.chromosome_name, pos, pos_end, ref, alt1, alt2, 'DEL']

    def parse_reads_and_select_candidates(self, reads):
        """
        Parse reads to aligned to a site to find variants
        :param reads: Set of reads aligned
        :return:
        """
        start_time = time.time()
        total_reads = 0
        for read in reads:
            # check if the mapping quality of the read is above threshold
            if read.mapping_quality > DEFAULT_MIN_MAP_QUALITY:
                self.find_read_candidates(read=read)
                total_reads += 1
        # print('Read processing time', time.time()-start_time)

        positional_selected_alleles = dict()

        for pos in range(self.region_start_position, self.region_end_position):
            if pos not in self.positional_allele_dictionary:
                continue

            ref = self.reference_dictionary[pos]

            # generate a 3 position list, 1 for each edit type
            position_data = [list() for type in [MISMATCH_ALLELE, INSERT_ALLELE, DELETE_ALLELE]]
            position_supported = False

            for type_of_record in self.positional_allele_dictionary[pos]:
                if type_of_record == MATCH_ALLELE:
                    continue
                all_allele_dictionary = self.positional_allele_dictionary[pos][type_of_record]

                # pick the top 2 most frequent allele
                allele_frequency_list = list(sorted(all_allele_dictionary.items(), key=operator.itemgetter(1), reverse=True))[:PLOIDY]

                allele_list = self._filter_alleles(pos, allele_frequency_list)

                # technically this should be expanded to scale with PLOIDY parameter...
                alt1 = allele_list[0] if len(allele_list) >= 1 else None
                alt2 = allele_list[1] if len(allele_list) >= 2 else '.'

                # if there are no edits of this type, skip
                if alt1 is None:
                    continue
                else:
                    # split allele strings and frequencies into their own single typed lists
                    alleles, frequencies = map(list, zip(*allele_frequency_list))
                    coverage = self.coverage[pos]
                    filtering_data = [alleles, frequencies, coverage]

                if type_of_record == MISMATCH_ALLELE:
                    record_data = self._get_substitution_record(pos, alt1, alt2, ref)
                    position_data[type_of_record-1] = record_data + filtering_data
                    position_supported = True

                elif type_of_record == INSERT_ALLELE:
                    record_data = self._get_insert_record(pos, alt1, alt2, ref)
                    position_data[type_of_record-1] = record_data + filtering_data
                    position_supported = True

                elif type_of_record == DELETE_ALLELE:
                    record_data = self._get_delete_record(pos, alt1, alt2, ref)
                    position_data[type_of_record-1] = record_data + filtering_data
                    position_supported = True

            if position_supported:
                positional_selected_alleles[pos] = position_data

        return positional_selected_alleles
