"""
Controller
"""
from collections import defaultdict

class CandidateFinder:
    def __init__(self, reads, fasta_handler, chromosome_name, window_ref_start_position):
        self.window_ref_start_position = window_ref_start_position
        self.chromosome_name = chromosome_name
        self.fasta_handler = fasta_handler
        self.reads = reads

        self.cigar_key = {0: "MATCH",
                          1: "INSERT",
                          2: "DELETE",
                          3: "REFSKIP",
                          4: "SOFTCLIP",
                          5: "HARDCLIP",
                          6: "PAD"}

        self.ref_advancement_by_cigar = {"MATCH": True,
                                         "INSERT": False,
                                         "DELETE": True,
                                         "REFSKIP": False,
                                         "SOFTCLIP": False,
                                         "HARDCLIP": False,
                                         "PAD": False}

        self.read_advancement_by_cigar = {"MATCH": True,
                                          "INSERT": True,
                                          "DELETE": False,
                                          "REFSKIP": True,
                                          "SOFTCLIP": True,
                                          "HARDCLIP": False,
                                          "PAD": True}

        self.candidates = defaultdict(list)

    def get_read_stop_position(self,read):
        ref_alignment_stop = read.reference_end
        if ref_alignment_stop == None:
            positions = read.get_reference_positions()

            print(read.mapping_quality)
            print(read.query_alignment_qualities)
            # find last entry that isn't None
            i = len(positions) - 1
            ref_alignment_stop = positions[-1]
            while i > 0 and ref_alignment_stop is None:
                i -= 1
                ref_alignment_stop = positions[i]

        return ref_alignment_stop

    def parse_reads(self, reads):
        for read in reads:
            if read.mapping_quality > 0:
                print(read.query_name)
                self.find_read_candidates(read=read)

    def parse_match(self, alignment_position, read_sequence, ref_sequence, read_name):
        for i in range(len(read_sequence)):
            if read_sequence[i] != ref_sequence[i]:
                self.candidates[alignment_position].append(read_name)

    def parse_delete(self, alignment_position, length, read_name):
        start = alignment_position-1
        stop = alignment_position + length + 1
        for i in range(start,stop):
            self.candidates[i].append(read_name)

    def parse_insert(self, alignment_position, length, read_name):
        self.candidates[alignment_position].append(read_name)

    def find_read_candidates(self, read):
        read_candidates = set()
        ref_alignment_start = read.reference_start
        ref_alignment_stop = self.get_read_stop_position(read)
        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence
        ref_sequence = self.fasta_handler.get_sequence(chromosome_name=self.chromosome_name,
                                                       start=ref_alignment_start,
                                                       stop=ref_alignment_stop)

        corresponding_ref_sequence = ''
        expanded_cigar = ''

        # read_index: index of read sequence
        # alignment_position: position that the current cigar tuple aligns to in reference
        read_index = 0
        ref_index = 0
        for cigar in cigar_tuples:
            cigar_code = cigar[0]
            length = cigar[1]

            ref_sequence_segment = ref_sequence[ref_index:ref_index+length]
            read_sequence_segment = read_sequence[read_index:read_index+length]

            ref_index_increment, read_index_increment, expanded_cigar_segment = \
                self.parse_cigar_tuple(cigar_code=cigar_code,
                                       length=length,
                                       alignment_position=ref_alignment_start+ref_index,
                                       ref_sequence=ref_sequence_segment,
                                       read_sequence=read_sequence_segment,
                                       read_name=read.query_name)

            read_index += read_index_increment
            ref_index += ref_index_increment

            expanded_cigar += expanded_cigar_segment
            corresponding_ref_sequence += ref_sequence_segment

            # print("----------")
            # print(ref_sequence_segment)
            # print(expanded_cigar_segment)
            # print(read_sequence_segment)

        # print("--------------------------")
        print(ref_alignment_start)
        print(corresponding_ref_sequence)
        print(expanded_cigar)
        print(read_sequence)
        print(read_candidates)
        print("--------------------------")


    def parse_cigar_tuple(self, cigar_code, length, alignment_position, ref_sequence, read_sequence, read_name):
        cigar_operation = self.cigar_key[cigar_code]
        ref_index_increment = 0
        read_index_increment = 0
        expanded_cigar_string = ''
        candidates = set()

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


