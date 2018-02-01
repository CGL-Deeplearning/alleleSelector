from collections import Counter

"""
Finds allele in a given window. The window is a candidate window where a variant is present.

- CandidateInformation:
    - Information of a candidate allele we need to save.
    - When we find a candidate allele we save:
        - read_id, allele_sequence, map_quality, base_quality, read_direction

- AlleleCandidateList:
    - Stores all candidate alleles present in that window

- AlleleFinder
    - Goes through all the reads and creates alleles that are possible candidates
"""

PLOIDY = 2


class CandidateInformation:
    """
    Creates a structure of candidate allele information

    For each candidate allele we save these information:
    - read_id: which read is proposing this candidate
    - allele_sequence: sequence of the allele
    - map_quality: map quality of the read
    - base_qualities: base qualities of each base in the allele sequence
    - read_direction: if the read is reverse this is true
    """
    def __init__(self, allele_sequence, map_quality, base_qualities, read_id, read_direction):
        """
        Initiate a candidate object
        :param allele_sequence: sequence of the allele
        :param map_quality: map quality of the read
        :param base_qualities: base qualities of each base in the allele sequence
        :param read_id: which read is proposing this candidate
        :param read_direction: if the read is reverse this is true
        """
        self.read_id = read_id
        self.allele_sequence = allele_sequence
        self.map_quality = map_quality
        self.base_qualities = base_qualities
        self.read_direction = read_direction

    def __str__(self):
        """
        Print the candidate
        :return: A string that contains the values of the attribute
        """
        print_str = str(self.read_id) + "\n" + str(self.allele_sequence) + " " \
                    + str(self.map_quality) + " " + str(self.base_qualities)
        return print_str


class AlleleFinder:
    """
    Find alleles in a given window.
    """
    def __init__(self, chromosome_name, window_start, window_end, column_pileups, reference_sequence):
        """
        Initialize AlleleFinder object
        :param chromosome_name: Name of chromosome
        :param window_start: Start position of the window
        :param window_end: End position of the window
        :param column_pileups: All column pileups that belong to the window
        :param reference_sequence: Reference sequence of the window
        """
        self.chromosome_name = chromosome_name
        self.window_start = window_start
        self.window_end = window_end
        self.column_pileups = column_pileups

        # [genomic_position] = [max_insert_length]
        self.insert_length_dictionary = {}
        # [read_id] = {{genomic_position}->base}
        self.read_dictionary = {}
        # [read_id] = {{genomic_position}->insert_bases}
        self.read_insert_dictionary = {}  # used
        # List of Read ids in a genomic position
        self.reads_aligned_to_pos = {}
        self.reference_sequence = reference_sequence

    @staticmethod
    def get_attributes_to_save_insert(pileupcolumn, pileupread):
        """
        If position has insert then return attributes to be saved
        :param pileupcolumn: Pileupcolumn of that position
        :param pileupread: Read that has the insert
        :return: set of attributes to save
        """
        insert_start = pileupread.query_position + 1
        insert_end = insert_start + pileupread.indel

        return pileupcolumn.pos, \
               pileupread.alignment.query_name, \
               pileupread.alignment.query_sequence[insert_start:insert_end], \
               pileupread.alignment.query_qualities[insert_start:insert_end], \
               pileupread.alignment.mapping_quality, \
               pileupread.alignment.is_reverse

    @staticmethod
    def get_attributes_to_save(pileupcolumn, pileupread):
        """
        If position has insert then return attributes to be saved
        :param pileupcolumn: Pileupcolumn of that position
        :param pileupread: Read that has the insert
        :return: set of attributes to save
        """
        if pileupread.is_del:
            return pileupcolumn.pos, \
                   pileupread.alignment.query_name, \
                   '*', \
                   0, \
                   pileupread.alignment.mapping_quality, \
                   pileupread.alignment.is_reverse
        else:
            return pileupcolumn.pos, \
                   pileupread.alignment.query_name, \
                   pileupread.alignment.query_sequence[pileupread.query_position], \
                   pileupread.alignment.query_qualities[pileupread.query_position], \
                   pileupread.alignment.mapping_quality, \
                   pileupread.alignment.is_reverse

    def initialize_dictionaries(self, genomic_position, read_id, is_insert):
        """
        Initialize all the dictionaries for a specific position
        :param genomic_position: Genomic position of interest
        :param read_id: Read id for which dictionaries should be initialized
        :param is_insert: If the position is an insert
        :return:
        """
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

    def save_info_of_a_position(self, genomic_position, read_id, base, base_qual, map_qual, is_rev, is_in):
        """
        Given the attributes of a base at a position
        :param genomic_position: Genomic position
        :param read_id: Read id of a read
        :param base: Base at the position
        :param base_qual: Base quality
        :param map_qual: Map quality
        :param is_rev: If read is reversed
        :param is_in: If position has insert
        :return:
        """
        self.initialize_dictionaries(genomic_position, read_id, is_in)

        if is_in is False:
            self.read_dictionary[read_id][genomic_position] = (base, base_qual, map_qual, is_rev)
        else:
            self.read_insert_dictionary[read_id][genomic_position] = (base, base_qual, map_qual, is_rev)
            self.insert_length_dictionary[genomic_position] = max(self.insert_length_dictionary[genomic_position],
                                                                  len(base))

    def generate_base_dictionaries(self):
        """
        Go through all the positions and update base and insert dictionary
        :return:
        """
        # in each column of pileup columns
        for pileupcolumn in self.column_pileups:
            # initialize the read aligned to position list
            self.reads_aligned_to_pos[pileupcolumn.pos] = []
            # iterate through each read in the pileup
            for pileupread in pileupcolumn.pileups:
                # add read to position
                self.reads_aligned_to_pos[pileupcolumn.pos].append(pileupread.alignment.query_name)
                # if there is an insert
                if pileupread.indel > 0:
                    gen_pos, read_id, base, base_qual, map_qual, is_rev = \
                        self.get_attributes_to_save_insert(pileupcolumn, pileupread)
                    self.save_info_of_a_position(gen_pos, read_id, base, base_qual, map_qual, is_rev, is_in=True)

                # or just take care of the base
                gen_pos, read_id, base, base_qual, map_qual, is_rev = \
                    self.get_attributes_to_save(pileupcolumn, pileupread)
                self.save_info_of_a_position(gen_pos, read_id, base, base_qual, map_qual, is_rev, is_in=False)

    @staticmethod
    def _select_alleles(allele_list, ref_sequence):
        """
        Given a dictionary of allele objects, naively find the top 2 represented sequences and return a dictionary
        that assigns read_ids to each of the 2 alleles.
        :param allele_list: the list of alleles found at a given site
        :param ref_sequence: reference sequence for the candidate window
        :return: top_2_alleles: the most frequent two alleles
        """

        alleles_sequences = [allele.allele_sequence for allele in allele_list]

        # find allele distribution by counting instances of each
        allele_counter = Counter(alleles_sequences)

        if ref_sequence in allele_counter:
            del allele_counter[ref_sequence]    # don't consider the reference allele

        # keep only the top n most frequent alleles (assuming diploidy)
        n = PLOIDY
        i = 0

        top_n_inserts = list()
        # SNPs = substitution + deletion
        top_n_snps = list()
        alleles = [entry[0] for entry in allele_counter.most_common()]
        while (len(top_n_inserts) < n or len(top_n_snps) < n) and i < len(alleles):
            if len(alleles[i]) > len(ref_sequence):         # insert
                if len(top_n_inserts) < n:
                    top_n_inserts.append(alleles[i])
            else:                                           # SNP or delete
                if len(top_n_snps) < n:
                    top_n_snps.append(alleles[i])
            i += 1

        return top_n_inserts, top_n_snps

    @staticmethod
    def _get_allele_frequency_vector(allele_list, ref_sequence, vector_length=8):    # , normalize_by_depth=True):
        """
        Find the sorted frequency vector: the # of alleles observed at a site, sorted in descending order
        :param allele_list: the list of alleles found at a given site
        :return: top_2_alleles: the most frequent two alleles
        """
        # print(allele_list)
        alleles_sequences = [allele.allele_sequence for allele in allele_list]

        # find allele distribution by counting instances of each
        allele_counter = Counter(alleles_sequences)

        try:
            ref_frequency = allele_counter[ref_sequence]
            del allele_counter[ref_sequence]
        except KeyError:
            ref_frequency = 0

        frequencies = sorted(dict(allele_counter).values(), reverse=True)[:vector_length-1]
        frequencies = [ref_frequency] + frequencies
        frequencies += [0] * (vector_length - len(frequencies))

        return frequencies

    def generate_candidate_allele_list(self):
        """
        Generate a list of candidates for the window
        :return: A list of candidate alleles for that window
        """
        reads = self.reads_aligned_to_pos[self.window_start]
        candidate_list = []
        for read in reads:

            candidate_allele = ''
            candidate_base_quality = []
            candidate_map_quality = 60
            candidate_read_is_rev = False
            read_covers_whole_window = True
            for pos in range(self.window_start, self.window_end+1):
                if read in self.read_dictionary and pos in self.read_dictionary[read]:
                    base, base_quality, map_quality, orientation = self.read_dictionary[read][pos]
                    candidate_allele += base
                    candidate_base_quality.append(base_quality)
                    candidate_map_quality = map_quality
                    candidate_read_is_rev = orientation
                else:
                    read_covers_whole_window = False
                    break

                if read in self.read_insert_dictionary and pos in self.read_insert_dictionary[read]:
                    base, base_quality, map_quality, orientation = self.read_insert_dictionary[read][pos]
                    candidate_allele += base
                    candidate_base_quality.extend(base_quality)
                    candidate_map_quality = map_quality
                    candidate_read_is_rev = orientation

            if read_covers_whole_window:
                candidate_info = CandidateInformation(candidate_allele, candidate_map_quality, candidate_base_quality,
                                                        read, candidate_read_is_rev)
                candidate_list.append(candidate_info)

        insert_allele, snp_allele = self._select_alleles(candidate_list, self.reference_sequence)

        return insert_allele, snp_allele

