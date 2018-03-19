import numpy

"""
            possible combos:
            gt1       gt2     Candidate validated?
            --------------------------------------
            hom       hom     no
            het       hom     yes
            het       het     yes
            hom_alt   hom     yes
            hom       None    no
            het       None    yes
            hom_alt   None    yes

            impossible combos:
            gt1       gt2     Candidate validated?
            --------------------------------------
            hom_alt   hom_alt NA
            het       hom_alt NA
            hom       hom_alt NA
            hom_alt   het     NA
            hom       het     NA
            None      None    NA

"""

GENOTYPE_DICT = {"Hom": 0, "Het": 1, "Hom_alt": 2}
GENOTYPE_NAMES = ["Hom", "Het", "Hom_alt"]
VCF_OFFSET = 1
PLOIDY = 8

# DEBUG_FREQUENCIES = False
DEBUG_PRINT_ALL = False

# Candidate data indexes
CHR_NAME = 0
START = 1
STOP = 2
REF_INDEX = 3
ALT1 = 4
ALT2 = 5
CANDIDATE_TYPE = 6
ALLELES = 7
FREQUENCIES = 8
COVERAGE = 9

# Positional vcf indexes
SNP, IN, DEL = 0, 1, 2
SNP_DEL = 3

# VCF record indexes
REF, ALT, GT = 0, 1, 2

# Genotype codes
HOM, HET, HOM_ALT = 0, 1, 2

class CandidateLabeler:
    def __init__(self, fasta_handler):
        """
        Initialize candidateLabeler object
        :param fasta_handler: module that fetches reference sequence substrings from a FASTA file
        """
        self.fasta_handler = fasta_handler      # unfortunately need to query reference sequence on a per-site basis

        self.vcf_offset = -1                    # pysam vcf coords are 1-based ... >:[ ... this is what Kishwar wanted
        self.delete_char = '*'

    def _handle_insert(self, rec):
        """
        Process a record that has an insert
        :param rec: VCF record
        :return: attributes of the record
        """
        ref_seq = rec.ref   # no change necessary
        alt_seq = rec.alt   # no change necessary

        pos = rec.pos + self.vcf_offset
        return pos, ref_seq, alt_seq, rec.type

    def _handle_delete(self, rec):
        """
        Process a record that has deletes.
        Deletes are usually grouped together, so we break each of the deletes to make a list.
        :param rec: VCF record containing a delete
        :return: A list of delete attributes
        """
        delete_list = []
        for i in range(0, len(rec.ref)):
            if i < len(rec.alt):
                continue
            ref_seq = rec.ref[i]
            alt_seq = '*'
            pos = rec.pos + i + self.vcf_offset
            genotype = rec.type
            delete_list.append((pos, ref_seq, alt_seq, genotype))
        return delete_list

    @staticmethod
    def get_label_of_allele(positional_vcf, variant_type, candidate_allele):
        """
        Given positional VCFs (IN, DEL, SNP), variant type and a candidate allele, return the try genotype.
        :param positional_vcf: Three dictionaries for each position
        :param variant_type: Type of variant (IN, DEL or SNP)
        :param candidate_allele: Candidate allele
        :return: genotype
        """
        # candidate attributes
        pos_start, pos_stop, ref, alts = candidate_allele
        if variant_type == DEL:
            del_alts = [ref]
            if len(alts) > 1:
                mx_alt = alts[0] if len(alts[0]) > len(alts[1]) else alts[1]
                del_alts = del_alts + [mx_alt]
            alts = del_alts
        # get all records of that position
        records = []
        if pos_start + VCF_OFFSET in positional_vcf.keys():
            records = positional_vcf[pos_start + VCF_OFFSET][variant_type]

        vcf_recs = []
        for record in records:
            # get the alt allele of the record
            rec_alt = record.alt
            if variant_type == DEL:
                rec_alt = record.ref

            if record.type == '':
                record.type = 'Hom'
            # if the alt allele of the record is same as candidate allele
            vcf_recs.append((rec_alt, record.type, record.filter))
        gts = []
        for alt in alts:
            gt = [0, '.']
            for vcf_alt in vcf_recs:
                if alt == vcf_alt[0]:
                    gt = [GENOTYPE_DICT[vcf_alt[1]], vcf_alt[2]]
            gts.extend(gt)
        return gts

    def get_index_by_rec_type(self, rec_type):
        if rec_type == 'SUB':
            return SNP
        if rec_type == 'IN':
            return IN
        if rec_type == 'DEL':
            return DEL

    def _get_genotype_label(self, positional_vcf, start, stop, ref_seq, alleles, rec_type):
        """
        Create a list of dictionaries of 3 types of alleles that can be in a position.

        In each position there can be Insert allele, SNP or Del alleles.
        For total 6 alleles, this method returns 6 genotypes
        :param positional_vcf: VCF records of each position
        :param start: Allele start position
        :param stop: Allele stop position
        :param ref_seq: Reference sequence
        :param alleles_snp: SNP Alleles
        :param alleles_insert: Insert Alleles
        :return: Dictionary of genotypes
        """
        gts = self.get_label_of_allele(positional_vcf, self.get_index_by_rec_type(rec_type), (start, stop, ref_seq, alleles))

        return gts

    @staticmethod
    def _is_supported(genotypes):
        """
        Check if genotype has anything other than Hom
        :param genotypes: Genotype tuple
        :return: Boolean [True if it has Het of Hom_alt]
        """
        supported = False

        if sum(genotypes) > 0:
            supported = True

        return supported

    def _generate_list(self, chromosome_name, start, stop, alleles_snp, alleles_in, ref_seq, genotypes):
        """
        Generate a list of attributes that can be saved of a labeled candidate
        :param chromosome_name: Name of chromosome
        :param start: Allele start position
        :param stop: Allele end position
        :param alleles: All alleles
        :param alleles_insert: Insert alleles
        :param ref_seq: reference Sequence
        :param genotypes: Genotypes
        :return: A list containing (chr start stop ref_seq alt1 alt2 gt1 gt2)
        """
        all_candidates = []
        for i, allele_tuple in enumerate(alleles_snp):
            allele, freq = allele_tuple
            gt = genotypes[SNP][i][0]
            gt_q = genotypes[SNP][i][1]
            gt_f = genotypes[SNP][i][2]
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, gt, gt_q, gt_f])

        for i, allele_tuple in enumerate(alleles_in):
            allele, freq = allele_tuple
            gt = genotypes[IN][i][0]
            gt_q = genotypes[IN][i][1]
            gt_f = genotypes[IN][i][2]
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, gt, gt_q, gt_f])

        return all_candidates

    def _generate_fixed_size_freq_list(self, site_frequencies, vector_length=PLOIDY):
        frequency_list = list()

        # concatenate frequencies for Match, Insert, Delete in order
        for i,frequencies in enumerate(site_frequencies):
            length_difference = vector_length-len(frequencies)
            frequencies.extend([0]*length_difference)

            frequency_list.extend(frequencies)

        return frequency_list

    def _get_chromosome_number(self, chromosome_name):
        number = int(chromosome_name.split("chr")[-1])

        return number

    def _generate_data_vector(self, chromosome_name, start, genotypes, frequencies, alleles, coverage, support):
        data_list = list()

        # convert list of frequencies into vector with length 3*PLOIDY
        site_frequency_list = self._generate_fixed_size_freq_list(frequencies)

        # normalize by coverage depth
        site_frequency_vector = numpy.array(site_frequency_list, dtype=numpy.float32)/coverage

        chromosome_number = self._get_chromosome_number(chromosome_name)
        label = int(support)

        # cap and normalize coverage (max = 1000)
        coverage = min(coverage, 1000)
        coverage = float(coverage)/1000

        data_list.append(chromosome_number) # 0
        data_list.append(start)             # 1
        data_list.extend(genotypes)         # 2-4
        data_list.extend([0]*PLOIDY*3)      # 5-28
        data_list.append(coverage)          # 29
        data_list.append(label)             # 30

        data_vector = numpy.array(data_list, dtype=numpy.float32).reshape((len(data_list),1))

        data_vector[5:29] = site_frequency_vector.reshape((PLOIDY*3,1))

        return data_vector

    def get_labeled_candidates(self, chromosome_name, positional_vcf, candidate_sites):
        """
        Label candidates given variants from a VCF
        :param positional_vcf: IN/DEL/SNPs separated into VCF records, expanded into a 1-to-1 ref_pos:variant allele
        :param candidate_sites: Candidates with the format list of lists, where each sublist is data pertaining to a
        type, MISMATCH, INSERT, or DELETE:
            [['chr1', 10400, 10400, 'T', 'A', 'G', 'SUB', ['A', 'G', 'C'], [11, 11, 7], 500],
            ['chr1', 10400, 10400, 'T', 'TA', '.', 'IN', ['TA'], [4], 500],
            ['chr1', 10400, 10401, 'TA', 'T', '.', 'DEL', ['TA'], [6], 500]]
        :return: List of labeled candidate sites
        """
        # list of all labeled training vectors
        all_labeled_vectors = list()

        # for each candidate
        for pos in candidate_sites:
            site_frequencies = list()
            site_alleles = list()
            site_genotypes = list()

            site_start = None
            site_coverage = None
            site_PASS = True

            for candidate in candidate_sites[pos]:
                if len(candidate) > 0:
                    allele_start = candidate[START]
                    allele_stop = candidate[STOP]
                    ref_sequence = candidate[REF_INDEX]
                    alt1 = candidate[ALT1]
                    alt2 = candidate[ALT2]
                    rec_type = candidate[CANDIDATE_TYPE]
                    all_allele_sequences = candidate[ALLELES]
                    frequencies = candidate[FREQUENCIES]
                    coverage = candidate[COVERAGE]

                    alleles = [alt1, alt2]

                    # test the alleles across IN, DEL, and SNP variant dictionaries
                    genotype = self._get_genotype_label(positional_vcf=positional_vcf,
                                                        start=allele_start,
                                                        stop=allele_stop,
                                                        ref_seq=ref_sequence,
                                                        alleles=alleles,
                                                        rec_type=rec_type)

                    vcf_quality_field = genotype[1]

                    if site_start is None:
                        site_start = allele_start

                    if site_coverage is None:
                        site_coverage = coverage

                    if vcf_quality_field != "PASS" and vcf_quality_field != '.':
                        site_PASS = False

                    site_frequencies.append(frequencies)
                    site_alleles.append(all_allele_sequences)
                    site_genotypes.append(genotype[0])


                else:
                    # generate null vectors, since the site contains no edits of this type
                    frequencies = [0]*PLOIDY
                    all_allele_sequences = ['']*PLOIDY

                    site_frequencies.append(frequencies)
                    site_alleles.append(all_allele_sequences)
                    site_genotypes.append(0)

            if site_PASS:
                support = self._is_supported(site_genotypes)

                vector = self._generate_data_vector(chromosome_name=chromosome_name,
                                                    start=site_start,
                                                    genotypes=site_genotypes,
                                                    frequencies=site_frequencies,
                                                    alleles=site_alleles,
                                                    coverage=site_coverage,
                                                    support=support)

                all_labeled_vectors.append(vector)

        # if there is no data for this region, append an empty vector
        if len(all_labeled_vectors) == 0:
            all_labeled_vectors = numpy.array([[]])
        else:
            numpy.concatenate(all_labeled_vectors, axis=1)

        return all_labeled_vectors
