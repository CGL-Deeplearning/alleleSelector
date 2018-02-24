import numpy
import sys
sys.path.insert(0, '..')
from modules.IntervalTree import IntervalTree
import time

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
PLOIDY = 10

# DEBUG_FREQUENCIES = False
DEBUG_PRINT_ALL = False

# Candidate data indexes
START = 0
REF_INDEX = 1
IN_ALLELES = 2
SNP_ALLELES = 3
COVERAGE = 4

# Positional vcf indexes
SNP, IN, DEL = 0, 1, 2
SNP_DEL = 3

# VCF record indexes
REF,ALT,GT = 0,1,2

# Genotype codes
HOM,HET,HOM_ALT = 0,1,2

# training set
DATA, LABEL = 0,1

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
        pos_start, pos_stop, ref, alt_tuple = candidate_allele
        alt, freq = alt_tuple
        # by default the genotype is Hom
        gt = HOM

        for i in range(pos_start, pos_stop+1):
            if i in positional_vcf.keys():
                # get all records of that position
                records = positional_vcf[i][variant_type]
                for record in records:
                    # get the alt allele of the record
                    rec_alt = record[ALT]

                    # print(candidate_allele, rec_alt, record[GT])

                    # if the alt allele of the record is same as candidate allele
                    if rec_alt == alt:
                        # then  we get the genotype
                        gt = record[GT]

        return gt

    def _get_all_genotype_labels(self, positional_vcf, start, stop, ref_seq, alleles_snp, alleles_insert):
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
        gts_in = list()
        gts_snp = list()

        for alt in alleles_insert:
            gt_in = self.get_label_of_allele(positional_vcf, IN, (start, stop, ref_seq, alt))
            gts_in.append(gt_in)

        for alt in alleles_snp:
            gt_snp = self.get_label_of_allele(positional_vcf, SNP, (start, stop, ref_seq, alt))
            gts_snp.append(gt_snp)

        return [gts_snp, gts_in]

    @staticmethod
    def _is_supported(genotypes):
        """
        Check if genotype has anything other than Hom
        :param genotypes: Genotype tuple
        :return: Boolean [True if it has Het of Hom_alt]
        """
        supported = False
        gt_set = set(genotypes)

        if len(gt_set) == 0:
            supported = False
        elif len(gt_set) == 1:
            if HOM in gt_set:
                supported = False
            else:
                supported = True
        elif len(gt_set) > 1:
            supported = True

        return supported

    def _is_position_supported(self, genotypes):
        """
        Check if a position has any genotype other than Hom
        :param genotypes: Genotypes list of that position
        :return: Boolean [True if it has Het or Hom_alt]
        """
        supported_types = list()
        for sublist in genotypes:
            supported_types.append(self._is_supported(sublist))

        return any(supported_types)

    def _generate_list(self, chromosome_name, start, stop, alleles_snp, alleles_in, ref_seq, genotypes):
        """
        Generate a list of attributes that can be saved of a labeled candidate
        :param chromosome_name: Name of chromosome
        :param start: Allele start position
        :param stop: Allele end position
        :param alleles_snp: All non-insert alleles
        :param alleles_in: Insert alleles
        :param ref_seq: reference Sequence
        :param genotypes: Genotypes
        :return: A list containing (chr start stop ref_seq alt1 alt2 gt1 gt2)
        """
        all_candidates = []
        for i, allele_tuple in enumerate(alleles_snp):
            allele, freq = allele_tuple
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, genotypes[SNP][i]])

        for i, allele_tuple in enumerate(alleles_in):
            allele, freq = allele_tuple
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, genotypes[IN][i]])

        return all_candidates

    @staticmethod
    def _convert_to_fixed_length_vectors(allele_frequencies, vector_length=PLOIDY):
        if len(allele_frequencies) == 0:
            return [0]*vector_length

        else:
            alleles, frequencies = zip(*allele_frequencies)

            frequencies = list(frequencies)
            frequencies.extend([0]*(vector_length-len(frequencies)))

            return frequencies

    def _get_chromosome_number(self, chromosome_name):
        number = int(chromosome_name.split("chr")[-1])

        return number

    def _get_separate_genotypes(self, genotypes):
        gt_snp, gt_in = genotypes

        if len(gt_snp) == 0:
            gt_snp = 0
        else:
            gt_snp = gt_snp[0]

        if len(gt_in) == 0:
            gt_in = 0
        else:
            gt_in = gt_in[0]

        return gt_snp, gt_in

    def _build_training_vector(self, chromosome_number, start, gt_insert, gt_snp, alleles_insert, alleles_snp, support, coverage_depth):
        freq_insert = self._convert_to_fixed_length_vectors(alleles_insert)
        freq_snp = self._convert_to_fixed_length_vectors(alleles_snp)

        # combined length of frequency vectors
        freq_length = len(freq_insert)+len(freq_snp)

        # initialize data vector for training with extra slots for:
        #   1 coverage depth
        #   2 training label
        #   3 het/hom/hom_alt for the SNPs/Deletes
        #   4 het/hom/hom_alt for the Inserts
        #   5 chromosome number
        #   6 start position in chromosome

        # create composite vector of snps and inserts and normalize by coverage depth
        freq_vector = numpy.array(freq_snp + freq_insert)
        freq_vector = freq_vector/coverage_depth

        # print(chromosome_number, start, gt_snp, gt_insert)

        # for troubleshooting
        metadata_vector = [chromosome_number, start, gt_snp, gt_insert]

        # the training label
        label = int(support)

        # allocate numpy vector and assign data to it
        vector = numpy.zeros((freq_length+6, 1))
        vector[0:4,0] = metadata_vector
        vector[4:4+freq_length,0] = freq_vector
        vector[-2,0] = float(coverage_depth)/1000
        vector[-1,0] = label

        # print(vector)

        return vector

    def get_labeled_candidates(self, chromosome_name, positional_vcf, candidate_sites, confident_intervals):
        """
        Label candidates given variants from a VCF
        :param positional_vcf: IN/DEL/SNPs separated into VCF records, expanded into a 1-to-1 ref_pos:variant allele
        :param candidate_sites: Candidates
        :return: List of labeled candidate sites
        """

        # list of all labeled candidates
        all_labeled_frequencies = list()

        interval_tree = IntervalTree(confident_intervals)

        # for each candidate
        for c,candidate_site in enumerate(sorted(candidate_sites)):
            allele_start = candidate_site[START]
            allele_stop = candidate_site[START]

            interval = [allele_start, allele_start]

            # print(allele_start, candidate_site[REF_INDEX], candidate_site[SNP_ALLELES], candidate_site[IN_ALLELES])

            if interval in interval_tree:
                ref_sequence = candidate_site[REF_INDEX]
                alleles_insert = candidate_site[IN_ALLELES]
                alleles_snp = candidate_site[SNP_ALLELES]
                coverage_depth = candidate_site[COVERAGE]

                chromosome_number = self._get_chromosome_number(chromosome_name)

                # test the alleles across IN, DEL, and SNP variant dictionaries
                genotypes = self._get_all_genotype_labels(positional_vcf=positional_vcf,
                                                          start=allele_start,
                                                          stop=allele_stop,
                                                          ref_seq=ref_sequence,
                                                          alleles_snp=alleles_snp,
                                                          alleles_insert=alleles_insert)

                gt_snp, gt_insert = self._get_separate_genotypes(genotypes)

                support = self._is_position_supported(genotypes)

                vector = self._build_training_vector(chromosome_number=chromosome_number,
                                                     start=allele_start,
                                                     gt_snp=gt_snp,
                                                     gt_insert=gt_insert,
                                                     alleles_insert=alleles_insert,
                                                     alleles_snp=alleles_snp,
                                                     support=support,
                                                     coverage_depth=coverage_depth)

                all_labeled_frequencies.append(vector)

                # print(allele_start,allele_stop,genotypes,support,coverage_depth)
                # print(numpy.array2string(vector, separator="\t", precision=2, max_line_width=400))

                if DEBUG_PRINT_ALL:
                    print()
                    print(allele_start, allele_stop, ref_sequence)
                    print("snp alleles:        ", alleles_snp)
                    print("insert alleles: ", alleles_insert)
                    print("insert gts:     ", genotypes[IN])
                    print("snp gts:        ", genotypes[SNP])

        # if there is no data for this region, append an empty vector
        if len(all_labeled_frequencies) == 0:
            all_labeled_frequencies.append(numpy.array(list()))

        try:
            numpy.concatenate(all_labeled_frequencies, axis=1)
        except:
            print("WARNING: Failed Concatenation: ", chromosome_name, positional_vcf.keys()[0], positional_vcf.keys()[-1])
            print("length: ",len(all_labeled_frequencies))
            if len(all_labeled_frequencies) > 0:
                print(all_labeled_frequencies[0])

        return
