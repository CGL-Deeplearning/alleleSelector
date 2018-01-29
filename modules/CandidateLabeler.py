from collections import Counter,defaultdict

'''
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

'''

# DEBUG_FREQUENCIES = False
DEBUG_PRINT_ALL = False

START = 0
STOP = 1
IN_ALLELES = 2
ALLELES = 3

SNP,IN,DEL = 0,1,2
SNP_DEL = 3

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

    def generate_position_based_vcf(self, variants):
        """
        Split a VCF to 3 different record types. SNP, IN or Delete.
        :param variants: Variant records from a VCF file
        :return:
        """
        position_based_vcf = dict()

        for variant_pos in variants:
            for variant in variants[variant_pos]:
                if variant.genotype_class == 'SNP':
                    pos, ref_seq, alt_seq, genotype = variant.pos + self.vcf_offset, variant.ref, variant.alt, variant.type

                    if pos not in position_based_vcf:
                        position_based_vcf[pos] = [[],[],[]]

                        position_based_vcf[pos][SNP].append((pos, ref_seq, alt_seq, genotype, variant.genotype_class))

                elif variant.genotype_class == 'DEL':
                    list_of_records = self._handle_delete(variant)
                    for record in list_of_records:
                        pos, ref_seq, alt_seq, genotype = record

                        if pos not in position_based_vcf:
                            position_based_vcf[pos] = [[], [], []]

                            position_based_vcf[pos][DEL].append((pos, ref_seq, alt_seq, genotype, variant.genotype_class))

                elif variant.genotype_class == 'IN':
                    pos, ref_seq, alt_seq, genotype = self._handle_insert(variant)

                    if pos not in position_based_vcf:
                        position_based_vcf[pos] = [[],[],[]]

                        position_based_vcf[pos][IN].append((pos, ref_seq, alt_seq, genotype, variant.genotype_class))

        return position_based_vcf

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
        pos_start, pos_stop, ref, alt = candidate_allele
        # by default the genotype is Hom
        gt = 'Hom'
        for i in range(pos_start, pos_stop+1):
            if i in positional_vcf.keys():
                # get all records of that position
                records = positional_vcf[i][variant_type]
                for record in records:
                    # get the alt allele of the record
                    rec_alt = record[2]
                    # if the alt allele of the record is same as candidate allele
                    if rec_alt == alt:
                        # then  we get the genotype
                        gt = record[3]

        return gt

    def _get_all_genotype_labels(self, positional_vcf, start, stop, ref_seq, alleles, alleles_insert):
        """
        Create a list of dictionaries of 3 types of alleles that can be in a position.

        In each position there can be Insert allele, SNP or Del alleles.
        For total 6 alleles, this method returns 6 genotypes
        :param positional_vcf: VCF records of each position
        :param start: Allele start position
        :param stop: Allele stop position
        :param ref_seq: Reference sequence
        :param alleles: Alleles
        :param alleles_insert: Insert Alleles
        :return: Dictionary of genotypes
        """
        gts_in = list()
        gts_del = list()
        gts_snp = list()

        for alt in alleles_insert:
            gt_in = self.get_label_of_allele(positional_vcf, IN, (start, stop, ref_seq, alt))

            gts_in.append(gt_in)

        for alt in alleles:
            gt_del = self.get_label_of_allele(positional_vcf, DEL, (start, stop, ref_seq, alt))
            gt_snp = self.get_label_of_allele(positional_vcf, SNP, (start, stop, ref_seq, alt))

            gts_snp.append(gt_snp)
            gts_del.append(gt_del)

        return [gts_snp, gts_in, gts_del]

    @staticmethod
    def _is_supported(genotypes):
        """
        Check if genotype has anything other than Hom
        :param genotypes: Genotype tuple
        :return: Boolean [True if it has Het of Hom_alt]
        """
        supported = False
        gt_set = set(genotypes)

        if len(gt_set) == 1:
            if "Hom" in gt_set:
                supported = False
            else:
                supported = True

        if len(gt_set) == 0:
            supported = False

        return supported

    def _is_position_supported(self, genotypes):
        """
        Check if a position has any genotype other than Hom
        :param genotypes: Genotypes list of that position
        :return: Boolean [True if it has Het or Hom_alt]
        """
        in_supported = self._is_supported(genotypes[IN])
        del_supported = self._is_supported(genotypes[IN])
        snp_supported = self._is_supported(genotypes[IN])

        return in_supported or del_supported or snp_supported

    def _generate_list(self, chromosome_name, start, stop, alleles, alleles_insert, ref_seq, genotypes):
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

        is_insert = None
        merged_genotype = ['Hom', 'Hom']
        supported = set()

        if self._is_supported(genotypes[SNP]) and self._is_supported(genotypes[DEL]):
            variant_type = "SNP_DEL"
            is_insert = False
            supported.add(variant_type)
            alleles = [allele if allele is not None else "None" for allele in alleles]
            merged_genotype = self._merge_snp_and_del(alleles, genotypes[SNP], genotypes[DEL])

        elif self._is_supported(genotypes[SNP]):
            variant_type = "SNP"
            is_insert = False
            supported.add(variant_type)
            alleles = [allele if allele is not None else "None" for allele in alleles]
            merged_genotype = genotypes[SNP]

        elif self._is_supported(genotypes[IN]):
            variant_type = "IN"
            is_insert = True
            supported.add(variant_type)
            alleles = [allele if allele is not None else "None" for allele in alleles_insert]
            merged_genotype = genotypes[IN]

        elif self._is_supported(genotypes[DEL]):
            variant_type = "DEL"
            is_insert = False
            supported.add(variant_type)
            alleles = [allele if allele is not None else "None" for allele in alleles]
            merged_genotype = genotypes[DEL]

        if len(supported) == 0:     # no true variant
            variant_type = "SNP"
            is_insert = False
            alleles = [None]
            merged_genotype = ["Hom"]

        data_list = [chromosome_name, start, stop, is_insert, ref_seq] + alleles + merged_genotype
        # bed conversion can't handle None type, so convert all None to string 'None'
        data_list = ['None' if v is None else v for v in data_list]

        print(data_list)

        return data_list

    @staticmethod
    def _merge_snp_and_del(alleles, snp_genotypes, del_genotypes):
        """
        Merge the SNP and DEL genotypes as they are interchangeable
        :param alleles: Alleles
        :param snp_genotypes: SNP genotypes
        :param del_genotypes: Del genotypes
        :return: Merged genotypes
        """
        output_genotypes = ["Hom", "Hom"]
        del_genotype_index = alleles.index('*')
        snp_genotype_index = 1 - del_genotype_index

        output_genotypes[del_genotype_index] = del_genotypes[del_genotype_index]
        output_genotypes[snp_genotype_index] = snp_genotypes[snp_genotype_index]

        return output_genotypes

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

        # keep only the top 2 most frequent alleles (assuming diploidy)
        i = 0
        top_2_inserts = list()
        top_2_non_inserts = list()
        alleles = [entry[0] for entry in allele_counter.most_common()]
        while (len(top_2_inserts) < 2 or len(top_2_non_inserts) < 2) and i < len(alleles):
            if len(alleles[i]) > len(ref_sequence):         # insert
                if len(top_2_inserts) < 2:
                    top_2_inserts.append(alleles[i])
            else:                                           # SNP or delete
                if len(top_2_non_inserts) < 2:
                    top_2_non_inserts.append(alleles[i])
            i += 1

        top_2_inserts += [None]*(2-len(top_2_inserts))
        top_2_non_inserts += [None]*(2-len(top_2_non_inserts))

        return top_2_inserts, top_2_non_inserts

    @staticmethod
    def _get_allele_frequency_vector(allele_list, ref_sequence, vector_length=8, normalize_by_depth=True):
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
        # frequencies = numpy.array(frequencies)/len(alleles_sequences)

        return frequencies

    def get_labeled_candidates(self, chromosome_name, variants, candidate_sites):
        """
        Label candidates given variants from a VCF
        :param variants: VCF records
        :param candidate_sites: Candidates
        :return: List of labeled candidate sites
        """
        # get separate positional variant dictionaries for IN, DEL, and SNP
        positional_vcf = self.generate_position_based_vcf(variants)


        # create a log of whether variant positions have been matched to a candidate
        validated_vcf_positions = {key: False for key in positional_vcf.keys()}

        # list of all labeled candidates
        all_labeled_candidates = []

        # for each candidate
        for candidate_site in candidate_sites:
            allele_start = candidate_site[START]
            allele_stop = candidate_site[STOP]
            alleles_insert = candidate_site[IN_ALLELES]
            alleles = candidate_site[ALLELES]

            ref_sequence = self.fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                           start=allele_start,
                                                           stop=allele_stop + 1)

            # test the alleles across IN, DEL, and SNP variant dictionaries
            genotypes = self._get_all_genotype_labels(positional_vcf=positional_vcf,
                                                      start=allele_start,
                                                      stop=allele_stop,
                                                      ref_seq=ref_sequence,
                                                      alleles=alleles,
                                                      alleles_insert=alleles_insert)

            # get a list of attributes that can be saved
            labeled_candidate_data = self._generate_list(chromosome_name, allele_start, allele_stop,
                                                         alleles, alleles_insert, ref_sequence, genotypes)

            all_labeled_candidates.append(labeled_candidate_data)

            if DEBUG_PRINT_ALL:
                print()
                print(allele_start,allele_stop,ref_sequence)
                print("alleles:        ", alleles)
                print("insert alleles: ", alleles_insert)
                print("insert gts:     ", genotypes[IN])
                print("del gts:        ", genotypes[DEL])
                print("snp gts:        ", genotypes[SNP])

            # test if the site has any true variants (not Hom)
            if self._is_position_supported(genotypes):

                if allele_start in validated_vcf_positions:
                    validated_vcf_positions[allele_start] = True

        # # see if any sites weren't supported by candidates
        # for pos in validated_vcf_positions:
        #     if not validated_vcf_positions[pos]:
        #         print("\nWARNING: Unsupported VCF position: ", positional_vcf[pos])
        #         print("\tRecord: ", positional_vcf[pos])

        return all_labeled_candidates