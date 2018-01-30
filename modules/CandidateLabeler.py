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
PLOIDY = 2

# DEBUG_FREQUENCIES = False
DEBUG_PRINT_ALL = False

# Candidate data indexes
START = 0
STOP = 1
IN_ALLELES = 2
SNP_ALLELES = 3
DEL_ALLELES = 4

# Positional vcf indexes
SNP,IN,DEL = 0,1,2
SNP_DEL = 3

# VCF record indexes
REF,ALT,GT = 0,1,2

# Genotype codes
HOM,HET,HOM_ALT = 0,1,2

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
        pos_start, pos_stop, ref, alt = candidate_allele
        # by default the genotype is Hom
        gt = HOM
        for i in range(pos_start, pos_stop+1):
            if i in positional_vcf.keys():
                # get all records of that position
                records = positional_vcf[i][variant_type]
                for record in records:
                    # get the alt allele of the record
                    rec_alt = record[ALT]
                    # if the alt allele of the record is same as candidate allele
                    if rec_alt == alt:
                        # then  we get the genotype
                        gt = record[GT]

        return gt

    def _get_all_genotype_labels(self, positional_vcf, start, stop, ref_seq, alleles_snp, alleles_insert, alleles_del):
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

        for alt in alleles_snp:
            gt_snp = self.get_label_of_allele(positional_vcf, SNP, (start, stop, ref_seq, alt))
            gts_snp.append(gt_snp)

        for alt in alleles_del:
            gt_del = self.get_label_of_allele(positional_vcf, DEL, (start, stop, ref_seq, alt))
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
        in_supported = self._is_supported(genotypes[IN])
        del_supported = self._is_supported(genotypes[DEL])
        snp_supported = self._is_supported(genotypes[SNP])

        return in_supported or del_supported or snp_supported

    def _generate_list(self, chromosome_name, start, stop, alleles_snp, alleles_in, alleles_del, ref_seq, genotypes):
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
        all_candidates= []
        for i, allele in enumerate(alleles_snp):
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, genotypes[SNP][i]])
        for i, allele in enumerate(alleles_in):
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, genotypes[IN][i]])
        for i, allele in enumerate(alleles_del):
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, genotypes[DEL][i]])

        return all_candidates

    def get_labeled_candidates(self, chromosome_name, positional_vcf, candidate_sites):
        """
        Label candidates given variants from a VCF
        :param positional_vcf: IN/DEL/SNPs separated into VCF records, expanded into a 1-to-1 ref_pos:variant allele
        :param candidate_sites: Candidates
        :return: List of labeled candidate sites
        """
        # create a log of whether variant positions have been matched to a candidate
        validated_vcf_positions = {key: 0 for key in positional_vcf.keys()}

        # list of all labeled candidates
        all_labeled_candidates = []

        # for each candidate
        for candidate_site in candidate_sites:
            allele_start = candidate_site[START]
            allele_stop = candidate_site[STOP]
            alleles_insert = candidate_site[IN_ALLELES]
            alleles_snp = candidate_site[SNP_ALLELES]
            alleles_del = candidate_site[DEL_ALLELES]

            ref_sequence = self.fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                           start=allele_start,
                                                           stop=allele_stop + 1)

            # test the alleles across IN, DEL, and SNP variant dictionaries
            genotypes = self._get_all_genotype_labels(positional_vcf=positional_vcf,
                                                      start=allele_start,
                                                      stop=allele_stop,
                                                      ref_seq=ref_sequence,
                                                      alleles_snp=alleles_snp,
                                                      alleles_insert=alleles_insert,
                                                      alleles_del=alleles_del)

            # get a list of attributes that can be saved
            all_labeled_candidates.extend(self._generate_list(chromosome_name, allele_start, allele_stop, alleles_snp,
                                                              alleles_insert, alleles_del, ref_sequence, genotypes))

            if DEBUG_PRINT_ALL:
                print()
                print(allele_start, allele_stop, ref_sequence)
                print("snp alleles:        ", alleles_snp)
                print("insert alleles: ", alleles_insert)
                print("delete alleles: ", alleles_del)
                print("insert gts:     ", genotypes[IN])
                print("del gts:        ", genotypes[DEL])
                print("snp gts:        ", genotypes[SNP])

            # test if the site has any true variants (not Hom)
            # if self._is_position_supported(genotypes):
            #     if allele_start in validated_vcf_positions:
            #         validated_vcf_positions[allele_start] = 1

        # see if any sites weren't supported by candidates
        # n_unsupported = 0
        # for pos in validated_vcf_positions:
        #     if validated_vcf_positions[pos] == 0:
        #         n_unsupported += 1
        #         print("\nWARNING: Unsupported VCF position: ", pos)
        #         print("\tRecord: ", positional_vcf[pos])

        # print("missed vcf positions: ", float(len(positional_vcf))/n_unsupported)

        return all_labeled_candidates
