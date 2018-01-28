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

DEBUG_FREQUENCIES = False
DEBUG_PRINT_ALL = False

class CandidateLabeler:
    def __init__(self, fasta_handler):
        '''
        :param allele_file_handler: module that loads candidate alleles from JSON file
        :param vcf_handler: module that fetches variant records from a VCF file
        :param fasta_handler: module that fetches reference sequence substrings from a FASTA file
        '''
        self.fasta_handler = fasta_handler      # unfortunately need to query reference sequence on a per-site basis

        self.vcf_offset = -1                    # pysam vcf coords are 1-based ... >:[ ... this is what Kishwar wanted
        self.delete_char = '*'

    def _handle_insert(self, rec):
        ref_seq = rec.ref   # no change necessary
        alt_seq = rec.alt   # no change necessary

        pos = rec.pos + self.vcf_offset
        return pos, ref_seq, alt_seq, rec.type

    def _handle_delete(self, rec):
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
        position_based_vcf_in = {}
        position_based_vcf_del = {}
        position_based_vcf_snp = {}
        for variant_pos in variants:
            for variant in variants[variant_pos]:
                if variant.genotype_class == 'IN':
                    pos, ref_seq, alt_seq, genotype = self._handle_insert(variant)

                    if pos not in position_based_vcf_in:
                        position_based_vcf_in[pos] = list()

                    position_based_vcf_in[pos].append((pos, ref_seq, alt_seq, genotype, variant.genotype_class))

                if variant.genotype_class == 'DEL':
                    list_of_records = self._handle_delete(variant)
                    for record in list_of_records:
                        pos, ref_seq, alt_seq, genotype = record

                        if pos not in position_based_vcf_del:
                            position_based_vcf_del[pos] = list()

                        position_based_vcf_del[pos].append((pos, ref_seq, alt_seq, genotype, variant.genotype_class))

                if variant.genotype_class == 'SNP':
                    pos, ref_seq, alt_seq, genotype = variant.pos + self.vcf_offset, variant.ref, variant.alt, variant.type

                    if pos not in position_based_vcf_snp:
                        position_based_vcf_snp[pos] = list()

                    position_based_vcf_snp[pos].append((pos, ref_seq, alt_seq, genotype, variant.genotype_class))

        return position_based_vcf_in, position_based_vcf_del, position_based_vcf_snp

    def get_label_of_allele(self, positional_vcfs, variant_type, candidate_allele):
        pos_start, pos_stop, ref, alt = candidate_allele
        gt = 'Hom'
        for i in range(pos_start, pos_stop+1):
            if i in positional_vcfs[variant_type].keys():
                records = positional_vcfs[variant_type][i]
                for record in records:

                    rec_alt = record[2]
                    if rec_alt == alt:
                        gt = record[3]

        return gt

    def _get_all_genotype_labels(self, positional_vcfs, start, stop, ref_seq, alleles, alleles_insert):
        alt1_seq, alt2_seq = alleles
        alt1_seq_insert, alt2_seq_insert = alleles_insert

        gt1_in = self.get_label_of_allele(positional_vcfs, "IN", (start,stop,ref_seq,alt1_seq_insert))
        gt2_in = self.get_label_of_allele(positional_vcfs, "IN", (start,stop,ref_seq,alt2_seq_insert))
        gt1_del = self.get_label_of_allele(positional_vcfs, "DEL", (start,stop,ref_seq,alt1_seq))
        gt2_del = self.get_label_of_allele(positional_vcfs, "DEL", (start,stop,ref_seq,alt2_seq))
        gt1_snp = self.get_label_of_allele(positional_vcfs, "SNP", (start,stop,ref_seq,alt1_seq))
        gt2_snp = self.get_label_of_allele(positional_vcfs, "SNP", (start,stop,ref_seq,alt2_seq))

        return {"IN":[gt1_in, gt2_in], "DEL":[gt1_del, gt2_del], "SNP":[gt1_snp, gt2_snp]}

    def _is_supported(self, genotypes):
        gt1,gt2 = genotypes
        return not (gt1 == "Hom" and (gt2 == "Hom" or gt2 is None))

    def _is_position_supported(self, genotypes):
        in_supported = self._is_supported(genotypes["IN"])
        del_supported = self._is_supported(genotypes["DEL"])
        snp_supported = self._is_supported(genotypes["SNP"])

        return in_supported or del_supported or snp_supported

    def get_labeled_candidates2(self, variants, candidate_sites):
        types = ["IN","DEL","SNP"]

        # get separate positional variant dictionaries for IN, DEL, and SNP
        positional_vcfs = {key:value for key,value in zip(types,self.generate_position_based_vcf(variants))}

        # get the unified set of all variant positions
        pos_in, pos_del, pos_snp = [list(vcf_dict.keys()) for vcf_dict in positional_vcfs.values()]
        all_positions = pos_in + pos_del + pos_snp
        validated_vcf_positions = {key:False for key in all_positions}

        # labeled_candidates_by_pos = dict()
        tsv_output = ''

        print(validated_vcf_positions.keys())

        # n_true_positive = 0
        # n_true_negative = 0
        # n_unsupported_vcf = 0

        for candidate_site in candidate_sites:
            chromosome_name = candidate_site.chromosome_name
            allele_start = candidate_site.window_start
            allele_stop = candidate_site.window_end

            ref_sequence = self.fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                           start=allele_start,
                                                           stop=allele_stop + 1)

            if DEBUG_FREQUENCIES:
                frequencies = self._get_allele_frequency_vector(candidate_site.candidate_alleles,ref_sequence)

            # get selected insert and non-insert alleles (alt1,alt2) (alt1_insert,alt2_insert)
            alleles_insert, alleles = self._select_alleles(allele_list=candidate_site.candidate_alleles,
                                                           ref_sequence=ref_sequence)

            # test the alleles across IN, DEL, and SNP variant dictionaries
            genotypes = self._get_all_genotype_labels(positional_vcfs=positional_vcfs,
                                                      start=allele_start,
                                                      stop=allele_stop,
                                                      ref_seq=ref_sequence,
                                                      alleles=alleles,
                                                      alleles_insert=alleles_insert)

            # create a record for this candidate window
            # labeled_candidates_by_pos[allele_start] = {"window_start": allele_start,
            #                                            "window_end": allele_stop,
            #                                            "ref_sequence": ref_sequence,
            #                                            "alleles":{"SNP": alleles,
            #                                                       "DEL": alleles,
            #                                                       "IN" : alleles_insert},
            #                                            "genotypes": genotypes}               # same format as alleles

            # get the BED format string for this candidate
            tsv_line = self._generate_tsv_string(chromosome_name,
                                                 allele_start,
                                                 allele_stop,
                                                 alleles,
                                                 alleles_insert,
                                                 ref_sequence,
                                                 genotypes)

            tsv_output += tsv_line + '\n'

            if DEBUG_PRINT_ALL:
                print()
                print(allele_start,allele_stop,ref_sequence)
                if DEBUG_FREQUENCIES:
                    print("frequencies:    ",frequencies)
                print("alleles:        ",alleles)
                print("insert alleles: ",alleles_insert)
                print("insert gts:     ", genotypes["IN"])
                print("del gts:        ", genotypes["DEL"])
                print("snp gts:        ", genotypes["SNP"])

            # test if the site has any true variants (not Hom)
            if self._is_position_supported(genotypes):
                validated_vcf_positions[allele_start] = True

        # see if any sites weren't supported by candidates
        for pos in validated_vcf_positions:
            if not validated_vcf_positions[pos]:
                print("\nWARNING: Unsupported VCF position: ", positional_vcfs[pos])
                print("\tRecord: ", positional_vcfs[pos])

        return tsv_output

    def _generate_tsv_string(self, chromosome_name, start, stop, alleles, alleles_insert, ref_seq, genotypes):
        # header: chr start stop ref_seq alt1 alt2 gt1 gt2

        allele_string = None
        gt_string = None
        variant_type = None
        is_insert = None
        supported = set()

        if self._is_supported(genotypes["SNP"]) and self._is_supported(genotypes["DEL"]):
            variant_type = "SNP/DEL"
            is_insert = False
            supported.add(variant_type)
            alleles = [allele if allele is not None else "None" for allele in alleles]
            allele_string = '\t'.join(alleles)
            merged_genotype = self._merge_snp_and_del(alleles,genotypes["SNP"],genotypes["DEL"])
            gt_string = '\t'.join(merged_genotype)

        elif self._is_supported(genotypes["SNP"]):
            variant_type = "SNP"
            is_insert = False
            supported.add(variant_type)
            alleles = [allele if allele is not None else "None" for allele in alleles]
            allele_string = '\t'.join(alleles)
            gt_string = '\t'.join(genotypes["SNP"])

        elif self._is_supported(genotypes["IN"]):
            variant_type = "IN"
            is_insert = True
            supported.add(variant_type)
            alleles_insert = [allele if allele is not None else "None" for allele in alleles_insert]
            allele_string = '\t'.join(alleles_insert)
            gt_string = '\t'.join(genotypes["IN"])

        elif self._is_supported(genotypes["DEL"]):
            variant_type = "DEL"
            is_insert = False
            supported.add(variant_type)
            alleles = [allele if allele is not None else "None" for allele in alleles]
            allele_string = '\t'.join(alleles)
            gt_string = '\t'.join(genotypes["DEL"])

        if len(supported) == 0:     # no true variant
            variant_type = "SNP"
            is_insert = False
            allele_string = '\t'.join(["None","None"])
            gt_string = '\t'.join(["Hom", "Hom"])

        # elif len(supported) > 1:
        #     print("WARNING: more than one supported variant type: ", chromosome_name, start, stop, list(supported))

        data = [chromosome_name, start, stop, is_insert, ref_seq, allele_string, gt_string]
        output_string = '\t'.join(list(map(str, data)))
        # print(output_string)

        return output_string

    def _merge_snp_and_del(self, alleles, snp_genotypes, del_genotypes):
        output_genotypes = ["Hom","Hom"]
        del_genotype_index = alleles.index('*')
        snp_genotype_index = 1 - del_genotype_index

        output_genotypes[del_genotype_index] = del_genotypes[del_genotype_index]
        output_genotypes[snp_genotype_index] = snp_genotypes[snp_genotype_index]

        return output_genotypes

    def test_merge_snp_and_del(self):
        alleles = ["A","*"]
        del_genotypes = ['Hom','Het']
        snp_genotypes = ['Het','Hom']

        print(self._merge_snp_and_del(alleles, snp_genotypes, del_genotypes))

    def _select_alleles(self, allele_list, ref_sequence):
        '''
        Given a dictionary of allele objects, naively find the top 2 represented sequences and return a dictionary
        that assigns read_ids to each of the 2 alleles.
        :param allele_list: the list of alleles found at a given site
        :param ref_sequence: reference sequence for the candidate window
        :return: top_2_alleles: the most frequent two alleles
        '''

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

    def _get_allele_frequency_vector(self, allele_list, ref_sequence, vector_length=8, normalize_by_depth=True):
        '''
        Find the sorted frequency vector: the # of alleles observed at a site, sorted in descending order
        :param allele_list: the list of alleles found at a given site
        :return: top_2_alleles: the most frequent two alleles
        '''
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

    def _get_merged_variants_for_site(self, variant, window_start, window_stop):
        """
        unfinished testing method ... do the vcf sequences actually match the candidate site alleles?

        Given a candidate window and any true variants that might be found within it, obtain the potential haplotypes
        that span the window and determine if they are supported by the alleles in this read.
        :param variant:
        :param window_start:
        :param window_stop:
        :return:
        """
        variant_start = variant.pos
        variant_stop = variant.pos + len(variant.ref)

        ref_sequence = self.fasta_handler.get_sequence(chromosome_name=self.chromosome_name,
                                                       start=window_start-1,
                                                       stop=window_stop+1)

        # if more than one variant, find window split points

        # build allele segments across split points

        # merge segments into full window alleles

        # compare vcf vs observed alleles

        return

    def _get_substituted_sequence(self, variant, ref_sequence, window_start, window_stop):
        '''
        unfinished testing method, to be used with get_merged_variants_for_site
        '''
        variant_class = variant.genotype_class
        variant_start = variant.pos

        if variant_class == "SNP":
            index = window_start - (variant_start + self.vcf_offset)
            alt_sequence = variant.alt

            return ref_sequence[:index] + alt_sequence + ref_sequence[index+1:]

        elif variant_class == "IN":
            # ignore anchor position
            index = window_start - (variant_start + self.vcf_offset) + 1
            alt_sequence = variant.alt[1:]

            return ref_sequence[:index] + alt_sequence + ref_sequence[index+1:]

        elif variant_class == "DEL":
            # ignore anchor position
            index = window_start - (variant_start + self.vcf_offset) + 1
            length = len(variant.ref) - len(variant.alt)
            alt_sequence = self.delete_char*length

            return ref_sequence[:index] + alt_sequence + ref_sequence[index+1:]

        else:
            raise ("ERROR: Invalid variant class: %s"%variant_class)
