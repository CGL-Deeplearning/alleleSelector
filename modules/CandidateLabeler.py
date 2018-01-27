from collections import Counter


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
        ref_seq = rec.ref[0]
        alt_seq = rec.alt[:-len(rec.ref)+1]
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
        position_based_vcf = {}
        for variant_pos in variants:
            for variant in variants[variant_pos]:
                if variant.genotype_class == 'IN':
                    pos, ref_seq, alt_seq, genotype = self._handle_insert(variant)
                    if pos not in position_based_vcf:
                        position_based_vcf[pos] = list()
                    position_based_vcf[pos].append((pos, ref_seq, alt_seq, genotype))

                if variant.genotype_class == 'DEL':
                    list_of_records = self._handle_delete(variant)
                    for record in list_of_records:
                        pos, ref_seq, alt_seq, genotype = record

                        if pos not in position_based_vcf:
                            position_based_vcf[pos] = list()

                        position_based_vcf[pos].append((pos, ref_seq, alt_seq, genotype))

                if variant.genotype_class == 'SNP':
                    pos, ref_seq, alt_seq, genotype = variant.pos + self.vcf_offset, variant.ref, variant.alt, variant.type
                    if pos not in position_based_vcf:
                        position_based_vcf[pos] = list()
                    position_based_vcf[pos].append((pos, ref_seq, alt_seq, genotype))

        return position_based_vcf

    def get_label_of_allele(self, position_based_vcf, candidate_allele):
        pos_start, pos_stop, ref, alt = candidate_allele
        gt = 'Hom'
        for i in range(pos_start, pos_stop+1):
            if i in position_based_vcf.keys():
                records = position_based_vcf[i]
                for record in records:
                    # print(record)
                    rec_alt = record[2]
                    if rec_alt == alt:
                        gt = record[3]
        return gt

    def get_labeled_candidates2(self, variants, candidate_sites):
        position_based_vcf = self.generate_position_based_vcf(variants)

        labled_candidates_list = list()
        for candidate_site in candidate_sites:
            chromosome_name = candidate_site.chromosome_name
            allele_start = candidate_site.window_start
            allele_stop = candidate_site.window_end

            ref_sequence = self.fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                           start=allele_start,
                                                           stop=allele_stop + 1)

            alt1_sequence, alt2_sequence = self._select_alleles(allele_list=candidate_site.candidate_alleles,
                                                                ref_sequence=ref_sequence)
            # print(alt1_sequence, alt2_sequence)
            gt1 = self.get_label_of_allele(position_based_vcf, (allele_start, allele_stop, ref_sequence, alt1_sequence))
            labled_candidates_list.append((allele_start, allele_stop, ref_sequence, alt1_sequence, gt1))
            if alt2_sequence is not None:
                gt2 = self.get_label_of_allele(position_based_vcf, (allele_start, allele_stop, ref_sequence, alt2_sequence))
                labled_candidates_list.append((allele_start, allele_stop, ref_sequence, alt2_sequence, gt2))
            if alt1_sequence is None and alt2_sequence is None:
                print(chromosome_name, allele_start, allele_stop, ref_sequence)

        return labled_candidates_list


    def get_labeled_candidates(self, variants, candidate_sites):
        """
        Finds the candidate allele window that contains each variant
        :param variants: dictionary of start:VariantRecord
        :param candidate_sites: dictionary of candidate sites
        :return:
        """
        labeled_sites = dict()

        v = 0
        a = 0

        variant_keys = sorted(list(variants.keys()))
        prev_start = None

        while v < len(variant_keys) and a < len(candidate_sites):
            chromosome_name = candidate_sites[a].chromosome_name
            allele_start = candidate_sites[a].window_start
            allele_stop = candidate_sites[a].window_end
            variant_start = variant_keys[v]
            site_variants = variants[variant_start]

            ref_sequence = self.fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                           start=allele_start,
                                                           stop=allele_stop+1)

            alt1_sequence, alt2_sequence = self._select_alleles(allele_list=candidate_sites[a].candidate_alleles,
                                                                ref_sequence=ref_sequence)

            # --- MATCH ---
            # if variant falls in the candidate allele's window, update the full and abridged dictionary
            if allele_start <= (variant_start + self.vcf_offset) <= allele_stop:
                output = self._update_candidate_site(allele=candidate_sites[a],
                                                     site_variants=site_variants,
                                                     alt1_sequence=alt1_sequence,
                                                     alt2_sequence=alt2_sequence)

                candidate_sites[a], labeled_sites[allele_start], is_multi_vcf_site = output

                # if a site already has a vcf entry, iterate through the entries and compare their genotypes
                if is_multi_vcf_site:
                    genotypes = list()
                    for record in candidate_sites[a].vcf_records:
                        for variant in record:
                            genotypes.append(variant.type)

                    # if there are non-matching genotypes in the same window, print site + warning
                    if len(set(genotypes)) > 1:
                        print("\nWARNING: Multiple non-matching genotypes in candidate window:")
                        self._print_site_data(candidate_site=candidate_sites[a],
                                              alt1_sequence=alt1_sequence,
                                              alt2_sequence=alt2_sequence,
                                              print_variants=True)
                v += 1

            # --- NO VARIANT ---
            # if variant comes after the candidate window, move to next candidate window
            elif allele_stop < (variant_start + self.vcf_offset):

                # update the candidate dictionaries as hom-ref if this candidate window has no true variants
                if prev_start != allele_start:
                    output = self._update_candidate_site(allele=candidate_sites[a],
                                                         site_variants=None,
                                                         alt1_sequence=alt1_sequence,
                                                         alt2_sequence=alt2_sequence)

                    candidate_sites[a], labeled_sites[allele_start], is_multi_vcf_site = output
                a += 1

            # --- NO ALLELE ---
            # if variant has been passed already, it was not supported in the reads as a "candidate"... shouldn't happen
            elif (variant_start + self.vcf_offset) < allele_start:
                print("\nWARNING: variant not supported by candidates: ")
                for rec in site_variants:
                    print(rec)
                self._print_site_data(candidate_site=candidate_sites[a],
                                      alt1_sequence=alt1_sequence,
                                      alt2_sequence=alt2_sequence,
                                      print_variants=False)
                v += 1
            else:
                print("INVALID COORDINATE RELATIONSHIP: ", variant_start, allele_start, allele_stop)    # ooh scandalous

            prev_start = allele_start

        return labeled_sites

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

        del allele_counter[ref_sequence] # don't consider the reference allele

        # keep only the top 2 most frequent alleles (assuming diploidy)
        top_2_alleles = list(dict(allele_counter.most_common(2)).keys())

        top_2_alleles += [None]*(2-len(top_2_alleles))

        return top_2_alleles

    def _print_site_data(self, candidate_site, alt1_sequence, alt2_sequence, print_variants):
        allele_start = candidate_site.window_start
        allele_stop = candidate_site.window_end
        ref_sequence = candidate_site.reference_sequence

        print("\nCandidate: ")
        print("window: ", allele_start, allele_stop)
        print("ref:    ", ref_sequence)
        print("alts:   ", alt1_sequence, alt2_sequence)
        if print_variants:
            for record in candidate_site.genotype:
                for variant in record:
                    print("\nVariant:  ")
                    print("position  ", variant.pos-1)
                    print("ref:      ", variant.ref)
                    print("alt:      ", variant.alt)
                    print("genotype: ", variant.genotype_class)
                    print("type      ", variant.type)

    def _update_candidate_site(self, allele, site_variants, alt1_sequence, alt2_sequence):
        if allele.genotype_reported is False:
            allele.genotype_reported = True
            allele.genotype = list()
            is_multi_vcf_site = False
        else:
            allele.vcf_updated = True
            is_multi_vcf_site = True

        if site_variants is not None:
            genotype = site_variants[0].type            # use the first genotype, assuming it is the same as any others
            allele.genotype.append(site_variants)    # each entry pertains to 1 or 2 records from a given coordinate
        else:
            genotype = "Hom"

        site = dict()
        site["chromosome_name"] = allele.chromosome_name
        site["window_start"] = allele.window_start
        site["window_end"] = allele.window_end
        site["ref_sequence"] = allele.reference_sequence
        site["alt1_sequence"] = alt1_sequence
        site["alt2_sequence"] = alt2_sequence
        site["genotype"] = genotype

        return allele, site, is_multi_vcf_site

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
