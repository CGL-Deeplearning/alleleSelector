from collections import Counter
from collections import defaultdict

class AlleleSelector:
    def __init__(self, allele_file_handler, vcf_handler, fasta_handler):
        '''
        :param allele_file_handler: module that loads candidate alleles from JSON file
        :param vcf_handler: module that fetches variant records from a VCF file
        :param fasta_handler: module that fetches reference sequence substrings from a FASTA file
        '''
        self.allele_file_handler = allele_file_handler
        self.vcf_handler = vcf_handler
        self.fasta_handler = fasta_handler

        self.vcf_offset = -1    # pysam vcf coords are 1-based ... >:[ ... this is what Kishwar wanted
        self.delete_char = '*'

        self.chromosome_name = None
        self.start = None
        self.stop = None

    def parse_candidate_sites(self):
        '''
        Iterate through allele site dictionary and find the top 2 alleles which may represent the reference, alt1, or
        alt2 allele
        :return:
        '''
        # get dictionary of sites containing observed alleles
        candidate_sites = self.allele_file_handler.get_allele_site_dictionary()["all_candidates"]

        # find start and stop position of region that covers all candidates
        start_record = candidate_sites[0]
        stop_record = candidate_sites[-1]

        self.chromosome_name = start_record["chromosome_name"]
        self.start = start_record["window_start"]
        self.stop = stop_record["window_end"]

        # iterate sites and return allele candidates
        # alleles_by_site = dict()
        for r,record in enumerate(candidate_sites):
            # start = record["window_start"]
            # stop = record["window_end"]
            alleles = self.select_alleles(record["candidate_alleles"])
            record["selected_alleles"] = alleles

            # add the list of sequences to the allele sequence dictionary with key=window_start
            candidate_sites[r] = record

            # print(candidate_sites[r])

        return candidate_sites

    def select_alleles(self, allele_list):
        '''
        Given a dictionary of allele objects, naively find the top 2 represented sequences and return a dictionary
        that assigns read_ids to each of the 2 alleles.
        :param allele_list: the list of alleles found at a given site
        :return: top_2_alleles: the most frequent two alleles
        '''
        alleles_sequences = [allele["allele_sequence"] for allele in allele_list]

        # find allele distribution by counting instances of each
        allele_counter = Counter(alleles_sequences)
        allele_frequencies = sorted(dict(allele_counter).items(), key=lambda x: x[1], reverse=True)

        # keep only the top 2 most frequent alleles (assuming diploidy)
        top_2_alleles = allele_counter.most_common(2)
        # print("----")
        # print(allele_frequencies)
        # print(top_2_alleles)
        return top_2_alleles

    def test_alleles(self, candidate_sites):
        '''
        Check the vcf at the given site and compare with top 2 alleles to determine whether verified alleles correspond
        with observed alleles
        :param alleles: dictionary of allele sites
        :return:
        '''

        self.vcf_handler.populate_dictionary(contig=self.chromosome_name,
                                             start_pos=self.start,
                                             end_pos=self.stop,
                                             hom_filter=True)

        variants = self.vcf_handler.get_variant_dictionary()

        self.match_vcf_variants_to_alleles(variants=variants,candidate_sites=candidate_sites)

    def match_vcf_variants_to_alleles(self,variants,candidate_sites):
        '''
        Finds the candidate allele window that contains each variant
        :param variants: dictionary of start:VariantRecord
        :param alleles: dictionary of start:allele_sequence
        :return:
        '''

        # expand variants into a single list
        # variants = [variant for variant_list in variants_dict.values() for variant in variant_list]
        v = 0
        a = 0

        prev_match = False
        vcf_matched_alleles = defaultdict(dict)

        variant_keys = list(variants.keys())

        print(len(candidate_sites))
        print(len(variant_keys))

        prev_variant_start = None
        prev_allele_start = None
        prev_allele_stop = None

        while v < len(variant_keys) and a < len(candidate_sites):
            variant_start = variant_keys[v]
            allele_start = candidate_sites[a]["window_start"]
            allele_stop = candidate_sites[a]["window_end"]

            # MATCH:
            # if variant falls in the candidate allele's window, add it to the unified dictionary
            if allele_start <= (variant_start + self.vcf_offset) < allele_stop:
                # print([[variant.alt, variant.ref] for variant in variants[variant_start]])
                v += 1

            # NO VARIANT:
            # if variant comes after the candidate window, move to next candidate window
            elif allele_start < allele_stop < variant_start:
                a += 1

            # NO ALLELE:
            # if variant has been passed already, it was not supported in the reads as a "candidate"... shouldn't happen
            elif variant_start < allele_start < allele_stop:
                print("WARNING: variant not supported: ",
                      variant_start,
                      [[variant.alt, variant.ref] for variant in variants[variant_start]])
                v += 1

            else:
                print("#################### ERROR: invalid vcf/allele position relationship #######################")

            # print(candidate_sites[a]["selected_alleles"])
            # print("a: ", allele_start, allele_stop, "\nv: ", (variant_start + self.vcf_offset))

            if prev_match:
                variant = variants[prev_variant_start]

                allele_observed = candidate_sites[a-1]["selected_alleles"]
                allele_vcf = self.get_merged_variants_for_site(variant, prev_allele_start, prev_allele_stop)

                vcf_matched_alleles[prev_allele_start]["observed"] = allele_observed
                vcf_matched_alleles[prev_allele_start]["vcf"] = allele_vcf

            prev_variant_start = variant_start
            prev_allele_start = allele_start
            prev_allele_stop = allele_stop



    def get_merged_variants_for_site(self, variant, window_start, window_stop):
        '''
        Given a candidate window and any true variants that might be found within it, obtain the potential haplotypes
        that span the window and determine if they are supported by the alleles in this read.
        :param variant:
        :param window_start:
        :param window_stop:
        :return:
        '''
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

    def get_substituted_sequence(self, variant, ref_sequence, window_start, window_stop):
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

# [('TAA', 13), ('TAAAAAA', 10), ('TAAAAA', 10), ('AAA', 3)]
