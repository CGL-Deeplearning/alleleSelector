from collections import Counter
from collections import defaultdict
import numpy
import argparse
from os import listdir, mkdir
from os.path import isfile, join, exists


class AlleleSelector:
    def __init__(self, fasta_handler):
        '''
        :param allele_file_handler: module that loads candidate alleles from JSON file
        :param vcf_handler: module that fetches variant records from a VCF file
        :param fasta_handler: module that fetches reference sequence substrings from a FASTA file
        '''
        self.fasta_handler = fasta_handler      # unfortunately need to query reference sequence on a per-site basis

        self.vcf_offset = -1                    # pysam vcf coords are 1-based ... >:[ ... this is what Kishwar wanted
        self.delete_char = '*'

        self.chromosome_name = None
        self.start = None
        self.stop = None

    def parse_candidate_sites(self, candidate_sites):
        '''
        Prediction pipeline:
        Iterate through allele site dictionary and find the top 2 alleles which may represent the reference, alt1, or
        alt2 allele
        :return:
        '''
        # get dictionary of sites containing observed alleles

        # find start and stop position of region that covers all candidates
        start_record = candidate_sites[0]
        stop_record = candidate_sites[-1]

        self.chromosome_name = start_record["chromosome_name"]
        self.start = start_record["window_start"]
        self.stop = stop_record["window_end"]

        # iterate sites and return allele candidates
        for r,record in enumerate(candidate_sites):
            # start = record["window_start"]
            # stop = record["window_end"]
            alleles = self.select_alleles(record["candidate_alleles"])
            record["selected_alleles"] = alleles

            # add the list of sequences to the allele sequence dictionary with key=window_start
            candidate_sites[r] = record

        return candidate_sites

    def get_allele_frequency_vector(self, allele_list, ref_sequence, vector_length=8, normalize_by_depth=True):
        '''
        Find the sorted frequency vector: the # of alleles observed at a site, sorted in descending order
        :param allele_list: the list of alleles found at a given site
        :return: top_2_alleles: the most frequent two alleles
        '''
        alleles_sequences = [allele["allele_sequence"] for allele in allele_list]

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
        frequencies = numpy.array(frequencies)/len(alleles_sequences)

        return frequencies

    def select_alleles(self, allele_list, ref_sequence):
        '''
        Given a dictionary of allele objects, naively find the top 2 represented sequences and return a dictionary
        that assigns read_ids to each of the 2 alleles.
        :param allele_list: the list of alleles found at a given site
        :param ref_sequence: reference sequence for the candidate window
        :return: top_2_alleles: the most frequent two alleles
        '''
        alleles_sequences = [allele["allele_sequence"] for allele in allele_list]

        # find allele distribution by counting instances of each
        allele_counter = Counter(alleles_sequences)

        del allele_counter[ref_sequence] # don't consider the reference allele

        # keep only the top 2 most frequent alleles (assuming diploidy)
        top_2_alleles = list(dict(allele_counter.most_common(2)).keys())

        top_2_alleles += [None]*(2-len(top_2_alleles))

        return top_2_alleles

    def test_allele_distributions(self, vcf_handler, candidate_sites):
        '''
        Check the vcf at the given site and compare with top 2 alleles to determine whether verified alleles correspond
        with observed alleles
        :param alleles: dictionary of allele sites
        :return:
        '''
        # find start and stop position of region that covers all candidates
        start_record = candidate_sites[0]
        stop_record = candidate_sites[-1]

        self.chromosome_name = start_record["chromosome_name"]
        self.start = start_record["window_start"]
        self.stop = stop_record["window_end"]

        vcf_handler.populate_dictionary(contig=self.chromosome_name,
                                        start_pos=self.start,
                                        end_pos=self.stop,
                                        hom_filter=True)

        variants = vcf_handler.get_variant_dictionary()

        self.bin_distributions_by_variant_status(variants=variants, candidate_sites=candidate_sites, vector_cutoff=8)

    def get_labeled_candidate_sites(self, vcf_handler, candidate_sites, filter_hom_ref=False):
        '''
        Check the vcf at the given site and compare with top 2 alleles to determine whether verified alleles correspond
        with observed alleles
        :param vcf_handler: vcf querying module (assumed to use 1-based coords)
        :param candidate_sites: dictionary of allele sites
        :param filter_hom_ref: turn on to exclude VCF hom-ref labels from output
        :return:
        '''
        # find start and stop position of region that covers all candidates
        start_record = candidate_sites[0]
        stop_record = candidate_sites[-1]

        self.chromosome_name = start_record["chromosome_name"]
        self.start = start_record["window_start"]
        self.stop = stop_record["window_end"]

        # get dictionary of variant records for full region
        vcf_handler.populate_dictionary(contig=self.chromosome_name,
                                        start_pos=self.start,
                                        end_pos=self.stop,
                                        hom_filter=filter_hom_ref)

        variants = vcf_handler.get_variant_dictionary()

        # iterate through variants and candidates to match them and compile training labels for candidates
        labeled_candidates = self.get_labeled_candidates(variants=variants, candidate_sites=candidate_sites)

        return labeled_candidates

    def bin_distributions_by_variant_status(self, variants, candidate_sites, vector_cutoff):
        '''
        Finds the candidate allele window that contains each variant
        :param variants: dictionary of start:VariantRecord
        :param candidate_sites: list of candidate site dictionaries containing data about the site
        :param vector_cutoff: cutoff for size of vector containing alleles at each site
        :return:
        '''

        print("\nregion:          ", self.chromosome_name, self.start, self.stop)
        print("candidate sites: ", len(candidate_sites))
        print("variant sites:   ", len(variants))

        numpy.set_printoptions(suppress=True, precision=2)

        variant_distribution = numpy.zeros(vector_cutoff)
        error_distribution = numpy.zeros(vector_cutoff)
        k_var = 0
        k_err = 0

        # variants = [variant for variant_list in variants_dict.values() for variant in variant_list]
        v = 0
        a = 0

        variant_keys = list(variants.keys())
        frequency_vector = None
        ref_sequence = None
        prev_start = None

        while v < len(variant_keys) and a < len(candidate_sites):
            variant_start = variant_keys[v]
            allele_start = candidate_sites[a]["window_start"]
            allele_stop = candidate_sites[a]["window_end"]

            if frequency_vector is None:
                ref_sequence = self.fasta_handler.get_sequence(chromosome_name=self.chromosome_name,
                                                               start=allele_start,
                                                               stop=allele_stop+1)

                frequency_vector = self.get_allele_frequency_vector(allele_list=candidate_sites[a]["candidate_alleles"],
                                                                    ref_sequence=ref_sequence)

            # MATCH:
            # if variant falls in the candidate allele's window, add it to the unified dictionary
            if allele_start <= (variant_start + self.vcf_offset) <= allele_stop:
                variant_distribution += frequency_vector
                k_var += 1
                v += 1

            # NO VARIANT:
            # if variant comes after the candidate window, move to next candidate window
            elif allele_start < allele_stop < (variant_start + self.vcf_offset):
                if prev_start != allele_start:
                    error_distribution += frequency_vector
                    frequency_vector = None
                    k_err += 1
                a += 1

            # NO ALLELE:
            # if variant has been passed already, it was not supported in the reads as a "candidate"... shouldn't happen
            elif (variant_start + self.vcf_offset) < allele_start < allele_stop:
                v += 1
            else:
                print("INVALID COORDINATE RELATIONSHIP: ")


            prev_start = allele_start

        print("error sites:     ", k_err, error_distribution/k_err)
        print("variant sites:   ", k_var, variant_distribution/k_var)

    def get_labeled_candidates(self, variants, candidate_sites):
        '''
        Finds the candidate allele window that contains each variant
        :param variants: dictionary of start:VariantRecord
        :param candidate_sites: dictionary of candidate sites
        :return:
        '''
        labeled_sites = dict()

        # variants = [variant for variant_list in variants_dict.values() for variant in variant_list]
        v = 0
        a = 0

        variant_keys = list(variants.keys())
        prev_start = None

        while v < len(variant_keys) and a < len(candidate_sites):
            allele_start = candidate_sites[a]["window_start"]
            allele_stop = candidate_sites[a]["window_end"]
            variant_start = variant_keys[v]
            site_variants = variants[variant_start]

            ref_sequence = self.fasta_handler.get_sequence(chromosome_name=self.chromosome_name,
                                                           start=allele_start,
                                                           stop=allele_stop+1)

            alt1_sequence, alt2_sequence = self.select_alleles(allele_list=candidate_sites[a]["candidate_alleles"],
                                                               ref_sequence=ref_sequence)

            # --- MATCH ---
            # if variant falls in the candidate allele's window, update the full and abridged dictionary
            if allele_start <= (variant_start + self.vcf_offset) <= allele_stop:
                output = self.update_candidate_site(allele=candidate_sites[a],
                                                    site_variants=site_variants,
                                                    alt1_sequence=alt1_sequence,
                                                    alt2_sequence=alt2_sequence)

                candidate_sites[a], labeled_sites[allele_start], is_multi_vcf_site = output

                # if a site already has a vcf entry, iterate through the entries and compare their genotypes
                if is_multi_vcf_site:
                    genotypes = list()
                    for record in candidate_sites[a]["vcf_records"]:
                        for variant in record:
                            genotypes.append(variant.type)

                    # if there are non-matching genotypes in the same window, print site + warning
                    if len(set(genotypes)) > 1:
                        print("\nWARNING: Multiple non-matching genotypes in candidate window:")
                        self.print_site_data(candidate_site=candidate_sites[a],
                                             alt1_sequence=alt1_sequence,
                                             alt2_sequence=alt2_sequence)
                v += 1

            # --- NO VARIANT ---
            # if variant comes after the candidate window, move to next candidate window
            elif allele_start < allele_stop < (variant_start + self.vcf_offset):

                # update the candidate dictionaries as hom-ref if this candidate window has no true variants
                if prev_start != allele_start:
                    output = self.update_candidate_site(allele=candidate_sites[a],
                                                        site_variants=None,
                                                        alt1_sequence=alt1_sequence,
                                                        alt2_sequence=alt2_sequence)

                    candidate_sites[a], labeled_sites[allele_start], is_multi_vcf_site = output
                a += 1

            # --- NO ALLELE ---
            # if variant has been passed already, it was not supported in the reads as a "candidate"... shouldn't happen
            elif (variant_start + self.vcf_offset) < allele_start < allele_stop:
                print("\nWARNING: variant not supported by candidates: ")
                self.print_site_data(candidate_site=candidate_sites[a],
                                     alt1_sequence=alt1_sequence,
                                     alt2_sequence=alt2_sequence)
                v += 1
            else:
                print("INVALID COORDINATE RELATIONSHIP: ", variant_start, allele_start, allele_stop)    # ooh scandalous

            prev_start = allele_start

        return labeled_sites

    def print_site_data(self, candidate_site, alt1_sequence, alt2_sequence):
        allele_start = candidate_site["window_start"]
        allele_stop = candidate_site["window_end"]
        ref_sequence = candidate_site["reference_sequence"]

        print("\nCandidate: ")
        print("window: ", allele_start, allele_stop)
        print("ref:    ", ref_sequence)
        print("alts:   ", alt1_sequence, alt2_sequence)
        for record in candidate_site["vcf_records"]:
            for variant in record:
                print("\nVariant:  ")
                print("position  ", variant.pos-1)
                print("ref:      ", variant.ref)
                print("alt:      ", variant.alt)
                print("genotype: ", variant.genotype_class)
                print("type      ", variant.type)

    def update_candidate_site(self, allele, site_variants, alt1_sequence, alt2_sequence):
        if "vcf_records" not in allele:
            allele["vcf_records"] = list()
            is_multi_vcf_site = False
        else:
            is_multi_vcf_site = True

        if site_variants is not None:
            genotype = site_variants[0].type            # use the first genotype, assuming it is the same as any others
            allele["vcf_records"].append(site_variants) # each entry pertains to 1 or 2 records from a given coordinate
        else:
            genotype = "Hom"

        site = dict()
        site["chromosome_name"] = allele["chromosome_name"]
        site["window_start"] = allele["window_start"]
        site["window_end"] = allele["window_end"]
        site["ref_sequence"] = allele["reference_sequence"]
        site["alt1_sequence"] = alt1_sequence
        site["alt2_sequence"] = alt2_sequence
        site["genotype"] = genotype

        return allele, site, is_multi_vcf_site

    def get_merged_variants_for_site(self, variant, window_start, window_stop):
        '''

        unfinished testing method ... do the vcf sequences actually match the candidate site alleles?

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
        '''
        unfinished testing method
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
#
# if __name__ == '__main__':
#     '''
#     Processes arguments and performs tasks to generate the pileup.
#     '''
#     parser = argparse.ArgumentParser()
#     parser.register("type", "bool", lambda v: v.lower() == "true")
#     parser.add_argument(
#         "--ref",
#         type=str,
#         required=True,
#         help="Reference corresponding to the BAM file."
#     )
#     parser.add_argument(
#         "--chromosome_name",
#         type=str,
#         default="3",
#         help="Desired chromosome number E.g.: 3"
#     )
#     parser.add_argument(
#         "--test_distributions",
#         type=bool,
#         default=False,
#         help="If true then a dry test is run."
#     )
#     parser.add_argument(
#         "--json_dir",
#         type=str,
#         help="Directory containing all json allele files to test."
#     )
#
#     FLAGS, unparsed = parser.parse_known_args()
#     # process the output directory
#     if FLAGS.output_dir[-1] != '/':
#         FLAGS.output_dir += '/'
#     if not exists(FLAGS.output_dir):
#         mkdir(FLAGS.output_dir)
#
#     view = View(chromosome_name=FLAGS.chromosome_name,
#                 bam_file_path=FLAGS.bam,
#                 reference_file_path=FLAGS.ref,
#                 output_file_path=FLAGS.output_dir)
#
#     if FLAGS.test is True:
#         view = View(chromosome_name=FLAGS.chromosome_name,
#                     bam_file_path=FLAGS.bam,
#                     reference_file_path=FLAGS.ref,
#                     output_file_path=FLAGS.output_dir)
#         view.test(FLAGS.json)
#     else:
#         do_parallel(FLAGS.chromosome_name, FLAGS.bam, FLAGS.ref, FLAGS.json, FLAGS.output_dir, FLAGS.max_threads)
