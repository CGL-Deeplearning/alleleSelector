from modules.AlleleFileHandler import AlleleFileHandler
from modules.CandidateLabeler import CandidateLabeler
from modules.vcf_handler import VCFFileProcessor
from modules.FastaHandler import FastaHandler
from os import listdir, mkdir
from os.path import isfile, join, exists
import argparse
import json


class ComplexEncoder(json.JSONEncoder):
    """
    JSON encoder for class attributes
    """
    def default(self, obj):
        if hasattr(obj, 'reprJSON'):
            return obj.reprJSON()
        else:
            return json.JSONEncoder.default(self, obj)


class View:
    '''
    Main class for candidate labeler module, can parse alleles from json file if run from command line, or directly
    using get_labeled_candidate_sites method.
    '''
    def __init__(self, fasta_file_path, vcf_file_path):
        self.fasta_handler = FastaHandler(reference_file_path=fasta_file_path)
        self.vcf_handler = VCFFileProcessor(file_path=vcf_file_path)

        self.chromosome_name = None
        self.start = None
        self.stop = None

    def get_labeled_candidate_sites(self, candidate_dictionary, filter_hom_ref=True):
        '''
        Takes a dictionary of allele data and compares with a VCF to determine which candidate alleles are supported.
        :param candidate_dictionary: dictionary with the list of allele sites under the top-level entry "all_candidates"
        :param filter_hom_ref: whether to ignore hom_ref VCF records during candidate validation
        :return: labeled_sites: the parsed candidate dictionary with the following structure:

        {window_start:
            {
            chromosome_name
            window_start
            window_end
            ref_sequence
            alt1_sequence
            alt2_sequence
            genotype
            }
        }
        '''
        candidate_sites = candidate_dictionary["all_candidates"]

        # find start and stop position of region that covers all candidates
        start_record = candidate_sites[0]
        stop_record = candidate_sites[-1]

        self.chromosome_name = start_record["chromosome_name"]
        self.start = start_record["window_start"]
        self.stop = stop_record["window_end"]

        # get dictionary of variant records for full region
        self.vcf_handler.populate_dictionary(contig=self.chromosome_name,
                                             start_pos=self.start,
                                             end_pos=self.stop,
                                             hom_filter=filter_hom_ref)

        variants = self.vcf_handler.get_variant_dictionary()

        allele_selector = CandidateLabeler(fasta_handler=self.fasta_handler)

        labeled_sites = allele_selector.get_labeled_candidates(variants=variants, candidate_sites=candidate_sites)

        return labeled_sites

    def geenerate_labeled_sites(self, start_position, end_position, json_out):
        pass

    def test_json_directory(self, json_directory_path):
        '''
        :param json_directory_path: path containing any number of json files with candidate site data
        :return:
        '''
        json_directory_path = json_directory_path
        json_file_paths = [join(json_directory_path, file) for file in listdir(json_directory_path) if
                           isfile(join(json_directory_path, file))]

        for file_path in json_file_paths:
            allele_handler = AlleleFileHandler(allele_file_path=file_path)
            candidate_dictionary = allele_handler.get_allele_site_dictionary()
            labeled_sites = self.get_labeled_candidate_sites(candidate_dictionary=candidate_dictionary,
                                                             filter_hom_ref=True)

            self.write_json(candidate_dict=labeled_sites)

    def write_json(self, candidate_dict, output_dir="output/labeled_candidates/"):
        """
        Create a json output of all candidates found in the region
        :param candidate_dict: Candidate dict to be saved
        :param output_dir: where to save the labeled candidates json file
        :return:
        """
        if not exists(output_dir):
            mkdir(output_dir)
        with open(output_dir + "labeled_candidates" + '_' + self.chromosome_name + '_' + str(self.start) + '_' + str(self.stop) + ".json", 'w') as json_file:
            json.dump(candidate_dict, json_file, indent=4, sort_keys=True)

if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference corresponding to the VCF file."
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF file containing all known variant records."
    )
    parser.add_argument(
        "--json_dir",
        type=str,
        default=None,
        required=True,
        help="Specify a directory containing any number of json allele files to parse"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/",
        help="If true then output will be in a json file in json folder."
    )

    FLAGS, unparsed = parser.parse_known_args()

    view = View(fasta_file_path=FLAGS.ref,
                vcf_file_path=FLAGS.vcf)

    view.test_json_directory(FLAGS.json_dir)
