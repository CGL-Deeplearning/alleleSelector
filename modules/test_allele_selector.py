from AlleleFileHandler import AlleleFileHandler
from AlleleSelector import AlleleSelector
from vcf_handler import VCFFileProcessor
from FastaHandler import FastaHandler

def get_labeled_sites(vcf_file_path, fasta_file_path, allele_file_path):
    fasta_handler = FastaHandler(reference_file_path=fasta_file_path)
    vcf_handler = VCFFileProcessor(file_path=vcf_file_path)
    allele_handler = AlleleFileHandler(allele_file_path=allele_file_path)

    candidate_dictionary = allele_handler.get_allele_site_dictionary()
    candidate_sites = candidate_dictionary["all_candidates"]

    allele_selector = AlleleSelector(fasta_handler=fasta_handler)
    labeled_sites = allele_selector.get_labeled_candidate_sites(vcf_handler=vcf_handler, candidate_sites=candidate_sites, filter_hom_ref=True)


fasta_file_path = "/Users/saureous/data/chr3.fa"
vcf_file_path = "/Users/saureous/data/NA12878_S1_confident.genome.vcf.gz"
allele_file_path = "/Users/saureous/data/json/Candidates_chr3_118813464_125415212.json"

get_labeled_sites(vcf_file_path, fasta_file_path, allele_file_path)