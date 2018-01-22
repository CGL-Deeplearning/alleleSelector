from AlleleFileHandler import AlleleFileHandler
from AlleleSelector import AlleleSelector
from vcf_handler import VCFFileProcessor
from FastaHandler import FastaHandler

fasta_file_path = "/Users/saureous/data/chr3.fa"
vcf_file_path = "/Users/saureous/data/NA12878_S1.genome.vcf.gz"
allele_file_path = "/Users/saureous/data/Candidates_chr3_100000_200000.json.txt"

fasta = FastaHandler(reference_file_path=fasta_file_path)
vcf = VCFFileProcessor(file_path=vcf_file_path)
alleles = AlleleFileHandler(allele_file_path=allele_file_path)

selector = AlleleSelector(allele_file_handler=alleles,vcf_handler=vcf,fasta_handler=fasta)

alleles = selector.parse_candidate_sites()

selector.test_alleles(alleles)
