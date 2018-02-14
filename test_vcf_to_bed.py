from modules.VcfHandler import VCFFileProcessor

vcf_path = "/users/saureous/data/NA12878_S1_confident.genome.vcf.gz"
output_path_name = "/Users/saureous/data/vcf_test.bed"
output_path_name_2 = "/Users/saureous/data/vcf_test2.bed"
chromosome_name = "chr3"

vcf_handler = VCFFileProcessor(vcf_path)

vcf_handler.populate_dictionary(contig=chromosome_name, start_pos=100000, end_pos=200000, hom_filter=False)

vcf_handler.save_positional_vcf_as_bed(output_path_name=output_path_name, chromosome_name=chromosome_name)

positional_vcf = vcf_handler.read_positional_vcf_from_bed(output_path_name)

vcf_handler.save_positional_vcf_as_bed(output_path_name=output_path_name_2, chromosome_name=chromosome_name)