from modules.VcfHandler import VCFFileProcessor
from modules.VcfHandlerNonPositional import VCFFileProcessor as VCFFileProcessorNP
from modules.FastaHandler import FastaHandler
from modules.BedHandler import BedHandler


fasta_path = "/Users/saureous/data/ref_hg19_WG.fa"
prefix = "chr"

chromosome_lengths = [249250621,
                      243199373,
                      198022430,
                      191154276,
                      180915260,
                      171115067,
                      159138663,
                      146364022,
                      141213431,
                      135534747,
                      135006516,
                      133851895,
                      115169878,
                      107349540,
                      102531392,
                      90354753,
                      81195210,
                      78077248,
                      59128983,
                      63025520,
                      48129895,
                      51304566]

vcf_path = "/users/saureous/data/NA12878_S1_confident.genome.vcf.gz"

for i in range(21,22):
    chr_name = prefix+str(i+1)

    print(chr_name)

    vcf_handler = VCFFileProcessor(file_path=vcf_path)
    vcf_handler.populate_dictionary(contig=chr_name, start_pos=1, end_pos=chromosome_lengths[i], hom_filter=True)

    chromosome_1_vcf = vcf_handler.get_variant_dictionary()

    print("positional: ", len(chromosome_1_vcf))

    # vcf_handler.save_positional_vcf_as_bed(chromosome_name="chr1", output_path_name="/users/saureous/data/vcf_chr1_full.bed")

    site = ":1-" + str(chromosome_lengths[i])

    vcf_handler_np = VCFFileProcessorNP(file_path=vcf_path)
    vcf_handler_np.populate_dictionary(contig=chr_name, site=site, hom_filter=True)

    chromosome_1_vcf_np = vcf_handler_np.get_variant_dictionary()

    print("condensed:  ", len(chromosome_1_vcf_np))


exit()

for i in range(1, 23):
    chr_name = prefix+str(i)
    fasta_handler = FastaHandler(fasta_path)
    print(fasta_handler.get_chr_sequence_length(chromosome_name=chr_name))


