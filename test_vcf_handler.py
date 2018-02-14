from modules.VcfHandler import VCFFileProcessor


vcf_path = "/Users/saureous/data/NA12878_S1.genome.vcf.gz"
contig = "chr1"
start = 0
stop = 1000000

print()
print(start,stop)
print()

handler = VCFFileProcessor(file_path=vcf_path)

handler.populate_dictionary(contig=contig, start_pos=start, end_pos=stop, hom_filter=True)

variants = handler.get_variant_dictionary()

for pos in variants:
    print(variants[pos])

#
# start = 0
# stop = 200000
#
# print()
# print(start,stop)
# print()
#
# handler = VCFFileProcessor(file_path=vcf_path)
#
# handler.populate_dictionary(contig=contig,start_pos=start,end_pos=stop,hom_filter=True)
#
# variants = handler.get_variant_dictionary()
#
# for pos in variants:
#     print(variants[pos])
