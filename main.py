import pysam
import modules.bam_handler as bam_processor

bam_object = bam_processor.BamProcessor("/Users/kishwar/Kishwar/Whole_chr3_data/illumina/vcf_whole_chr/chr3.bam")
reads = bam_object.get_reads_from_region("chr3", 100000, 100100)
for read in reads:
    print(read)
