from VcfHandler import VCFFileProcessor
from BedHandler import BedHandler
import pybedtools
import argparse

"""
EXAMPLE USAGE:

python3 modules/FilterTest.py --vcf /users/saureous/data/NA12878_S1.genome.vcf.gz --allele_bed /users/saureous/data/bed_alleles/Label_chr3_120000_140010.bed --confident_bed /users/saureous/data/ConfidentRegions.bed.gz

"""

# VCF record indexes
REF,ALT,GT = 0,1,2

# BED record indexes
CHROMOSOME_NAME,START,STOP,REF_BED,ALT_BED,GT_BED = 0,1,2,3,4,5

class View:
    """
    I/O manager for FilterTest
    """
    def __init__(self, vcf_file_path, confident_bed_file_path, allele_bed_file_path):
        self.vcf_handler = VCFFileProcessor(vcf_file_path)
        self.bed_handler_confident = BedHandler(confident_bed_file_path)
        self.bed_handler_allele = BedHandler(allele_bed_file_path)
        self.tester = FilterTest()

    def get_allele_validation_stats(self, deallocate=False):
        confident_alleles = self.bed_handler_allele.intersect(self.bed_handler_confident)

        if deallocate:
            del self.bed_handler_allele
            del self.bed_handler_confident

        length = len(confident_alleles)
        chromosome_name = confident_alleles[0].chrom
        allele_start = confident_alleles[0].start
        allele_stop = confident_alleles[length-1].stop

        print("Intersection region:       ", chromosome_name, allele_start, allele_stop)
        print("Items before intersection: ",len(self.bed_handler_allele))
        print("Items after intersection:  ",len(confident_alleles))

        # get dictionary of variant records for full region
        self.vcf_handler.populate_dictionary(contig=chromosome_name,
                                             start_pos=allele_start,
                                             end_pos=allele_stop,
                                             hom_filter=True)

        # get separate positional variant dictionaries for IN, DEL, and SNP
        positional_variants = self.vcf_handler.get_variant_dictionary()

        self.tester.validate_alleles(positional_variants, confident_alleles)


class FilterTest:
    """
    Performs summary statistics on B
    """

    def validate_alleles(self, variants, alleles):
        """

        :param variants:
        :return:
        """
        validated_vcf_positions = {key: 0 for key in variants.keys()}
        n_false_negative = 0
        n_false_positive = 0
        n_true_positive = 0

        for entry in alleles:
            genotype = int(entry[GT_BED])
            allele_start = int(entry[START])

            # test if the site has any true variants (not Hom)
            if self._is_supported(genotype):
                if allele_start in validated_vcf_positions:
                    validated_vcf_positions[allele_start] = 1
                    n_true_positive += 1
            else:
                n_false_positive += 1

        for pos in validated_vcf_positions:
            if validated_vcf_positions[pos] == 0:
                n_false_negative += 1
                print("\nWARNING: Unsupported VCF position: ", pos)
                print("\tRecord: ", variants[pos])

        print("False Negative: ", n_false_negative)
        print("False Positive: ", n_false_positive)
        print("True Positive:  ", n_true_positive)

    def _is_supported(self,genotype):
        supported = False

        if genotype > 0:
            supported = True

        return supported


if __name__ == '__main__':
    """
    Processes arguments and performs tasks to generate the pileup.
    """
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF file containing all known variant records."
    )
    parser.add_argument(
        "--confident_bed",
        type=str,
        required=True,
        help="VCF file containing all known variant records."
    )
    parser.add_argument(
        "--allele_bed",
        type=str,
        required=True,
        help="VCF file containing all known variant records."
    )

    FLAGS, unparsed = parser.parse_known_args()

    view = View(vcf_file_path=FLAGS.vcf,
                confident_bed_file_path=FLAGS.confident_bed,
                allele_bed_file_path=FLAGS.allele_bed)

    view.get_allele_validation_stats()
