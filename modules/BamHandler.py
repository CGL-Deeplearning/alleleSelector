import pysam

"""
This class handles bam files using pysam API.
"""


class BamHandler:
    """
    Handles bam files using pysam API
    """
    def __init__(self, bam_file_path):
        """
        create AlignmentFile object given file path to a bam file
        :param bam_file_path: full path to a bam file
        """
        self.bam_file_path = bam_file_path
        try:
            self.bamFile = pysam.AlignmentFile(self.bam_file_path, "rb")
        except:
            raise IOError("BAM FILE READ ERROR")

    def get_reads(self, chromosome_name, start, stop):
        """
        Return reads that map to a given site
        :param chromosome_name: Chromosome name. Ex: chr3
        :param start: Site start in the chromosome
        :param stop: Site end in the chromosome
        :return: Reads that align to that site
        """
        return self.bamFile.fetch(contig=chromosome_name, start=start, end=stop)

