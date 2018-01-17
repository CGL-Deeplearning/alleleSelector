"""
Implemented by: Kishwar Shafin
Date: 02/01/2018
"""
import pysam

"""
This class handles bam files using pysam API. Two methods mostly handles the BAM interaction:
- __init__: create AlignmentFile object given file path to a bam file
- get_pileup_columns_aligned_to_a_site: return a AlignmentFile.pileup object given a site
"""


class BamProcessor:
    """
    Handles bam files using pysam API
    """
    def __init__(self, file_path):
        """
        create AlignmentFile object given file path to a bam file
        :param file_path: file_path to a bam file
        """
        self.file_path = file_path
        try:
            self.bamFile = pysam.AlignmentFile(self.file_path, "rb")
        except:
            raise IOError("BAM FILE READ ERROR")

    def get_reads_from_region(self, contig, start_pos, end_pos):
        return self.bamFile.fetch(contig, start_pos, end_pos)

