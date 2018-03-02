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
            raise IOError("BAM FILE ERROR")

    def get_pileupcolumns_aligned_to_a_site(self, contig, pos):
        """
        Return a AlignmentFile.pileup object given a site
        :param contig: Contig [ex. chr3]
        :param pos: Position [ex 100001]
        :return: pysam.AlignmentFile.pileup object
        """
        # get pileup columns
        pileup_columns = self.bamFile.pileup(contig, pos, pos+1)
        # return pileup columns
        return pileup_columns

    def get_pileup_of_a_site(self, contig, pos):
        """
        ***USED FOR TESTING***

        This method returns a string containing all bases of a site.
        :param contig: Contig [ex. chr3]
        :param pos: Position [ex 100001]
        :return: String containing the pileup of the site
        """
        # initialize pileup
        pileup_str = ""
        for pileupcolumn in self.bamFile.pileup(contig, pos, pos+1, truncate=True):
            if pileupcolumn.pos < pos:
                continue
            elif pileupcolumn.pos >= pos+1:
                break
            # the position
            pileup_str += str(pileupcolumn.pos) + " "
            # each read in the pileup
            for pileupread in pileupcolumn.pileups:
                # query position is None if is_del or is_refskip is set.
                if not pileupread.is_del and not pileupread.is_refskip:
                    pileup_str += pileupread.alignment.query_sequence[pileupread.query_position]
                if pileupread.is_del:
                    pileup_str += '*'

        # return pileup string
        return pileup_str

