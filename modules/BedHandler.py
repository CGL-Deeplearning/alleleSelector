from pybedtools import BedTool

"""
This class handles bed files using pybedtools.
"""


class BedHandler:
    """
    Handles bed files using pybedtools API
    """
    def __init__(self, bed_file_path):
        """
        create AlignmentFile object given file path to a bam file
        :param bam_file_path: full path to a bam file
        """
        self.bed_file_path = bed_file_path
        try:
            self.bed_file = BedTool(self.bed_file_path)
        except:
            raise IOError("BED FILE READ ERROR")

    @staticmethod
    def list_to_bed(list_bed_format):
        """
        Create a bedtools object from a list in bedtools format
        :param list_bed_format: List
        :return: BedTools object
        """
        bedtool_obj = BedTool(list_bed_format)
        return bedtool_obj

    def intersect(self, bed_object):
        """
        Intersect with another BedHandler instance
        :return:
        """
        return self.bed_file.intersect(bed_object.bed_file)

    def __getitem__(self, index):
        return self.bed_file[index]

    def __len__(self):
        return len(self.bed_file)