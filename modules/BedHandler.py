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

        if bed_file_path is not None:
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

    def intersect(self, bed_object, u=False):
        """
        Intersect with another BedHandler instance, return BedHandler wrapper for new BED
        :return: BedHandler object for intersected BED
        """
        out_bed_object = BedHandler(None)
        out_bed_object.bed_file = self.bed_file.intersect(bed_object.bed_file,u=True)

        return out_bed_object

    def save(self, output_file_path):
        """
        Write the BEDTools object to hard drive
        :param output_file_path:
        :return:
        """
        self.bed_file.saveas(output_file_path)

    def __getitem__(self, index):
        return self.bed_file[index]

    def __len__(self):
        return len(self.bed_file)
