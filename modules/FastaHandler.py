from pyfaidx import Fasta

"""
This class handles fasta reference files, ensuring that the sequence is not a terminal 'N' and that the end of the
sequence has not been reached
"""


class FastaHandler:
    """
    Handles fasta files using pyfaidx API
    """
    def __init__(self, reference_file_path):
        """
        create fasta file object given file path to a fasta reference file
        :param fasta_file_path: full path to a fasta reference file
        """

        self.fasta_file_path = reference_file_path

        try:
            self.fasta = Fasta(self.fasta_file_path, as_raw=True, sequence_always_upper=True)
        except:
            raise IOError("BAM FILE READ ERROR")

    def get_sequence(self, chromosome_name, start, stop):
        return self.fasta.get_seq(name=chromosome_name, start=start+1, end=stop+1)  # must be 1-based

