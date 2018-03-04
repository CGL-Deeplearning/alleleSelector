from modules.KalleR.ImageCreator import ImageCreator
from modules.KalleR.BamHandler import BamHandler
from modules.KalleR.FastaHandler import FastaHandler
import os

class TrainBed2ImageAPI:
    """
    Works as a main class and handles user interaction with different modules.
    """
    def __init__(self, bam_file_path, reference_file_path):
        # --- initialize handlers ---
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)

    @staticmethod
    def create_image(bam_handler, fasta_handler, record, output_dir, file_name):
        """
        Create an image from a bed record
        :param bam_handler: Handles bam file
        :param fasta_handler: Handles fasta file
        :param bed_record: Bed record
        :return: Imagearray, label
        """
        chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type, label = tuple(record)

        alts = [alt1]
        if alt2 != '.':
            alts.append(alt2)

        start_position = int(pos_start)

        pileups = bam_handler.get_pileupcolumns_aligned_to_a_site(chr_name, start_position)
        image_creator = ImageCreator(fasta_handler, pileups, chr_name, start_position, alts)

        image_array, image_shape = image_creator.create_image(start_position, ref, alts)
        image_creator.save_image_as_png(image_array, output_dir, file_name)

        return image_array, image_shape