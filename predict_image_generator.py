import sys
import time
import random
import multiprocessing

from image_analyzer import *
from modules.KalleR.Bed2Image_API import Bed2ImageAPI
from modules.TextColor import TextColor
from modules.FileManager import FileManager
from image_analyzer import *

CLASS_BY_INDEX = ["HOM", "HET", "HOM_ALT"]
class_count = {}


def generate_pileup(contig, bam_file, ref_file, records, output_dir, thread_name):
    """
    Generate pileup images from a vcf file
    :param contig: Which contig to fetch ("chr3")
    :param site: Which site to fetch (":100000-200000")
    :param bam_file: Path to the bam alignment file
    :param ref_file: Path to the reference file
    :param vcf_file: Path to the vcf file
    :param output_dir: Output directory, where the image will be saved
    :return:
    """
    # create a summary file
    smry = open(output_dir + "summary/" + "summary" + '_' + contig + "_" + thread_name + ".csv", 'w')
    st_time = time.time()
    for rec in records:
        rec_id, chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type = tuple(rec)
        if alt2 == '.':
            alt = alt1
        else:
            alt = alt1 + "_" + alt2

        # file name for the image and save the image
        file_name = str(rec_id) + "_" + chr_name + "_" + str(pos_start) + "_" + str(ref) + "_" + str(alt)

        api_object = Bed2ImageAPI(bam_file, ref_file)
        img, img_shape = api_object.create_image(api_object.bam_handler, api_object.fasta_handler, rec, output_dir, file_name)

        # label of the image and save the image
        rec_str = ' '.join(str(e) for e in rec)

        smry.write(os.path.abspath(output_dir + file_name) + ".png," + ','.join(
            map(str, img_shape)) + ',' + rec_str + ',' + str(rec_id) + '\n')

    # sys.stderr.write(TextColor.PURPLE + "FINISHED: " + thread_name + " TIME: " + str(time.time()-st_time) + "\n" + TextColor.END)


def get_images_for_two_alts(record, rec_id):
    chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type = record.rstrip().split('\t')[0:7]

    rec_1 = [rec_id, chr_name, pos_start, pos_end, ref, alt1, '.', rec_type]
    rec_2 = [rec_id, chr_name, pos_start, pos_end, ref, alt2, '.', rec_type]
    rec_3 = [rec_id, chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type]

    return [rec_1, rec_2, rec_3]


def get_prediction_set_from_bed(candidate_bed):
    with open(candidate_bed) as bed_file:
        bed_records = bed_file.readlines()

    train_set = {}

    rec_id = 1
    for record in bed_records:

        chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type = record.rstrip().split('\t')[0:7]
        if chr_name not in train_set.keys():
            train_set[chr_name] = []

        if alt2 != '.':
            train_set[chr_name].extend(get_images_for_two_alts(record, rec_id))
        else:
            train_set[chr_name].append([rec_id, chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type])
        rec_id += 1

    return train_set


def create_output_dir_for_chromosome(output_dir, chr_name):
    """
    Create an internal directory inside the output directory to dump choromosomal bed files
    :param output_dir: Path to output directory
    :param chr_name: chromosome name
    :return: New directory path
    """
    path_to_dir = output_dir + chr_name + "/"
    if not os.path.exists(path_to_dir):
        os.mkdir(path_to_dir)

    summary_path = path_to_dir + "summary" + "/"
    if not os.path.exists(summary_path):
        os.mkdir(summary_path)

    return path_to_dir


def chromosome_level_parallelization(chr_name, bam_file, ref_file, records, output_dir, max_threads):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param chr_name: Chromosome name
    :param bam_file: Bam file
    :param ref_file: Ref file
    :param vcf_file: VCF file
    :param output_dir: Output directory
    :param max_threads: Maximum number of threads
    :return: A list of results returned by the processes
    """
    chunks = int(math.ceil(len(records) / max_threads))

    for i in range(max_threads):
        # parse window of the segment. Use a 1000 overlap for corner cases.
        start_position = i * chunks
        end_position = min((i + 1) * chunks, len(records))
        subset_records = records[start_position:end_position]
        args = (chr_name, bam_file, ref_file, subset_records, output_dir, str(i))

        p = multiprocessing.Process(target=generate_pileup, args=args)
        p.start()

        while True:
            if len(multiprocessing.active_children()) < max_threads:
                break


def generate_images(bam_file, ref_file, candidate_bed, output_dir, max_threads):
    train_set = get_prediction_set_from_bed(candidate_bed)
    program_start_time = time.time()

    chr_list = train_set.keys()
    for chr_name in train_set:
        sys.stderr.write(TextColor.BLUE + "STARTING " + str(chr_name) + " PROCESSES" + "\n")
        sys.stderr.write(TextColor.BLUE + "TOTAL " + str(len(train_set[chr_name])) + " RECORDS" + "\n")

        start_time = time.time()

        # create dump directory inside output directory
        chr_output_dir = create_output_dir_for_chromosome(output_dir, chr_name)

        # do a chromosome level parallelization
        chromosome_level_parallelization(chr_name, bam_file, ref_file, train_set[chr_name], chr_output_dir, max_threads)

        end_time = time.time()
        sys.stderr.write(TextColor.PURPLE + "FINISHED " + str(chr_name) + " PROCESSES" + "\n")
        sys.stderr.write(TextColor.CYAN + "TIME ELAPSED: " + str(end_time - start_time) + "\n" + TextColor.END)

        # wait for the last process to end before file processing
    while True:
        if len(multiprocessing.active_children()) == 0:
            break

    for chr in chr_list:

        # here we dumped all the bed files
        path_to_dir = output_dir + chr + '/' + 'summary' + "/"

        concatenated_file_name = output_dir + chr + "_train.csv"

        filemanager_object = FileManager()
        # get all bed file paths from the directory
        file_paths = filemanager_object.get_file_paths_from_directory(path_to_dir)
        # dump all bed files into one
        filemanager_object.concatenate_files(file_paths, concatenated_file_name)
        # delete all temporary files
        filemanager_object.delete_files(file_paths)
        os.rmdir(path_to_dir)

    program_end_time = time.time()
    sys.stderr.write(TextColor.RED + "PROCESSED FINISHED SUCCESSFULLY" + "\n")
    sys.stderr.write(
        TextColor.CYAN + "TOTAL TIME FOR GENERATING ALL RESULTS: " + str(program_end_time - program_start_time) + "\n")


def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # create an internal directory so we don't overwrite previous runs
    timestr = time.strftime("%m%d%Y_%H%M%S")
    internal_directory = "run_" + timestr + "/"
    output_dir = output_dir + internal_directory
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return output_dir


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file containing reads of interest."
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--candidate_bed",
        type=str,
        required=True,
        help="Bed file containing all candidates"
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="predict_image_output/",
        help="Path to output directory."
    )
    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = handle_output_directory(FLAGS.output_dir)

    generate_images(FLAGS.bam, FLAGS.ref, FLAGS.candidate_bed, FLAGS.output_dir, FLAGS.max_threads)