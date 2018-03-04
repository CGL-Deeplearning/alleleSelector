import sys
import time
import random
import multiprocessing

from image_analyzer import *
from modules.KalleR.TrainBed2Image_API import TrainBed2ImageAPI
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
        chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type, label = tuple(rec)
        if alt2 == '.':
            alt = alt1
        else:
            alt = alt1 + "_" + alt2

        # file name for the image and save the image
        file_name = chr_name + "_" + str(pos_start) + "_" + str(ref) + "_" + str(alt) + "_" + str(label)

        api_object = TrainBed2ImageAPI(bam_file, ref_file)
        img, img_shape = api_object.create_image(api_object.bam_handler, api_object.fasta_handler, rec, output_dir, file_name)

        # label of the image and save the image
        smry.write(os.path.abspath(output_dir + file_name) + ".png," + str(label) + ',' + ','.join(
            map(str, img_shape)) + ',' + str(rec_type) + '\n')

    # sys.stderr.write(TextColor.PURPLE + "FINISHED: " + thread_name + " TIME: " + str(time.time()-st_time) + "\n" + TextColor.END)


def get_combined_gt(gt1, gt2):
    if gt1 == '0':
        return gt2
    if gt2 == '0':
        return gt1
    if gt1 == '0' and gt2 == '0':
        return '0'
    if gt1 == '1' and gt2 == '1':
        return '2'
    return None


def _initialize_class_count_dictionary(chr_name):
    if chr_name not in class_count.keys():
        class_count[chr_name] = {'0': 0, '1': 0, '2': 0}


def get_images_for_two_alts(record):
    chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type = record.rstrip().split('\t')[0:7]
    gt1, filter1, gt2, filter2 = record.rstrip().split('\t')[7:11]
    _initialize_class_count_dictionary(chr_name)

    class_count[chr_name][gt1] += 1
    class_count[chr_name][gt2] += 1
    gt3 = get_combined_gt(gt1, gt2)
    if gt3 is None:
        sys.stderr.write(TextColor.RED + "WEIRD RECORD: " + str(record) + "\n")
    rec_1 = [chr_name, pos_start, pos_end, ref, alt1, '.', rec_type, gt1]
    rec_2 = [chr_name, pos_start, pos_end, ref, alt2, '.', rec_type, gt2]
    if gt3 is not None:
        class_count[chr_name][str(gt3)] += 1
        rec_3 = [chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type, gt3]
        return [rec_1, rec_2, rec_3]

    return [rec_1, rec_2]


def select_or_not(downsample_rate):
    """
    Determines if a bed record should be selected given a downsampling rate
    :param bed_record: A bed record
    :param downsample_rate: A downsampling probability
    :return: Boolean
    """
    # else do a sampling based on probability
    random_chance = random.uniform(0, 1)
    if random_chance <= downsample_rate:
        return True
    return False


def get_downsample_rate(total_hom, total_het, total_hom_alt):
    """
    Downsample the bed file
    :return:
    """
    # calculate the downsample rate based on distribution of three classes
    downsample_rate = max(total_het, total_hom_alt) / total_hom
    # we want the homozygous to be twice the size of the next most frequent class.
    downsample_rate = 2 * downsample_rate

    return downsample_rate


def get_prediction_set_from_bed(candidate_bed):
    with open(candidate_bed) as bed_file:
        bed_records = bed_file.readlines()

    train_set = {}

    rec_id = 1
    for record in bed_records:

        chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type = record.rstrip().split('\t')[0:7]
        gt1, filter1, gt2, filter2 = record.rstrip().split('\t')[7:11]
        if chr_name not in train_set.keys():
            train_set[chr_name] = []
        _initialize_class_count_dictionary(chr_name)

        if alt2 != '.':
            train_set[chr_name].extend(get_images_for_two_alts(record))
        else:
            train_set[chr_name].append([chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type, gt1])
            class_count[chr_name][gt1] += 1
        rec_id += 1

    for chr in train_set:
        downsample_rate = get_downsample_rate(class_count[chr]['0'], class_count[chr]['1'], class_count[chr]['2'])
        if downsample_rate >= 1:
            continue
        downsampled_list = []
        for prediction in train_set[chr]:
            alt2 = prediction[5]
            gt = prediction[7]
            if gt == '0' and alt2 == '.' and select_or_not(downsample_rate) is False:
                continue
            downsampled_list.append(prediction)
        train_set[chr] = downsampled_list

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
        default="train_image_output/",
        help="Path to output directory."
    )
    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = handle_output_directory(FLAGS.output_dir)

    generate_images(FLAGS.bam, FLAGS.ref, FLAGS.candidate_bed, FLAGS.output_dir, FLAGS.max_threads)