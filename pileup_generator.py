"""
Implemented by: Kishwar Shafin
Date: 02/01/2018
"""
import time
import sys
import random
import multiprocessing
from image_analyzer import *
from modules.TextColor import TextColor
from modules.genomePileup.vcf_handler import VCFFileProcessor
from modules.genomePileup.ref_handler import RefFileProcessor
from modules.genomePileup.bam_handler_mpileup import BamProcessor
from modules.genomePileup.pileup_creator import PileupProcessor
from modules.FileManager import FileManager

TEST_CHR = "chr3"
TEST_REGION = "100000-200000"


def handle_directory(directory_path):
    """
    Create a directory if doesn't exist
    :param directory_path: path to the directory
    :return: desired directory name
    """
    # if directory has no trailing '/' then add it
    if directory_path[-1] != '/':
        directory_path += '/'
    # if directory doesn't exist then create it
    if not os.path.exists(directory_path):
        os.mkdir(directory_path)

    return directory_path


def get_label(genotype_type):
    """
    Get genotype label for a type of genotype
    :param genotype_type: Genotype in string
    :return: Integer label for the type
    """
    if genotype_type == "Hom":
        return 0
    elif genotype_type == "Het":
        return 1
    elif genotype_type == "Hom_alt":
        return 2


def chunk_it(seq, num):
    """
    Chunk a sequence in N equal segments
    :param seq: Sequence of numbers
    :param num: Number of chunks
    :return: chunked start and end positions
    """
    # find average chunk size
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    # until the end of sequence
    while last < len(seq):
        # append the value to a bin
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


def get_alts_in_hom_pileup(pileup_str, ref_base):
    """
    Return possible alts in homozygous cases.
    :param pileup_str: Pileup of a homozygous base.
    :param ref_base: Reference base of that position.
    :return:
    """
    alts = {'A':0, 'C':0, 'G':0, 'T':0}
    for base in pileup_str:
        if base != ref_base and base in alts.keys():
            alts[base] += 1

    return max(alts, key=alts.get), alts[max(alts, key=alts.get)]


def get_odds_for_hom(total_hom, total_het, total_homalt):
    """
    This class will return the odds of generating an image each time we see a hom case
    :param total_hom: Total hom cases present
    :param total_het: Total het cases present
    :param total_homalt: Total homalt cases present
    :return:
    """
    probability_of_seeing_hom = total_hom / (total_hom + total_het + total_homalt)
    odds_of_selecting_hom = 1.0 - probability_of_seeing_hom

    return odds_of_selecting_hom


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


def test(contig, site, bam_file, ref_file, vcf_file, output_dir):
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

    # create the vcf handler
    vcf_handler = VCFFileProcessor(vcf_file)
    # generate dictionary of the region
    vcf_handler.populate_dictionary(contig, site, hom_filter=False)

    # create ref and bam files handler
    ref_handler = RefFileProcessor(ref_file)
    bam_handler = BamProcessor(bam_file)

    # create a summary file
    smry = open(output_dir + "summary" + '_' + contig + site.replace(':', '_').replace('-', '_') + ".csv", 'w')

    # get the vcf dictionary of that region
    vcf_dict = vcf_handler.get_variant_dictionary()

    # get the odds of selecting a homozygous case
    total_hom, total_het, total_homalt = vcf_handler.get_genotype_counts()
    downsample_rate = get_downsample_rate(total_hom, total_het, total_homalt)

    # keep count of how many images of each type is generated
    total_generated_hom, total_generated_het, total_generated_hom_alt = 0, 0, 0

    start_time = time.time()
    for pos in vcf_dict.keys():
        for rec in vcf_dict[pos]:
            alt = '.'
            if rec.type == 'Hom' and select_or_not(downsample_rate) is False:
                continue
            elif rec.type == 'Hom':
                rec.alt = alt

            total_generated_hom += 1 if rec.type == 'Hom' else 0
            total_generated_het += 1 if rec.type == 'Het' else 0
            total_generated_hom_alt += 1 if rec.type == 'Hom_alt' else 0

            # get pileup columns from bam file
            pileup_columns = bam_handler.get_pileupcolumns_aligned_to_a_site(contig, pos-1)

            # create the pileup processor object
            pileup_object = PileupProcessor(ref_handler, pileup_columns, contig, pos-1, rec.type, rec.alt)

            # create the image
            image_array, array_shape = pileup_object.create_image(pos-1, image_height=300, image_width=300,
                                                                  ref_band=5, alt=rec.alt, ref=rec.ref)
            # file name for the image and save the image
            file_name = contig + "_" + str(rec.pos) + "_" + str(rec.ref) + "_" + str(rec.alt) + "_" + str(rec.type)
            pileup_object.save_image_as_png(image_array, output_dir, file_name)

            # label of the image and save the image
            label = get_label(rec.type)
            smry.write(os.path.abspath(output_dir + file_name) + ".png," + str(label) + ',' + ','.join(
                map(str, array_shape)) + ',' + str(rec.genotype_class) + '\n')

            # report progress
            if (total_generated_hom_alt+total_generated_hom+total_generated_het) % 100 == 0:
                total = (total_generated_hom_alt+total_generated_hom+total_generated_het)
                sys.stderr.write(str(total) + ' variants processed in region ' + str(contig) + " " + str(site) + "\n")

    # print some stats
    sys.stderr.write('IN REGION: ' + str(contig) + ' ' + site + '\n')
    sys.stderr.write('TOTAL IN RECORDS:\n' + 'HOM\t' + 'HET\t' + 'HOM_ALT\t' + '\n')
    sys.stderr.write(str(total_hom) + '\t' + str(total_het) + '\t' + str(total_homalt) + '\n')

    sys.stderr.write('TOTAL GENERATED:\n' + 'HOM\t' + 'HET\t' + 'HOM_ALT' + '\n')
    sys.stderr.write(str(total_generated_hom) + '\t' + str(total_generated_het) + '\t'
                     + str(total_generated_hom_alt) + '\n')
    sys.stderr.write('TOTAL TIME: ' + str(time.time()-start_time) + '\n')


def get_records_given_contig(contig, vcf_file):
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
    # create the vcf handler
    vcf_handler = modules.vcf_handler.VCFFileProcessor(vcf_file)
    # generate dictionary of the region
    vcf_handler.populate_dictionary(contig, hom_filter=False)

    # get the vcf dictionary of that region
    vcf_dict = vcf_handler.get_variant_dictionary()

    # get the odds of selecting a homozygous case
    total_hom, total_het, total_homalt = vcf_handler.get_genotype_counts()
    downsample_rate = get_downsample_rate(total_hom, total_het, total_homalt)

    # keep count of how many images of each type is generated
    total_generated_hom, total_generated_het, total_generated_hom_alt = 0, 0, 0
    # generate_count = 0

    selected_records = []
    for pos in vcf_dict.keys():
        for rec in vcf_dict[pos]:
            alt = '.'
            if rec.type == 'Hom' and select_or_not(downsample_rate) is False:
                continue
            elif rec.type == 'Hom':
                rec.alt = alt

            total_generated_hom += 1 if rec.type == 'Hom' else 0
            total_generated_het += 1 if rec.type == 'Het' else 0
            total_generated_hom_alt += 1 if rec.type == 'Hom_alt' else 0
            # selected regions
            selected_records.append(rec)

    # print some stats
    sys.stderr.write('IN CHROMOSOME: ' + str(contig) + '\n')
    sys.stderr.write('TOTAL IN RECORDS:\n' + 'HOM\t' + 'HET\t' + 'HOM_ALT\t' + '\n')
    sys.stderr.write(str(total_hom) + '\t' + str(total_het) + '\t' + str(total_homalt) + '\n')

    sys.stderr.write('TOTAL PICKED:\n' + 'HOM\t' + 'HET\t' + 'HOM_ALT' + '\n')
    sys.stderr.write(str(total_generated_hom) + '\t' + str(total_generated_het) + '\t'
                     + str(total_generated_hom_alt) + '\n')

    return selected_records, total_generated_het+total_generated_hom+total_generated_hom_alt


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
    # create ref and bam files handler
    ref_handler = modules.ref_handler.RefFileProcessor(ref_file)
    bam_handler = modules.bam_handler_mpileup.BamProcessor(bam_file)

    # create a summary file
    smry = open(output_dir + "summary/" + "summary" + '_' + contig + "_" + thread_name + ".csv", 'w')
    counter = 0
    st_time = time.time()
    for rec in records:
        pos = rec.pos
        alt = '.'
        if rec.type == 'Hom':
            rec.alt = alt

        # get pileup columns from bam file
        pileup_columns = bam_handler.get_pileupcolumns_aligned_to_a_site(contig, pos-1)
        # create the pileup processor object
        pileup_object = modules.pileup_creator.PileupProcessor(ref_handler, pileup_columns, contig, pos-1,
                                                               rec.type, rec.alt)

        # create the image
        image_array, array_shape = pileup_object.create_image(pos-1, image_height=300, image_width=300, ref_band=5,
                                                              alt=rec.alt, ref=rec.ref)

        # file name for the image and save the image
        file_name = contig + "_" + str(rec.pos) + "_" + str(rec.ref) + "_" + str(rec.alt) + "_" + str(rec.type)

        pileup_object.save_image_as_png(image_array, output_dir, file_name)
        # label of the image and save the image
        label = get_label(rec.type)
        smry.write(os.path.abspath(output_dir + file_name) + ".png," + str(label) + ',' + ','.join(
            map(str, array_shape)) + ',' + str(rec.genotype_class) + '\n')
        print('Counter: ', counter, time.time()-st_time)
        counter += 1


def chromosome_level_parallelization(chr_name, bam_file, ref_file, vcf_file, output_dir, max_threads):
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
    sys.stderr.write(TextColor.BLUE + "PROCESSING VCF FILE\n" + TextColor.END)
    all_records, total_records = get_records_given_contig(chr_name, vcf_file)
    sys.stderr.write(TextColor.BLUE + "STARTING TO GENERATE IMAGES\n" + TextColor.END)

    chunk_size = int(math.ceil(total_records/max_threads)) + 2
    chunk_size = 100

    for i in range(max_threads):
        st = time.time()
        start_index = i * chunk_size
        end_index = min((i + 1) * chunk_size - 1, len(all_records)-1)
        records = all_records[start_index:end_index]
        print('Starting: ', len(records))
        args = (chr_name, bam_file, ref_file, records, output_dir, str(i))

        p = multiprocessing.Process(target=generate_pileup, args=args)
        p.start()
        p.join()
        while True:
            if len(multiprocessing.active_children()) < max_threads:
                break
        print('TIME PROFILE: ', time.time()-st)
        exit()


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

    path_to_dir_smry = path_to_dir + "summary/"
    if not os.path.exists(path_to_dir_smry):
        os.mkdir(path_to_dir_smry)

    return path_to_dir


def genome_level_parallelization(bam_file, ref_file, vcf_file, output_dir_path, max_threads):
    """
    This method calls chromosome_level_parallelization for each chromosome.
    :param bam_file: BAM file path
    :param ref_file: Reference file path
    :param vcf_file: VCF file path
    :param output_dir: Output directory
    :param max_threads: Maximum number of threads to create in chromosome level
    :return: Saves a bed file
    """
    # chr_list = ["chr1", "chr2", "chr3", "chr4", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    #             "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    program_start_time = time.time()

    chr_list = ["chr3"]

    # each chormosome in list
    for chr in chr_list:
        sys.stderr.write(TextColor.BLUE + "STARTING " + str(chr) + " PROCESSES" + "\n" + TextColor.END)
        start_time = time.time()

        # create dump directory inside output directory
        output_dir = create_output_dir_for_chromosome(output_dir_path, chr)

        # do a chromosome level parallelization
        chromosome_level_parallelization(chr, bam_file, ref_file, vcf_file, output_dir, max_threads)

        end_time = time.time()
        sys.stderr.write(TextColor.PURPLE + "FINISHED " + str(chr) + " PROCESSES" + "\n")
        sys.stderr.write(TextColor.CYAN + "TIME ELAPSED: " + str(end_time - start_time) + "\n")

    # wait for the last process to end before file processing
    while True:
        if len(multiprocessing.active_children()) == 0:
            break

    for chr in chr_list:
        # here we dumped all the bed files
        path_to_dir = output_dir_path + chr + "/" + "summary/"

        concatenated_file_name = output_dir_path + chr + ".csv"

        filemanager_object = FileManager()
        # get all csv file paths from the directory
        file_paths = filemanager_object.get_file_paths_from_directory(path_to_dir)
        # dump all csv files into one
        filemanager_object.concatenate_files(file_paths, concatenated_file_name)
        # delete all temporary files
        filemanager_object.delete_files(file_paths)
        os.rmdir(path_to_dir)

    program_end_time = time.time()
    sys.stderr.write(TextColor.RED + "PROCESSED FINISHED SUCCESSFULLY" + "\n")
    sys.stderr.write(TextColor.CYAN + "TOTAL TIME FOR GENERATING ALL RESULTS: " + str(program_end_time-program_start_time) + "\n")


if __name__ == '__main__':
    """
    Processes arguments and performs tasks to generate the pileup.
    """

    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file with alignments."
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF file containing SNPs and SVs."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/",
        help="Name of output directory"
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    parser.add_argument(
        "--parallel",
        type=bool,
        default=False,
        help="If true, will use multiple threads."
    )
    parser.add_argument(
        "--test",
        type=bool,
        default=False,
        help="If true, will only test the model."
    )
    FLAGS, not_parsed_flags = parser.parse_known_args()
    # make output directory if not already created
    FLAGS.output_dir = handle_directory(FLAGS.output_dir)

    if FLAGS.parallel is True:
        genome_level_parallelization(bam_file=FLAGS.bam,
                                     ref_file=FLAGS.ref,
                                     vcf_file=FLAGS.vcf,
                                     output_dir_path=FLAGS.output_dir,
                                     max_threads=FLAGS.max_threads)
    elif FLAGS.test is True:
        test(contig=TEST_CHR,
             site=TEST_REGION,
             bam_file=FLAGS.bam,
             ref_file=FLAGS.ref,
             vcf_file=FLAGS.vcf,
             output_dir=FLAGS.output_dir)
