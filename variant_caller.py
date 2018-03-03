import argparse
import sys
import os
import time

import torch
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils
from torch.autograd import Variable
from modules.deepore.dataloader_pred import DataSetLoader
from modules.TextColor import TextColor
from image_analyzer import *

CLASS_BY_INDEX = ["HOM", "HET", "HOM_ALT"]

def get_images_for_two_alts(rec_id, record):
    chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type = record.rstrip().split('\t')[0:7]
    rec_1 = [rec_id, chr_name, pos_start, pos_end, ref, alt1, '.', rec_type]
    rec_2 = [rec_id, chr_name, pos_start, pos_end, ref, alt2, '.', rec_type]
    rec_3 = [rec_id, chr_name, pos_start, pos_end, ref, alt1, alt2,rec_type]
    return [rec_1, rec_2, rec_3]


def get_prediction_set_from_bed(candidate_bed):
    with open(candidate_bed) as bed_file:
        bed_records = bed_file.readlines()

    prediction_set = []
    rec_id = 1
    for record in bed_records:
        chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type = record.rstrip().split('\t')[0:7]
        if alt2 != '.':
            prediction_set.extend(get_images_for_two_alts(rec_id, record))
        else:
            prediction_set.append([rec_id, chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type])
        rec_id += 1

    return prediction_set


def predict(bam_file, ref_file, prediction_set, batch_size, model_path, gpu_mode):
    transformations = transforms.Compose([transforms.ToTensor()])

    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    predict_dset = DataSetLoader(bam_file, ref_file, prediction_set, transformations)
    predict_loader = DataLoader(predict_dset,
                                batch_size=batch_size,
                                shuffle=False,
                                num_workers=1,
                                pin_memory=gpu_mode)

    model = torch.load(model_path)
    if gpu_mode:
        model = model.cuda()

    for i, (images, record) in enumerate(predict_loader):
        rec_ids, chr_names, pos_starts, pos_ends, refs, alt1s, alt2s, rec_types = tuple(record)

        images = Variable(images, volatile=True)
        if gpu_mode:
            images = images.cuda()

        preds = model(images).cpu()
        preds_numpy = preds.cpu().data.topk(1)[1].numpy().ravel().tolist()

        for j in range(0, images.size()[0]):
            pred_array = preds[j].data
            print(rec_ids[j], chr_names[j], pos_starts[j], pos_ends[j], refs[j], alt1s[j], alt2s[j], rec_types[j],
                  CLASS_BY_INDEX[preds_numpy[j]], pred_array[0], pred_array[1], pred_array[2])


def call_variants(bam_file, ref_file, candidate_bed, model_path, gpu_mode, batch_size, output_dir):
    prediction_set = get_prediction_set_from_bed(candidate_bed)

    predict(bam_file, ref_file, prediction_set, batch_size, model_path, gpu_mode)


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
        "--model_path",
        type=str,
        required=True,
        help="Path to a trained model to be used for prediction"
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    parser.add_argument(
        "--gpu_mode",
        type=bool,
        default=False,
        help="If true then model will be loaded in the GPU."
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=64,
        help="Batch size for model prediction."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="variant_call_output/",
        help="Path to output directory."
    )
    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = handle_output_directory(FLAGS.output_dir)

    call_variants(FLAGS.bam, FLAGS.ref, FLAGS.candidate_bed, FLAGS.model_path, FLAGS.gpu_mode, FLAGS.batch_size,
                  FLAGS.output_dir)