import argparse
import sys
import torch
import numpy as np
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils
from torch.autograd import Variable
from pysam import VariantFile, VariantHeader, VariantRecord
from modules.deepore.inception import Inception3
from modules.deepore.dataset_prediction import PileupDataset, TextColor
from collections import defaultdict
import operator
import math


def predict(test_file, batch_size, model_path, gpu_mode):
    prediction_dict = defaultdict(list)
    transformations = transforms.Compose([transforms.ToTensor()])

    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    test_dset = PileupDataset(test_file, transformations)
    testloader = DataLoader(test_dset,
                            batch_size=batch_size,
                            shuffle=False,
                            num_workers=72,
                            pin_memory=gpu_mode # CUDA only
                            )

    sys.stderr.write(TextColor.PURPLE + 'Data loading finished\n' + TextColor.END)
    if gpu_mode is False:
        checkpoint = torch.load(model_path, map_location = 'cpu')
        state_dict = checkpoint['state_dict']
        # print('loaded state dict:', state_dict.keys())
        # print('\nIn state dict keys there is an extra word inserted by model parallel: "module.". We remove it here:')
        from collections import OrderedDict
        new_state_dict = OrderedDict()

        for k, v in state_dict.items():
            name = k[7:]  # remove `module.`
            new_state_dict[name] = v

        model = Inception3()
        model.load_state_dict(new_state_dict)
        model.cpu()
    else:
        model = torch.load(model_path)
        model = model.cuda()

    model.eval()  # Change model to 'eval' mode (BN uses moving mean/var).

    for counter, (images, image_name, records) in enumerate(testloader):
        images = Variable(images, volatile=True)

        if gpu_mode:
            images = images.cuda()

        preds = model(images).cpu()
        for i in range(0, preds.size(0)):
            rec = records[i]
            rec_id, chr_name, pos_st, pos_end, ref, alt1, alt2, rec_type = rec.rstrip().split(' ')
            probs = preds[i].data.numpy()
            prob_hom, prob_het, prob_hom_alt = probs
            prediction_dict[rec_id].append((chr_name, pos_st, pos_end, ref, alt1, alt2, rec_type, prob_hom, prob_het, prob_hom_alt))
        sys.stderr.write(TextColor.BLUE+ " BATCHES DONE: " + str(counter+1) + "/" + str(len(testloader)) + "\n" + TextColor.END)

    return prediction_dict


def get_genotype_for_multiple_allele(records):
    ref = '.'
    st_pos = 0
    end_pos = 0
    chrm = ''
    rec_alt1 = '.'
    rec_alt2 = '.'
    alt_probs = {}
    for record in records:
        chrm = record[0]
        ref = record[3]
        st_pos = record[1]
        end_pos = record[2]
        alt1 = record[4]
        alt2 = record[5]
        if alt1 != '.' and alt2 != '.':
            rec_alt1 = alt1
            rec_alt2 = alt2
            alt_probs['both'] = (record[7:])
        else:
            alt_probs[alt1] = (record[7:])
    p00 = min(alt_probs[rec_alt1][0], alt_probs[rec_alt2][0], alt_probs['both'][0])
    p01 = min(alt_probs[rec_alt1][1], alt_probs['both'][1])
    p11 = min(alt_probs[rec_alt1][2], alt_probs['both'][2])
    p02 = min(alt_probs[rec_alt2][1], alt_probs['both'][1])
    p22 = min(alt_probs[rec_alt2][2], alt_probs['both'][2])
    p12 = alt_probs['both'][2]
    prob_list = [p00, p01, p11, p02, p22, p12]
    genotype_list = ['0/0', '0/1', '1/1', '0/2', '2/2', '1/2']
    val, index = max([(v, i) for i, v in enumerate(prob_list)])
    return chrm, st_pos, end_pos, ref, [rec_alt1, rec_alt2], genotype_list[index], val


def get_genotype_for_single_allele(records):
    for record in records:
        probs = [record[7], record[8], record[9]]
        genotype_list = ['0/0', '0/1', '1/1']
        val, index = max([(v, i) for i, v in enumerate(probs)])
        return record[0], record[1], record[2], record[3], [record[4]], genotype_list[index], val


def get_vcf_header():
    header = VariantHeader()
    items = [('ID', "PASS"),
             ('Description', "All filters passed")]
    header.add_meta(key='FILTER', items=items)
    items = [('ID', "GT"),
             ('Number', 1),
             ('Type', 'String'),
             ('Description', "Genotype")]
    header.add_meta(key='FORMAT', items=items)
    items = [('ID', "chr19"),
             ('length', 198022430)]
    header.add_meta(key='contig', items=items)
    items = [('ID', "chr3"),
             ('length', 198022430)]
    header.add_meta(key='contig', items=items)
    header.add_sample('NA12878')
    return header


def get_genotype_tuple(genotype):
    split_values = genotype.split('/')
    split_values = [int(x) for x in split_values]
    return tuple(split_values)


def get_vcf_record(vcf_file, chrm, st_pos, end_pos, ref, alts, genotype, phred_qual):
    alleles = tuple([ref]) + tuple(alts)
    genotype = get_genotype_tuple(genotype)
    end_pos = int(end_pos)+1
    st_pos = int(st_pos)
    vcf_record = vcf_file.new_record(contig=chrm, start=st_pos, stop=end_pos, id='.', qual=phred_qual,
                                     filter='PASS', alleles=alleles, GT=genotype)
    return vcf_record


def produce_vcf(prediction_dict):
    header = get_vcf_header()
    vcf = VariantFile('out.vcf', 'w', header=header)
    all_calls = []
    for rec_id in sorted(prediction_dict.keys()):
        records = prediction_dict[rec_id]

        if len(records) > 1:
            chrm, st_pos, end_pos, ref, alt_field, genotype, val = get_genotype_for_multiple_allele(records)
        else:
            chrm, st_pos, end_pos, ref, alt_field, genotype, val = get_genotype_for_single_allele(records)
        if genotype == '0/0':
             continue
        phred_qual = min(60, -10 * np.log10(1 - val) if 1-val >= 0.0000000001 else 60)
        phred_qual = math.ceil(phred_qual * 100.0) / 100.0
        all_calls.append((chrm, st_pos, end_pos, ref, alt_field, genotype, phred_qual))

    all_calls.sort(key=operator.itemgetter(1))
    for record in all_calls:
        chrm, st_pos, end_pos, ref, alt_field, genotype, phred_qual = record
        # print(chrm, pos, ref, alt_field, genotype, val)
        vcf_rec = get_vcf_record(vcf, chrm, st_pos, end_pos, ref, alt_field, genotype, phred_qual)
        # print(vcf_rec)
        vcf.write(vcf_rec)


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--test_file",
        type=str,
        required=True,
        help="Testing data description csv file.."
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        required=False,
        default=100,
        help="Batch size for testing, default is 100."
    )
    parser.add_argument(
        "--model_path",
        type=str,
        default='./CNN.pkl',
        help="Saved model path."
    )
    parser.add_argument(
        "--gpu_mode",
        type=bool,
        default=False,
        help="If true then cuda is on."
    )
    FLAGS, unparsed = parser.parse_known_args()

    prediction_dict = predict(FLAGS.test_file, FLAGS.batch_size, FLAGS.model_path, FLAGS.gpu_mode)
    produce_vcf(prediction_dict)

