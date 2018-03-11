import argparse
import sys
import torch
import numpy as np

from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils
from torch.autograd import Variable
import torchnet.meter as meter

from modules.deepore.dataset_prediction import PileupDataset, TextColor


def predict(test_file, batch_size, model_path, gpu_mode):
    transformations = transforms.Compose([transforms.ToTensor()])

    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    test_dset = PileupDataset(test_file, transformations)
    testloader = DataLoader(test_dset,
                            batch_size=batch_size,
                            shuffle=False,
                            num_workers=16,
                            pin_memory=gpu_mode # CUDA only
                            )

    sys.stderr.write(TextColor.PURPLE + 'Data loading finished\n' + TextColor.END)

    model = torch.load(model_path)
    if gpu_mode:
        model = model.cuda()
    model.eval()  # Change model to 'eval' mode (BN uses moving mean/var).
    smry = open("predictions_" + test_file.split('/')[-1], 'w')

    for counter, (images, image_name, records) in enumerate(testloader):
        images = Variable(images, volatile=True)

        if gpu_mode:
            images = images.cuda()

        preds = model(images)
        for i in range(0, preds.size(0)):
            print(records[i], preds[i])


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

    predict(FLAGS.test_file, FLAGS.batch_size, FLAGS.model_path, FLAGS.gpu_mode)


