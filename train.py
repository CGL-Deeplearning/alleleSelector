import argparse
import os
import sys
import math
import time
import numpy as np

import torch
import torchnet.meter as meter
import torch.nn.parallel
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils
from torch.autograd import Variable

from modules.inception import Inception3
from modules.deepore.dataset import PileupDataset, TextColor
np.set_printoptions(threshold=np.nan)


def test(data_file, batch_size, gpu_mode, trained_model, num_classes):
    transformations = transforms.Compose([transforms.ToTensor()])

    validation_data = PileupDataset(data_file, transformations)
    validation_loader = DataLoader(validation_data,
                                   batch_size=batch_size,
                                   shuffle=False,
                                   num_workers=16,
                                   pin_memory=gpu_mode
                                   )
    sys.stderr.write(TextColor.PURPLE + 'Data loading finished\n' + TextColor.END)

    model = trained_model.eval()
    if gpu_mode:
        model = model.cuda()

    # Loss
    criterion = nn.CrossEntropyLoss()

    # Test the Model
    sys.stderr.write(TextColor.PURPLE + 'Test starting\n' + TextColor.END)
    total_loss = 0
    total_images = 0
    batches_done = 0
    confusion_matrix = meter.ConfusionMeter(num_classes)
    for i, (images, labels, image_name, type_class) in enumerate(validation_loader):
        if gpu_mode is True and images.size(0) % 8 != 0:
            continue

        images = Variable(images, volatile=True)
        labels = Variable(labels, volatile=True)
        if gpu_mode:
            images = images.cuda()
            labels = labels.cuda()

        # Forward + Backward + Optimize
        outputs = model(images)
        confusion_matrix.add(outputs.data.squeeze(), labels.data.type(torch.LongTensor))
        loss = criterion(outputs.contiguous().view(-1, num_classes), labels.contiguous().view(-1))
        # Loss count
        total_images += images.size(0)
        total_loss += loss.data[0]

        batches_done += 1
        sys.stderr.write(str(confusion_matrix.conf)+"\n")
        sys.stderr.write(TextColor.BLUE+'Batches done: ' + str(batches_done) + " / " + str(len(validation_loader)) +
                         "\n" + TextColor.END)

    print('Test Loss: ' + str(total_loss/total_images))
    print('Confusion Matrix: \n', confusion_matrix.conf)

    sys.stderr.write(TextColor.YELLOW+'Test Loss: ' + str(total_loss/total_images) + "\n"+TextColor.END)
    sys.stderr.write("Confusion Matrix \n: " + str(confusion_matrix.conf) + "\n" + TextColor.END)


def save_checkpoint(state, filename):
    torch.save(state, filename)


def get_base_color(base):
    if base == 'A':
        return 250.0
    if base == 'C':
        return 100.0
    if base == 'G':
        return 180.0
    if base == 'T':
        return 30.0
    if base == '*' or 'N':
        return 5.0


def get_base_by_color(color):
    if color == 250:
        return 'A'
    if color == 100:
        return 'C'
    if color == 180:
        return 'G'
    if color == 30:
        return 'T'
    if color == 5:
        return '*'
    if color == 0:
        return ' '


def get_match_by_color(color):
    if color == 0:
        return ' '
    if color <= 50: #match
        return '.'
    else:
        return 'x' #mismatch

def get_support_by_color(color):
    if color == 254:
        return '.'
    if color == 0:
        return ' '
    if color == 152:
        return 'x'


def test_image(image, img_name):
    # base_color, base_quality_color, map_quality_color, strand_color, match_color, support_color, cigar_color
    image *= 254
    # print(image.size())
    for i in range(0,image.size(1)):
        for j in range(0, image.size(2)):
            print(get_base_by_color(math.ceil(image[0][i][j])), end='')
        print()

    for i in range(0,image.size(1)):
        for j in range(0, image.size(2)):
            print(get_support_by_color(math.ceil(image[5][i][j])), end='')
        print()


def train(train_file, validation_file, batch_size, epoch_limit, file_name, gpu_mode, num_classes=3):

    transformations = transforms.Compose([transforms.ToTensor()])

    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)
    train_data_set = PileupDataset(train_file, transformations)
    train_loader = DataLoader(train_data_set,
                              batch_size=batch_size,
                              shuffle=True,
                              num_workers=16,
                              pin_memory=gpu_mode
                              )
    sys.stderr.write(TextColor.PURPLE + 'Data loading finished\n' + TextColor.END)

    # model = Inception3()
    model = Inception3()
    # Loss and Optimizer
    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001, weight_decay=0.0001)
    start_epoch = 0

    if gpu_mode:
        model = torch.nn.DataParallel(model).cuda()

    # Train the Model
    sys.stderr.write(TextColor.PURPLE + 'Training starting\n' + TextColor.END)
    for epoch in range(start_epoch, epoch_limit, 1):
        total_loss = 0
        total_images = 0
        start_time = time.time()
        batches_done = 0
        for i, (images, labels, image_name, type) in enumerate(train_loader):
            # print(image_name[0], labels[0])
            # test_image(images[0], image_name)
            # exit()

            if gpu_mode is True and images.size(0) % 8 != 0:
                continue

            images = Variable(images)
            labels = Variable(labels)
            if gpu_mode:
                images = images.cuda()
                labels = labels.cuda()

            x = images
            y = labels

            # Forward + Backward + Optimize
            optimizer.zero_grad()
            outputs = model(x)
            loss = criterion(outputs.contiguous().view(-1, num_classes), y.contiguous().view(-1))
            loss.backward()
            optimizer.step()

            # loss count
            total_images += (x.size(0))
            total_loss += loss.data[0]
            batches_done += 1

            if batches_done % 10 == 0:
                avg_loss = total_loss / total_images if total_images else 0
                print(str(epoch + 1) + "\t" + str(i + 1) + "\t" + str(avg_loss))
                sys.stderr.write(TextColor.BLUE + "EPOCH: " + str(epoch+1) + " Batches done: " + str(batches_done)
                                 + " / " + str(len(train_loader)) + "\n" + TextColor.END)
                sys.stderr.write(TextColor.YELLOW + " Loss: " + str(avg_loss) + "\n" + TextColor.END)
                sys.stderr.write(TextColor.DARKCYAN + "Time Elapsed: " + str(time.time() - start_time) +
                                 "\n" + TextColor.END)
                start_time = time.time()

        avg_loss = total_loss/total_images if total_images else 0
        sys.stderr.write(TextColor.BLUE + "EPOCH: " + str(epoch+1)
                         + " Batches done: " + str(i+1) + "/" + str(len(train_loader)) + "\n" + TextColor.END)
        sys.stderr.write(TextColor.YELLOW + " Loss: " + str(avg_loss) + "\n" + TextColor.END)
        print(str(epoch+1) + "\t" + str(i + 1) + "\t" + str(avg_loss))

        if (i+1) % 1000 == 0:
            torch.save(model, file_name + '_checkpoint_' + str(epoch+1) + '_model.pkl')
            save_checkpoint({
                'epoch': epoch + 1,
                'state_dict': model.state_dict(),
                'optimizer': optimizer.state_dict(),
            }, file_name + '_checkpoint_' + str(epoch+1) + "." + str(i+1) + "_params.pkl")
            sys.stderr.write(TextColor.RED+" MODEL SAVED \n" + TextColor.END)

        avg_loss = total_loss / total_images if total_images else 0
        sys.stderr.write(TextColor.YELLOW + 'EPOCH: ' + str(epoch+1))
        sys.stderr.write(' Loss: ' + str(avg_loss) + "\n" + TextColor.END)

        torch.save(model, file_name + '_checkpoint_' + str(epoch+1) + '_model.pkl')
        save_checkpoint({
            'epoch': epoch + 1,
            'state_dict': model.state_dict(),
            'optimizer': optimizer.state_dict(),
        }, file_name + '_checkpoint_' + str(epoch+1) + "_params.pkl")

        # After each epoch do validation
        test(validation_file, batch_size, gpu_mode, model, num_classes)

    sys.stderr.write(TextColor.PURPLE + 'Finished training\n' + TextColor.END)

    torch.save(model, file_name+'_final_model.pkl')
    save_checkpoint({
        'epoch': epoch_limit,
        'state_dict': model.state_dict(),
        'optimizer': optimizer.state_dict(),
    }, file_name + '_final_params.pkl')
    sys.stderr.write(TextColor.PURPLE + 'Model saved as:' + file_name + '_final.pkl\n' + TextColor.END)
    sys.stderr.write(TextColor.PURPLE + 'Model parameters saved as:' + file_name + '_final_params.pkl\n' + TextColor.END)


def directory_control(file_path):
    directory = os.path.dirname(file_path)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--train_file",
        type=str,
        required=True,
        help="Training data description csv file."
    )
    parser.add_argument(
        "--validation_file",
        type=str,
        required=True,
        help="Training data description csv file."
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        required=False,
        default=100,
        help="Batch size for training, default is 100."
    )
    parser.add_argument(
        "--epoch_size",
        type=int,
        required=False,
        default=10,
        help="Epoch size for training iteration."
    )
    parser.add_argument(
        "--model_out",
        type=str,
        required=False,
        default='./model',
        help="Path and file_name to save model, default is ./model"
    )
    parser.add_argument(
        "--gpu_mode",
        type=bool,
        default=False,
        help="If true then cuda is on."
    )

    FLAGS, unparsed = parser.parse_known_args()

    directory_control(FLAGS.model_out.rpartition('/')[0]+"/")
    train(FLAGS.train_file, FLAGS.validation_file, FLAGS.batch_size, FLAGS.epoch_size, FLAGS.model_out, FLAGS.gpu_mode)


