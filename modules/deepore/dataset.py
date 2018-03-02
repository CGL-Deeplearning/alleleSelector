import os
from PIL import Image, ImageOps
import numpy as np
import torch
import pandas as pd
from torch.utils.data import Dataset
from torchvision import transforms, utils
from sklearn.preprocessing import MultiLabelBinarizer
from torch.autograd import Variable

class TextColor:
    """
    Defines color codes for text used to give different mode of errors.
    """
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


class PileupDataset(Dataset):
    """
    Arguments:
        A CSV file path
    """

    def __init__(self, csv_path, transform=None):
        tmp_df = pd.read_csv(csv_path, header=None)
        assert tmp_df[0].apply(lambda x: os.path.isfile(x)).all(), \
            "Some images referenced in the CSV file were not found"

        self.mlb = MultiLabelBinarizer()
        self.transform = transform

        self.X_train = tmp_df[0]
        labelLists = []
        for label in tmp_df[1]:
            labelList = [int(x) for x in str(label)]
            labelLists.append(np.array(labelList, dtype=np.long))
        self.y_train = np.array(labelLists)
        self.coverage = tmp_df[2]
        self.window = tmp_df[3]
        self.channels = tmp_df[4]
        self.type_str = tmp_df[5]

    def __getitem__(self, index):
        file = self.X_train[index]
        img = Image.open(file)

        shape = (self.coverage[index], self.window[index], self.channels[index])
        np_array_of_img = np.array(img.getdata())

        img = np.reshape(np_array_of_img, shape)
        img = np.transpose(img, (0, 1, 2))
        if self.transform is not None:
            img = self.transform(img)
        label = torch.from_numpy(self.y_train[index])
        pic_name = self.X_train[index]
        type_str = self.type_str[index] # SNP, INDEL
        return img, label, pic_name, type_str

    def __len__(self):
        return len(self.X_train.index)