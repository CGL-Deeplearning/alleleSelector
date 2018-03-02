import torch
import os
import numpy as np
from PIL import Image, ImageOps
from torch.utils.data import Dataset

from modules.KalleR.Bed2Image_API import Bed2ImageAPI


class DataSetLoader(Dataset):
    def __init__(self, bam_file_path, fasta_file_path, predictions, transform):
        self.bam_file_path = bam_file_path
        self.fasta_file_path = fasta_file_path
        self.prediction_records = predictions
        self.transform = transform

    def __getitem__(self, index):
        record = self.prediction_records[index]

        api_object = Bed2ImageAPI(self.bam_file_path, self.fasta_file_path)
        img, img_shape = api_object.create_image(api_object.bam_handler, api_object.fasta_handler, record)

        if self.transform is not None:
            img = self.transform(img)

        return img, list(record)

    def __len__(self):
        return len(self.prediction_records)
