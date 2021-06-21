import numpy as np
from skimage.util import img_as_float
from skimage.exposure import rescale_intensity
from keras.models import load_model
from tkinter import filedialog as fd

class UnetSegmentationClassifier(object):

    def __init__(self, model_path=None):
        self.model = None

    def segment_image(self, img, model):
        if model == "SIM Unet NR":
            self.model = load_model("model_sim_nr_1.h5")
        elif model == "WF Unet BF":
            self.model = load_model("model_wf_bf_1.h5")
        elif model == "WF Unet NR":
            self.model = load_model("model_wf_nr_1.h5")
        else:
            print("Not a valid model")

        fullmask = np.zeros(img.shape)

        print(model)

        return fullmask

    def create_mask(self, img):
        img = img_as_float(rescale_intensity(np.array([img.reshape(256, 256, 1)])))

        mask = self.model.predict(img)[0] >= 0.5

        return mask[:, :, 1]

