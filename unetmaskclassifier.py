import os
import numpy as np
from skimage.util import img_as_float
from skimage.exposure import rescale_intensity
from keras.models import load_model
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

class UnetSegmentationClassifier(object):

    def __init__(self, model_path=None):
        self.model = None

    def segment_image(self, img, model):

        max_y = img.shape[0]
        max_x = img.shape[1]
        img = img_as_float(rescale_intensity(np.array(img.reshape(max_y, max_x, 1))))

        if model == "SIM Unet NR":
            self.model = load_model("model_sim_nr_1.h5")
        elif model == "WF Unet BF":
            self.model = load_model("model_wf_bf_1.h5")
        elif model == "WF Unet NR":
            self.model = load_model("model_wf_nr_1.h5")
        else:
            print("Not a valid model")

        fullmask = np.zeros((max_y, max_x))

        # starts with bottom right corner
        tmp = np.array(img[max_y-256:, max_x-256:]).reshape([1, 256, 256, 1])
        tmp = self.model.predict(tmp)[0] > 0.5
        tmp = tmp[:, :, 1]
        fullmask[max_y-256:, max_x-256:] = tmp

        # followed by bottom crops
        for i in range(int(img.shape[1]/256.0)):
            tmp = np.array(img[max_y-256:, i*256:(i+1)*256]).reshape([1, 256, 256, 1])
            tmp = self.model.predict(tmp)[0] > 0.5
            tmp = tmp[:, :, 1]
            fullmask[max_y-256:, i*256:(i+1)*256] = tmp

        # followed by right edge crops
        for i in range(int(img.shape[0]/256.0)):
            tmp = np.array(img[i*256:(i+1)*256, max_x-256:]).reshape([1, 256, 256, 1])
            tmp = self.model.predict(tmp)[0] > 0.5
            tmp = tmp[:, :, 1]
            fullmask[i*256:(i+1)*256, max_x-256:] = tmp

        # followed by everything else
        for i in range(int(img.shape[0]/256)):
            for ii in range(int(img.shape[1]/256)):
                tmp = np.array(img[i*256:(i+1)*256, ii*256:(ii+1)*256]).reshape([1, 256, 256, 1])
                tmp = self.model.predict(tmp)[0] > 0.5
                tmp = tmp[:, :, 1]
                fullmask[i*256:(i+1)*256, ii*256:(ii+1)*256] = tmp

        return img_as_float(fullmask)

    def create_mask(self, img):
        img = img_as_float(rescale_intensity(np.array([img.reshape(256, 256, 1)])))

        mask = self.model.predict(img)[0] >= 0.5

        return mask[:, :, 1]

