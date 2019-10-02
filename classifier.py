import numpy as np
from skimage.util import img_as_float
from skimage.exposure import rescale_intensity
from skimage.transform import resize as skresize
from keras.models import load_model

class CellCycleClassifier(object):

    def __init__(self):

        self.model = load_model("model")

    def preprocess_image(self, image, microscope):

        if microscope == "Epifluorescence":
            max_dim = 50

        elif microscope == "SIM":
            max_dim = 100

        h, w = image.shape
        image = image[:, int(w/2):w+1]
        h, w = image.shape

        max_h, max_w = max_dim, max_dim

        lines_to_add = max_h - h
        columns_to_add = max_w - w

        if lines_to_add%2 == 0:
            new_line = np.zeros((int(lines_to_add/2), w))
            image = np.concatenate((new_line, image, new_line), axis=0)
        else:
            new_line_top = np.zeros((int(lines_to_add/2)+1, w))
            new_line_bot = np.zeros((int(lines_to_add/2), w))
            image = np.concatenate((new_line_top, image, new_line_bot), axis=0)

        if columns_to_add%2 == 0:
            columns_to_add = np.zeros((max_dim, int(columns_to_add/2)))
            image = np.concatenate((columns_to_add, image, columns_to_add), axis=1)
        else:
            columns_to_add_left = np.zeros((max_dim, int(columns_to_add/2)+1))
            columns_to_add_right = np.zeros((max_dim, int(columns_to_add/2)))
            image = np.concatenate((columns_to_add_left, image, columns_to_add_right), axis=1)

        image = img_as_float(image)
        image = image.reshape(max_dim, max_dim, 1)

        return image

    def classify_cell(self, fluor, optional, microscope):

        fluor_img = skresize(self.preprocess_image(fluor, microscope),
                             (100, 100),
                             order=0,
                             preserve_range=True,
                             anti_aliasing=False,
                             anti_aliasing_sigma=None)

        optional_img = skresize(self.preprocess_image(optional, microscope),
                                (100, 100),
                                order=0,
                                preserve_range=True,
                                anti_aliasing=False,
                                anti_aliasing_sigma=None)

        pred = self.model.predict_classes(np.concatenate((fluor_img, optional_img), axis=1).reshape(-1, 100, 200, 1))

        return pred[0] + 1

    def classify_cells(self, image_manager, cell_manager, microscope):
        fluor = image_manager.fluor_image

        if image_manager.optional_image is not None:
            optional = image_manager.optional_image
        else:
            print("No optional image provided, using dummy optional image")
            optional = np.ones(fluor.shape)

        for k in cell_manager.cells.keys():
            cell = cell_manager.cells[k]

            x0, y0, x1, y1 = cell.box

            cell_fluor = rescale_intensity(fluor[x0:x1 + 1, y0:y1 + 1] * cell.cell_mask)
            cell_optional = rescale_intensity(optional[x0:x1 + 1, y0:y1 + 1] * cell.cell_mask)

            cell_manager.cells[k].stats["Cell Cycle Phase"] = self.classify_cell(cell_fluor, cell_optional, microscope)
