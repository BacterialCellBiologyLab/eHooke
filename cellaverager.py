import math
import numpy as np

from skimage.transform import rotate, resize
from skimage.morphology import binary_erosion
from sklearn.decomposition import PCA

from cells import CellManager


class CellAverager:
    """
    Class in charge of building an average heatmap of the fluorescence of an array of cells
    """

    def __init__(self, imgman, cellman):
        self.imgman = imgman
        self.cellman = cellman

    def process(self):

        # 1. Align cells
        self.align()

        # 2. Average them
        self.average()

    def align(self):
        for key in self.cellman.cells:

            if np.count_nonzero(self.cellman.cells[key].cell_mask) == 1:
                self.cellman.cells[key].selection_state = 0
                continue

            angle = self.calculate_rotation_angle(key)
            self.cellman.cells[key].aligned_cell_mask = rotate(self.cellman.cells[key].cell_mask, angle)
            self.cellman.cells[key].aligned_fluor_mask = rotate(
                self.cellman.cells[key].fluor * self.cellman.cells[key].cell_mask, angle)

    def average(self):
        selected_cells = [self.cellman.cells[key] for key in self.cellman.cells if
                          self.cellman.cells[key].selection_state == 1]

        mean_x = int(np.median([s.aligned_fluor_mask.shape[0] for s in selected_cells]))
        mean_y = int(np.median([s.aligned_fluor_mask.shape[1] for s in selected_cells]))

        fluor_crops_array = [resize(cell.aligned_fluor_mask, (mean_x, mean_y)) for cell in selected_cells]

        model_cell = np.zeros((mean_x, mean_y))
        for cell in fluor_crops_array:
            model_cell += cell
        model_cell /= float(len(fluor_crops_array))

        self.cellman.model_cell = model_cell

    def calculate_rotation_angle(self, ids):
        binary = self.cellman.cells[ids].fluor * self.cellman.cells[ids].cell_mask
        outline = self.calculate_cell_outline(binary)
        major_axis = self.calculate_major_axis(outline)
        return self.calculate_axis_angle(major_axis)

    @staticmethod
    def calculate_cell_outline(binary):
        outline = binary * (1 - binary_erosion(binary))

        return outline

    @staticmethod
    def calculate_major_axis(outline):
        x, y = np.nonzero(outline)
        x = [[val] for val in x]
        y = [[val] for val in y]
        coords = np.concatenate((x, y), axis=1)

        pca = PCA(n_components=1)
        pca.fit(coords)

        pos_x, pos_y = pca.mean_
        eigenvector_x, eigenvector_y = pca.components_[0]
        eigenval = pca.explained_variance_[0]

        return [[pos_x - eigenvector_x * eigenval, pos_y - eigenvector_y * eigenval],
                [pos_x + eigenvector_x * eigenval, pos_y + eigenvector_y * eigenval]]

    @staticmethod
    def calculate_axis_angle(major_axis):
        # TODO refactor, atan2 should pick correct quadrant
        x0, y0 = major_axis[0]
        x1, y1 = major_axis[1]

        if x0 - x1 == 0:
            angle = 0.0

        elif y0 - y1 == 0:
            angle = 90.0

        else:
            if y1 > y0:
                if x1 > x0:
                    direction = -1
                    opposite = x1 - x0
                    adjacent = y1 - y0
                else:
                    direction = 1
                    opposite = x0 - x1
                    adjacent = y1 - y0

            elif y0 > y1:
                if x1 > x0:
                    direction = 1
                    opposite = x1 - x0
                    adjacent = y0 - y1
                else:
                    direction = -1
                    opposite = x0 - x1
                    adjacent = y0 - y1

            angle = math.degrees(math.atan(opposite / adjacent)) * direction

        if angle != 0:
            angle = 90.0 - angle
        else:
            angle = 90

        return angle
