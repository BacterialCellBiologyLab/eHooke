"""Module containg the classes needed for the computation of each
cell, their stats and the selection/filterion of the computed cells.
Contains a class Cell that works as the template for each cell object
A class RegionManager that is responsible for the computation of the
cell from the regions given by the segments module. A class
SelectionManager that takes care of the selection and filtering of the
cells. A CellStatsManager class responsible for the computation of the
stats of each cell; and a CellManager class that works as the interface
of the module."""

from math import atan2, degrees, pi
import numpy as np
from scipy.optimize import curve_fit
from skimage.io import imsave
from skimage.draw import line
from skimage.measure import label
from skimage.filters import threshold_isodata
from skimage.segmentation import mark_boundaries
from skimage.util import img_as_float, img_as_int
from skimage import morphology, color, exposure

class Cell(object):
    """Template for each cell object."""
    def __init__(self, cell_id):
        self.cell_id = cell_id
        self.box = None
        self.stats = {"Area": None, "Perimeter": None, "Length": None,
                      "Width": None, "Eccentricity": None, "Irregularity": None,
                      "Neighbours": None}
        self.fluor_stats = {"Baseline": None, "Cell Median": None,
                            "Membrane Median": None, "Septum Median": None,
                            "Cytoplasm Median": None, "Fluor Ratio": None,
                            "Fluor Ratio 75%": None, "Fluor Ratio 25%": None,
                            "Fluor Ratio 10%": None}


class RegionManager(object):
    """ """
    def __init__(self):
        pass

    def cells_from_labels(self):
        pass

    def merge_cells(self):
        pass

    def split_cells(self):
        pass


class SelectionManager(object):
    """ """
    def __init__(self):
        pass

    def filter_cells(self):
        pass

    def change_selection_state(self):
        pass


class CellStatsManager(object):
    """ """
    def __init__(self):
        pass


class CellManager(object):
    """ """
    def __init__(self):
        pass

    def compute_cells(self):
        pass

    def merge_cells(self):
        pass

    def split_cells(self):
        pass

    def filter_cells(self):
        pass

    def compute_cell_regions(self):
        pass

    def compute_cell_stats(self):
        pass
