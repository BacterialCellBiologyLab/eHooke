"""Module containg the classes needed for the computation of each
cell, their stats and the selection/filterion of the computed cells.
Contains a class Cell that works as the template for each cell object
and a CellManager class that controls the different steps of the cell
processing."""

from math import atan2, degrees, pi
from collections import OrderedDict
import numpy as np
import matplotlib as plt
import cellprocessing as cp
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
        self.label = cell_id
        self.merged_with = 0
        self.box = None
        self.box_margin = 5
        self.lines = []
        self.outline = []
        self.neighbours = {}
        self.color_i = -1
        self.long_axis = []
        self.short_axis = []
        self.length = 0
        self.width = 0

        self.cell_mask = None
        self.perim_mask = None
        self.sept_mask = None
        self.cyto_mask = None

        self.fluor = None
        self.local_baseline = 0
        self.image = None

        self.stats = OrderedDict([("Area", 0),
                                  ("Perimeter", 0),
                                  ("Length", 0),
                                  ("Width", 0),
                                  ("Eccentricity", 0),
                                  ("Irregularity", 0),
                                  ("Neighbours", 0),
                                  ("Baseline", 0),
                                  ("Cell Median", 0),
                                  ("Membrane Median", 0),
                                  ("Septum Median", 0),
                                  ("Cytoplasm Median", 0),
                                  ("Fluor Ratio", 0),
                                  ("Fluor Ratio 75%", 0),
                                  ("Fluor Ratio 25%", 0),
                                  ("Fluor Ratio 10%", 0)])

        self.selection_state = 1

    def clean_cell(self):
        self.label = 0
        self.merged_with = 0
        self.box = None
        self.box_margin = 5
        self.lines = []
        self.outline = []
        self.color_i = -1
        self.long_axis = []
        self.short_axis = []
        self.length = 0
        self.width = 0

        self.cell_mask = None
        self.perim_mask = None
        self.sept_mask = None
        self.cyto_mask = None

        self.fluor = None
        self.local_baseline = 0
        self.image = None

        self.stats = {"Area": None, "Perimeter": None, "Length": None,
                      "Width": None, "Eccentricity": None, "Irregularity": None,
                      "Neighbours": None}
        self.fluor_stats = {"Baseline": None, "Cell Median": None,
                            "Membrane Median": None, "Septum Median": None,
                            "Cytoplasm Median": None, "Fluor Ratio": None,
                            "Fluor Ratio 75%": None, "Fluor Ratio 25%": None,
                            "Fluor Ratio 10%": None}

        self.selection_state = 1


    def add_line(self, y, x1, x2):
        """
        Adds a line to the cell region and updates area
        """

        self.lines.append((y, x1, x2))
        self.stats["Area"] = self.stats["Area"]+x2-x1+1

    def add_frontier_point(self, x, y, neighs):
        """
        Adds an external point. Neighs is the neighbourhood labels
        """
        # check if any neighbour not in labels
        # nlabels=np.unique(neighs[neighs <> self.label])

        nlabels = []
        notzero = []
        for line in neighs:
            for p in line:
                if p != self.label and not p in nlabels:
                    nlabels.append(p)
                    if p > 0:
                        notzero.append(p)

        if nlabels != []:
            self.outline.append((x, y))

        if notzero != []:
            for l in notzero:
                if l in self.neighbours.keys():
                    count = self.neighbours[l]
                else:
                    count = 0
                self.neighbours[l] = count+1

class CellManager(object):
    """Main class of the module. Should be used to interact with the rest of
    the modules."""
    def __init__(self, params):
        self.cells = {}
        self.merged_cells = []

        spmap = plt.cm.get_cmap("hsv", params.cellprocessingparams.cell_colors)
        self.cell_colors = spmap(np.arange(
            params.cellprocessingparams.cell_colors))

        self.base_w_cells = None
        self.fluor_w_cells = None

    def cell_regions_from_labels(self, labels):
        """creates a list of N cells assuming self.labels has consecutive values from 1 to N
        create cell regions, frontiers and neighbours from labeled regions
        presumes that cell list is created and has enough elements for all
        different labels. Each cell is at index label-1
        """

        difLabels = []
        for line in labels:
            difLabels.extend(set(line))
        difLabels = sorted(set(difLabels))[1:]

        cells = {}

        for f in difLabels:
            cells[str(int(f))] = Cell(f)

        for y in range(1, len(labels[0, :])-1):
            old_label = 0
            x1 = -1
            x2 = -1

            for x in range(1, len(labels[:, 0])-1):
                l = int(labels[x, y])

                # check if line began or ended, add line
                if l != old_label:
                    if x1 > 0:
                        x2 = x-1
                        cells[str(old_label)].add_line(y, x1, x2)
                        x1 = -1
                    if l > 0:
                        x1 = x
                    old_label = l

                # check neighbours
                if l > 0:
                    square = labels[x-1:x+2, y-1:y+2]
                    cells[str(l)].add_frontier_point(x, y, square)

        for key in cells.keys():
            cells[key].stats["Perimeter"] = len(cells[key].outline)
            cells[key].stats["Neighbours"] = len(cells[key].neighbours)

        self.cells = cells

    def overlay_cells_w_base(self, base_image):
        pass

    def overlay_cells_w_fluor(self, fluor_image):
        pass

    def compute_cells(self, params, image_manager, segments_manager):
        self.cell_regions_from_labels(segments_manager.labels)

    def merge_cells(self, c1, c2, image_manager):
        pass

    def split_cells(self, c1, image_manager):
        pass

    def compute_cell_regions(self, params, image_manager):
        pass

    def compute_cell_fluor_stats(self, params, image_manager):
        pass


    def filter_cells(self):
        pass
