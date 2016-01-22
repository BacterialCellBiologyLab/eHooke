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

    def compute_box(self, maskshape):
        """ computes the box
        """

        points = np.asarray(self.outline)  # in two columns, x, y
        bm = self.box_margin
        w, h = maskshape
        self.box = (max(min(points[:, 0])-bm, 0),
                    max(min(points[:, 1])-bm, 0),
                    min(max(points[:, 0])+bm, w-1),
                    min(max(points[:, 1])+bm, h-1))

    def axes_from_rotation(self, x0, y0, x1, y1, rotation):
        """ sets the cell axes from the box and the rotation
        """

        # midpoints
        mx = (x1+x0) / 2
        my = (y1+y0) / 2

        # assumes long is X. This duplicates rotations but simplifies
        # using different algorithms such as brightness
        self.long_axis = [[x0, my], [x1, my]]
        self.short_axis = [[mx, y0], [mx, y1]]
        self.short_axis = \
            np.asarray(np.dot(self.short_axis, rotation.T), dtype=np.int32)
        self.long_axis = \
            np.asarray(np.dot(self.long_axis, rotation.T), dtype=np.int32)

        # check if axis fall outside area due to rounding errors
        bx0, by0, bx1, by1 = self.box
        self.short_axis[0] = \
            cp.bounded_point(bx0, bx1, by0, by1, self.short_axis[0])
        self.short_axis[1] = \
            cp.bounded_point(bx0, bx1, by0, by1, self.short_axis[1])
        self.long_axis[0] = \
            cp.bounded_point(bx0, bx1, by0, by1, self.long_axis[0])
        self.long_axis[1] = \
            cp.bounded_point(bx0, bx1, by0, by1, self.long_axis[1])

        self.stats["Length"] = \
            np.linalg.norm(self.long_axis[1]-self.long_axis[0])
        self.stats["Width"] = \
            np.linalg.norm(self.short_axis[1]-self.short_axis[0])

    def compute_axes(self, rotations, maskshape):
        """ scans rotation matrices for the narrowest rectangle
        stores the result in self.long_axis and self.short_axis, each a 2,2 array
        with one point per line (coords axes in columns)

        also computes the box for masks and images
        WARNING: Rotations cannot be empty and must include a null rotation
        """

        self.compute_box(maskshape)
        points = np.asarray(self.outline)  # in two columns, x, y
        width = len(points)+1

        for rix in range(len(rotations)/2 + 1):  # no need to do more rotations, due to symmetry
            r = rotations[rix]
            nx0, ny0, nx1, ny1, nwidth = cp.bound_rectangle(np.asarray(np.dot(points, r)))

            if nwidth < width:
                width = nwidth
                x0 = nx0
                x1 = nx1
                y0 = ny0
                y1 = ny1
                angle = rix

        self.axes_from_rotation(x0, y0, x1, y1, rotations[angle])

        if self.stats["Length"] < self.stats["Width"]:
            dum = self.stats["Length"]
            self.stats["Length"] = self.stats["Width"]
            self.stats["Width"] = dum
            dum = self.short_axis
            self.short_axis = self.long_axis
            self.long_axis = dum

        self.stats["Eccentricity"] = \
            ((self.stats["Length"]-self.stats["Width"])/(self.stats["Length"]+
                                                         self.stats["Width"]))
        self.stats["Irregularity"] = \
            (len(self.outline)/(self.stats["Area"]**0.5))

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
        """creates a list of N cells assuming self.labels has consecutive
        values from 1 to N create cell regions, frontiers and neighbours from
        labeled regions presumes that cell list is created and has enough
        elements for all different labels. Each cell is at index label-1
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

    def overlay_cells_w_base(self, base_image, clip):
        """Creates an overlay of the cells over the base image.
        Besides the base image this method also requires the clipping
        coordinates for the image"""
        x0, y0, x1, y1 = clip

        base = color.rgb2gray(img_as_float(base_image[x0:x1, y0:y1]))
        base = exposure.rescale_intensity(base)

        self.base_w_cells = cp.overlay_cells(self.cells, base,
                                             self.cell_colors)

    def overlay_cells_w_fluor(self, fluor_image):
        """Creates na overlay of the cells over the fluor image)"""
        fluor = color.rgb2gray(img_as_float(fluor_image))
        fluor = exposure.rescale_intensity(fluor)
        self.fluor_w_cells = cp.overlay_cells(self.cells, fluor,
                                              self.cell_colors)

    def overlay_cells(self, image_manager):
        """Calls the methods used to create an overlay of the cells
        over the base and fluor images"""
        self.overlay_cells_w_base(image_manager.base_image, image_manager.clip)
        self.overlay_cells_w_fluor(image_manager.fluor_image)

    def compute_box_axes(self, rotations, maskshape):
        for k in self.cells.keys():
            if self.cells[k].stats["Area"] > 0:
                self.cells[k].compute_axes(rotations, maskshape)


    def compute_cells(self, params, image_manager, segments_manager):
        """Creates a cell list that is stored on self.cells as a dict, where
        each cell id is a key of the dict.
        Also creates an overlay of the cells edges over both the base and
        fluor image.
        Requires the loading of the images and the computation of the
        segments"""

        self.cell_regions_from_labels(segments_manager.labels)
        rotations = cp.rotation_matrices(params.axial_step)

        self.compute_box_axes(rotations, image_manager.mask.shape)

        #TODO: handle auto cell merging

        for k in self.cells.keys():
            cp.assign_cell_color(self.cells[k], self.cells, self.cell_colors)

        self.overlay_cells(image_manager)

    def merge_cells(self, c1, c2, image_manager):
        pass

    def split_cells(self, c1, image_manager):
        pass

    def compute_cell_regions(self, params, image_manager):
        pass

    def compute_cell_fluor_stats(self, params, image_manager):
        pass

    def process_cells(self, params, image_manager):
        """Calls the compute_cell_regions and compute_cell_fluor_stats
        methods."""
        self.compute_cell_regions(params, image_manager)
        self.compute_cell_fluor_stats(params, image_manager)

    def filter_cells(self):
        pass
