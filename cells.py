"""Module containg the classes needed for the computation of each
cell, their stats and the selection/filterion of the computed cells.
Contains a class Cell that works as the template for each cell object
and a CellManager class that controls the different steps of the cell
processing."""

from collections import OrderedDict
import numpy as np
import matplotlib as plt
from copy import deepcopy
from skimage.draw import line
from skimage.measure import label
from skimage.filters import threshold_isodata
from skimage.util import img_as_float, img_as_int
from skimage import morphology, color, exposure
import cellprocessing as cp


class Cell(object):
    """Template for each cell object."""

    def __init__(self, cell_id):
        self.label = cell_id
        self.merged_with = "No"
        self.merged_list = []
        self.marked_as_noise = "No"
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
        self.membsept_mask = None

        self.fluor = None
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
                                  ("Fluor Ratio 10%", 0),
                                  ("Cell Cycle Phase", 0)])

        self.selection_state = 1

    def clean_cell(self):
        """Resets the cell to an empty instance.
        Can be used to mark the cell to discard"""
        self.label = 0
        self.merged_with = "No"
        self.merged_list = []
        self.marked_as_noise = "No"
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
                                  ("Fluor Ratio 10%", 0),
                                  ("Memb+Sept Median", 0),
                                  ("Cell Cycle Phase", 0)],)

        self.selection_state = 1

    def add_line(self, y, x1, x2, pixel_size):
        """
        Adds a line to the cell region and updates area
        """

        self.lines.append((y, x1, x2))
        self.stats["Area"] = self.stats["Area"] + (x2 - x1 + 1) * float(pixel_size) * float(pixel_size)

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
                self.neighbours[l] = count + 1

    def compute_box(self, maskshape):
        """ computes the box
        """

        points = np.asarray(self.outline)  # in two columns, x, y
        bm = self.box_margin
        w, h = maskshape
        self.box = (max(min(points[:, 0]) - bm, 0),
                    max(min(points[:, 1]) - bm, 0),
                    min(max(points[:, 0]) + bm, w - 1),
                    min(max(points[:, 1]) + bm, h - 1))

    def axes_from_rotation(self, x0, y0, x1, y1, rotation, pixel_size):
        """ sets the cell axes from the box and the rotation
        """

        # midpoints
        mx = (x1 + x0) / 2
        my = (y1 + y0) / 2

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
            np.linalg.norm(self.long_axis[1] - self.long_axis[0]) * float(pixel_size)
        self.stats["Width"] = \
            np.linalg.norm(self.short_axis[1] - self.short_axis[0]) * float(pixel_size)

    def compute_axes(self, rotations, maskshape, pixel_size):
        """ scans rotation matrices for the narrowest rectangle
        stores the result in self.long_axis and self.short_axis, each a 2,2 array
        with one point per line (coords axes in columns)

        also computes the box for masks and images
        WARNING: Rotations cannot be empty and must include a null rotation
        """

        self.compute_box(maskshape)
        points = np.asarray(self.outline)  # in two columns, x, y
        width = len(points) + 1

        # no need to do more rotations, due to symmetry
        for rix in range(int(len(rotations) / 2 + 1)):
            r = rotations[rix]
            nx0, ny0, nx1, ny1, nwidth = cp.bound_rectangle(
                np.asarray(np.dot(points, r)))

            if nwidth < width:
                width = nwidth
                x0 = nx0
                x1 = nx1
                y0 = ny0
                y1 = ny1
                angle = rix

        self.axes_from_rotation(x0, y0, x1, y1, rotations[angle], pixel_size)

        if self.stats["Length"] < self.stats["Width"]:
            dum = self.stats["Length"]
            self.stats["Length"] = self.stats["Width"]
            self.stats["Width"] = dum
            dum = self.short_axis
            self.short_axis = self.long_axis
            self.long_axis = dum

        self.stats["Eccentricity"] = \
            ((1 - ((self.stats["Width"]/2.0)**2/(self.stats["Length"]/2.0)**2))**0.5)
        self.stats["Irregularity"] = \
            (len(self.outline) * float(pixel_size) / (self.stats["Area"] ** 0.5))

    def fluor_box(self, fluor):
        """ returns box of flurescence from fluor image """

        x0, y0, x1, y1 = self.box

        return fluor[x0:x1 + 1, y0:y1 + 1]

    def compute_cell_mask(self):
        x0, y0, x1, y1 = self.box
        mask = np.zeros((x1 - x0 + 1, y1 - y0 + 1))
        for lin in self.lines:
            y, st, en = lin
            mask[st - x0:en - x0 + 1, y - y0] = 1.0
        return mask

    def compute_perim_mask(self, mask, thick):
        """returns mask for perimeter
        needs cell mask
        """
        # create mask

        eroded = morphology.binary_erosion(mask, np.ones(
            (thick * 2 - 1, thick - 1))).astype(float)
        perim = mask - eroded

        return perim

    def compute_sept_mask(self, mask, thick, septum_base, algorithm):
        """ returns mask for axis.
        needs cell mask
        """

        if algorithm == "Isodata":
            return self.compute_sept_isodata(mask, thick, septum_base)

        elif algorithm == "Box":
            return self.compute_sept_box(mask, thick)

        else:
            print("Not a a valid algorithm")

    def compute_opensept_mask(self, mask, thick, septum_base, algorithm):
        """ returns mask for axis.
        needs cell mask
        """

        if algorithm == "Isodata":
            return self.compute_opensept_isodata(mask, thick, septum_base)

        elif algorithm == "Box":
            return self.compute_sept_box(mask, thick)

        else:
            print("Not a a valid algorithm")

    def compute_sept_isodata(self, mask, thick, septum_base):
        """Method used to create the cell sept_mask using the threshold_isodata
        to separate the cytoplasm from the septum"""
        cell_mask = mask
        if septum_base:
            fluor_box = 1 - self.base_box
        else:
            fluor_box = self.fluor
        perim_mask = self.compute_perim_mask(cell_mask, thick)
        inner_mask = cell_mask - perim_mask
        inner_fluor = (inner_mask > 0) * fluor_box

        threshold = threshold_isodata(inner_fluor[inner_fluor > 0])
        interest_matrix = inner_mask * (inner_fluor > threshold)

        label_matrix = label(interest_matrix, connectivity=2)
        interest_label = 0
        interest_label_sum = 0

        for l in range(np.max(label_matrix)):
            if np.sum(img_as_float(label_matrix == l + 1)) > interest_label_sum:
                interest_label = l + 1
                interest_label_sum = np.sum(
                    img_as_float(label_matrix == l + 1))

        return img_as_float(label_matrix == interest_label)

    def compute_opensept_isodata(self, mask, thick, septum_base):
        """Method used to create the cell sept_mask using the threshold_isodata
        to separate the cytoplasm from the septum"""
        cell_mask = mask
        if septum_base:
            fluor_box = 1 - self.base_box
        else:
            fluor_box = self.fluor
        perim_mask = self.compute_perim_mask(cell_mask, thick)
        inner_mask = cell_mask - perim_mask
        inner_fluor = (inner_mask > 0) * fluor_box

        threshold = threshold_isodata(inner_fluor[inner_fluor > 0])
        interest_matrix = inner_mask * (inner_fluor > threshold)

        label_matrix = label(interest_matrix, connectivity=2)
        label_sums = []

        for l in range(np.max(label_matrix)):
            label_sums.append(np.sum(img_as_float(label_matrix == l + 1)))

        print(label_sums)

        sorted_label_sums = sorted(label_sums)

        first_label = 0
        second_label = 0

        for i in range(len(label_sums)):
            if label_sums[i] == sorted_label_sums[-1]:
                first_label = i+1
                label_sums.pop(i)
                break

        for i in range(len(label_sums)):
            if label_sums[i] == sorted_label_sums[-2]:
                second_label = i+2
                label_sums.pop(i)
                break

        if second_label != 0:
            return img_as_float((label_matrix == first_label) + (label_matrix == second_label))
        else:
            return img_as_float((label_matrix == first_label))

    def compute_sept_box(self, mask, thick):
        """Method used to create a mask of the septum based on creating a box
        around the cell and then defining the septum as being the dilated short
        axis of the box."""
        x0, y0, x1, y1 = self.box
        lx0, ly0 = self.short_axis[0]
        lx1, ly1 = self.short_axis[1]
        x, y = line(lx0 - x0, ly0 - y0, lx1 - x0, ly1 - y0)

        linmask = np.zeros((x1 - x0 + 1, y1 - y0 + 1))
        linmask[x, y] = 1
        linmask = morphology.binary_dilation(
            linmask, np.ones((thick, thick))).astype(float)

        if mask is not None:
            linmask = mask * linmask

        return linmask

    def get_outline_points(self, data):
        """Method used to obtain the outline pixels of the septum"""
        outline = []
        for x in range(0, len(data)):
            for y in range(0, len(data[x])):
                if data[x, y] == 1:
                    if x == 0 and y == 0:
                        neighs_sum = data[x, y] + data[x + 1, y] + \
                            data[x + 1, y + 1] + data[x, y + 1]
                    elif x == len(data) - 1 and y == len(data[x]) - 1:
                        neighs_sum = data[x, y] + data[x, y - 1] + \
                            data[x - 1, y - 1] + data[x - 1, y]
                    elif x == 0 and y == len(data[x]) - 1:
                        neighs_sum = data[x, y] + data[x, y - 1] + \
                            data[x + 1, y - 1] + data[x + 1, y]
                    elif x == len(data) - 1 and y == 0:
                        neighs_sum = data[x, y] + data[x - 1, y] + \
                            data[x - 1, y + 1] + data[x, y + 1]
                    elif x == 0:
                        neighs_sum = data[x, y] + data[x, y - 1] + data[x, y + 1] + \
                            data[x + 1, y - 1] + \
                            data[x + 1, y] + data[x + 1, y + 1]
                    elif x == len(data) - 1:
                        neighs_sum = data[x, y] + data[x, y - 1] + data[x, y + 1] + \
                            data[x - 1, y - 1] + \
                            data[x - 1, y] + data[x - 1, y + 1]
                    elif y == 0:
                        neighs_sum = data[x, y] + data[x - 1, y] + data[x + 1, y] + \
                            data[x - 1, y + 1] + \
                            data[x, y + 1] + data[x + 1, y + 1]
                    elif y == len(data[x]) - 1:
                        neighs_sum = data[x, y] + data[x - 1, y] + data[x + 1, y] + \
                            data[x - 1, y - 1] + \
                            data[x, y - 1] + data[x + 1, y - 1]
                    else:
                        neighs_sum = data[x, y] + data[x - 1, y] + data[x + 1, y] + data[x - 1, y - 1] + data[
                            x, y - 1] + data[x + 1, y - 1] + data[x - 1, y + 1] + data[x, y + 1] + data[x + 1, y + 1]
                    if neighs_sum != 9:
                        outline.append((x, y))
        return outline

    def compute_sept_box_fix(self, outline, maskshape):
        """Method used to create a box aroung the septum, so that the short
        axis of this box can be used to choose the pixels of the membrane
        mask that need to be removed"""
        points = np.asarray(outline)  # in two columns, x, y
        bm = self.box_margin
        w, h = maskshape
        box = (max(min(points[:, 0]) - bm, 0),
               max(min(points[:, 1]) - bm, 0),
               min(max(points[:, 0]) + bm, w - 1),
               min(max(points[:, 1]) + bm, h - 1))

        return box

    def remove_sept_from_membrane(self, maskshape):
        """Method used to remove the pixels of the septum that were still in
        the membrane"""

        # get outline points of septum mask
        septum_outline = []
        septum_mask = self.sept_mask
        septum_outline = self.get_outline_points(septum_mask)

        # compute box of the septum
        septum_box = self.compute_sept_box_fix(septum_outline, maskshape)

        # compute axis of the septum
        rotations = cp.rotation_matrices(5)
        points = np.asarray(septum_outline)  # in two columns, x, y
        width = len(points) + 1

        # no need to do more rotations, due to symmetry
        for rix in range(int(len(rotations) / 2) + 1):
            r = rotations[rix]
            nx0, ny0, nx1, ny1, nwidth = cp.bound_rectangle(
                np.asarray(np.dot(points, r)))

            if nwidth < width:
                width = nwidth
                x0 = nx0
                x1 = nx1
                y0 = ny0
                y1 = ny1
                angle = rix

        rotation = rotations[angle]

        # midpoints
        mx = (x1 + x0) / 2
        my = (y1 + y0) / 2

        # assumes long is X. This duplicates rotations but simplifies
        # using different algorithms such as brightness
        long = [[x0, my], [x1, my]]
        short = [[mx, y0], [mx, y1]]
        short = np.asarray(np.dot(short, rotation.T), dtype=np.int32)
        long = np.asarray(np.dot(long, rotation.T), dtype=np.int32)

        # check if axis fall outside area due to rounding errors
        bx0, by0, bx1, by1 = septum_box
        short[0] = cp.bounded_point(bx0, bx1, by0, by1, short[0])
        short[1] = cp.bounded_point(bx0, bx1, by0, by1, short[1])
        long[0] = cp.bounded_point(bx0, bx1, by0, by1, long[0])
        long[1] = cp.bounded_point(bx0, bx1, by0, by1, long[1])

        length = np.linalg.norm(long[1] - long[0])
        width = np.linalg.norm(short[1] - short[0])

        if length < width:
            dum = length
            length = width
            width = dum
            dum = short
            short = long
            long = dum

        # expand long axis to create a linmask
        bx0, by0 = long[0]
        bx1, by1 = long[1]

        h, w = self.sept_mask.shape
        linmask = np.zeros((h, w))

        h, w = self.sept_mask.shape[0] - 2, self.sept_mask.shape[1] - 2
        bin_factor = int(width)

        if bx1 - bx0 == 0:
            x, y = line(bx0, 0, bx0, w)
            linmask[x, y] = 1
            try:
                linmask = morphology.binary_dilation(
                    linmask, np.ones((bin_factor, bin_factor))).astype(float)
            except RuntimeError:
                bin_factor = 4
                linmask = morphology.binary_dilation(
                    linmask, np.ones((bin_factor, bin_factor))).astype(float)

        else:
            m = ((by1 - by0) / (bx1 - bx0))
            b = by0 - m * bx0

            if b < 0:
                l_y0 = 0
                l_x0 = int(-b / m)

                if h * m + b > w:
                    l_y1 = w
                    l_x1 = int((w - b) / m)

                else:
                    l_x1 = h
                    l_y1 = int(h * m + b)

            elif b > w:
                l_y0 = w
                l_x0 = int((w - b) / m)

                if h * m + b < 0:
                    l_y1 = 0
                    l_x1 = int(-b / m)

                else:
                    l_x1 = h
                    l_y1 = int((h - b) / m)

            else:
                l_x0 = 0
                l_y0 = int(b)

                if m > 0:
                    if h * m + b > w:
                        l_y1 = w
                        l_x1 = int((w - b) / m)
                    else:
                        l_x1 = h
                        l_y1 = int(h * m + b)

                elif m < 0:
                    if h * m + b < 0:
                        l_y1 = 0
                        l_x1 = int(-b / m)
                    else:
                        l_x1 = h
                        l_y1 = int(h * m + b)

                else:
                    l_x1 = h
                    l_y1 = int(b)

            x, y = line(l_x0, l_y0, l_x1, l_y1)
            linmask[x, y] = 1
            try:
                linmask = morphology.binary_dilation(
                    linmask, np.ones((bin_factor, bin_factor))).astype(float)
            except RuntimeError:
                bin_factor = 4
                linmask = morphology.binary_dilation(
                    linmask, np.ones((bin_factor, bin_factor))).astype(float)
        return img_as_float(linmask)

    def recursive_compute_sept(self, cell_mask, inner_mask_thickness,
                               septum_base, algorithm):
        try:
            self.sept_mask = self.compute_sept_mask(cell_mask,
                                                    inner_mask_thickness,
                                                    septum_base,
                                                    algorithm)
        except IndexError:
            try:
                self.recursive_compute_sept(cell_mask, inner_mask_thickness - 1,
                                        septum_base,
                                        algorithm)
            except RuntimeError:
                    self.recursive_compute_sept(cell_mask, inner_mask_thickness-1, septum_base, "Box")

    def recursive_compute_opensept(self, cell_mask, inner_mask_thickness,
                                   septum_base, algorithm):
        try:
            self.sept_mask = self.compute_opensept_mask(cell_mask,
                                                        inner_mask_thickness,
                                                        septum_base,
                                                        algorithm)
        except IndexError:
            try:
                self.recursive_compute_opensept(cell_mask, inner_mask_thickness - 1,
                                                septum_base,
                                                algorithm)
            except RuntimeError:
                self.recursive_compute_opensept(cell_mask, inner_mask_thickness-1, septum_base, "Box")

    def compute_regions(self, params, image_manager):
        """Computes each different region of the cell (whole cell, membrane,
        septum, cytoplasm) and creates their respectives masks."""
        if params.look_for_septum_in_base:
            x0, y0, x1, y1 = image_manager.clip
            self.base_box = self.fluor_box(
                image_manager.base_image[x0:x1, y0:y1])
        self.fluor = self.fluor_box(image_manager.fluor_image)

        self.cell_mask = self.compute_cell_mask()

        if params.find_septum:
            self.recursive_compute_sept(self.cell_mask,
                                        params.inner_mask_thickness,
                                        params.look_for_septum_in_base,
                                        params.septum_algorithm)

            if params.septum_algorithm == "Isodata":
                self.perim_mask = self.compute_perim_mask(self.cell_mask,
                                                          params.inner_mask_thickness)
                
                self.membsept_mask = (self.perim_mask + self.sept_mask) > 0
                linmask = self.remove_sept_from_membrane(
                    image_manager.mask.shape)
                self.cyto_mask = (self.cell_mask - self.perim_mask -
                                      self.sept_mask) > 0
                if linmask is not None:
                    old_membrane = self.perim_mask
                    self.perim_mask = (old_membrane - linmask) > 0
            else:
                self.perim_mask = (self.compute_perim_mask(self.cell_mask,
                                                           params.inner_mask_thickness) -
                                   self.sept_mask) > 0
                self.membsept_mask = (self.perim_mask + self.sept_mask) > 0
                self.cyto_mask = (self.cell_mask - self.perim_mask -
                                  self.sept_mask) > 0
        elif params.find_openseptum:
            self.recursive_compute_opensept(self.cell_mask,
                                            params.inner_mask_thickness,
                                            params.look_for_septum_in_base,
                                            params.septum_algorithm)

            if params.septum_algorithm == "Isodata":
                self.perim_mask = self.compute_perim_mask(self.cell_mask,
                                                          params.inner_mask_thickness)

                self.membsept_mask = (self.perim_mask + self.sept_mask) > 0
                linmask = self.remove_sept_from_membrane(
                    image_manager.mask.shape)
                self.cyto_mask = (self.cell_mask - self.perim_mask -
                                  self.sept_mask) > 0
                if linmask is not None:
                    old_membrane = self.perim_mask
                    self.perim_mask = (old_membrane - linmask) > 0
            else:
                self.perim_mask = (self.compute_perim_mask(self.cell_mask,
                                                           params.inner_mask_thickness) -
                                   self.sept_mask) > 0
                self.membsept_mask = (self.perim_mask + self.sept_mask) > 0
                self.cyto_mask = (self.cell_mask - self.perim_mask -
                                  self.sept_mask) > 0
        else:
            self.sept_mask = None
            self.perim_mask = self.compute_perim_mask(self.cell_mask,
                                                      params.inner_mask_thickness)
            self.cyto_mask = (self.cell_mask - self.perim_mask) > 0

    def compute_fluor_baseline(self, mask, fluor, margin):
        """mask and fluor are the global images
           NOTE: mask is 0 (black) at cells and 1 (white) outside
        """

        x0, y0, x1, y1 = self.box
        wid, hei = mask.shape
        x0 = max(x0 - margin, 0)
        y0 = max(y0 - margin, 0)
        x1 = min(x1 + margin, wid - 1)
        y1 = min(y1 + margin, hei - 1)
        mask_box = mask[x0:x1, y0:y1]

        count = 0

        inverted_mask_box = 1 - mask_box

        while count < 5:
            inverted_mask_box = morphology.binary_dilation(inverted_mask_box)
            count += 1

        mask_box = 1 - inverted_mask_box

        fluor_box = fluor[x0:x1, y0:y1]
        self.stats["Baseline"] = np.median(
            mask_box[mask_box > 0] * fluor_box[mask_box > 0])

    def measure_fluor(self, fluorbox, roi, fraction=1.0):
        """returns the median and std of  fluorescence in roi
        fluorbox has the same dimensions as the roi mask
        """
        fluorbox = fluorbox
        if roi is not None:
            bright = fluorbox * roi
            bright = bright[roi > 0.5]
            # check if not enough points

            if (len(bright) * fraction) < 1.0:
                return 0.0

            if fraction < 1:
                sortvals = np.sort(bright, axis=None)[::-1]
                sortvals = sortvals[np.nonzero(sortvals)]
                sortvals = sortvals[:int(len(sortvals) * fraction)]
                return np.median(sortvals)

            else:
                return np.median(bright)
        else:
            return 0

    def compute_fluor_stats(self, params, image_manager):
        """Computes the cell stats related to the fluorescence"""
        self.compute_fluor_baseline(image_manager.mask,
                                    image_manager.original_fluor_image,
                                    params.baseline_margin)

        fluorbox = self.fluor_box(image_manager.original_fluor_image)

        self.stats["Cell Median"] = \
            self.measure_fluor(fluorbox, self.cell_mask) - \
            self.stats["Baseline"]

        self.stats["Membrane Median"] = \
            self.measure_fluor(fluorbox, self.perim_mask) - \
            self.stats["Baseline"]

        self.stats["Cytoplasm Median"] = \
            self.measure_fluor(fluorbox, self.cyto_mask) - \
            self.stats["Baseline"]

        if params.find_septum:
            self.stats["Septum Median"] = self.measure_fluor(
                fluorbox, self.sept_mask) - self.stats["Baseline"]

            self.stats["Fluor Ratio"] = (self.measure_fluor(fluorbox, self.sept_mask) - self.stats[
                                        "Baseline"]) / (self.measure_fluor(fluorbox, self.perim_mask) - self.stats["Baseline"])

            self.stats["Fluor Ratio 75%"] = (self.measure_fluor(fluorbox, self.sept_mask, 0.75) - self.stats[
                                            "Baseline"]) / (self.measure_fluor(fluorbox, self.perim_mask) - self.stats["Baseline"])

            self.stats["Fluor Ratio 25%"] = (self.measure_fluor(fluorbox, self.sept_mask, 0.25) - self.stats[
                                            "Baseline"]) / (self.measure_fluor(fluorbox, self.perim_mask) - self.stats["Baseline"])

            self.stats["Fluor Ratio 10%"] = (self.measure_fluor(fluorbox, self.sept_mask, 0.10) - self.stats[
                                            "Baseline"]) / (self.measure_fluor(fluorbox, self.perim_mask) - self.stats["Baseline"])
            self.stats["Memb+Sept Median"] = self.measure_fluor(fluorbox, self.membsept_mask) - self.stats["Baseline"]

        else:
            self.stats["Septum Median"] = 0

            self.stats["Fluor Ratio"] = 0

            self.stats["Fluor Ratio 75%"] = 0

            self.stats["Fluor Ratio 25%"] = 0

            self.stats["Fluor Ratio 10%"] = 0

            self.stats["Memb+Sept Median"] = 0

    def set_image(self, params, images, background):
        """ creates a strip with the cell in different images
            images is a list of rgb images
            background is a grayscale image to use for the masks
        """

        x0, y0, x1, y1 = self.box
        img = color.gray2rgb(
            np.zeros((x1 - x0 + 1, (len(images) + 4) * (y1 - y0 + 1))))
        bx0 = 0
        bx1 = x1 - x0 + 1
        by0 = 0
        by1 = y1 - y0 + 1

        for im in images:
            img[bx0:bx1, by0:by1] = im[x0:x1 + 1, y0:y1 + 1]
            by0 = by0 + y1 - y0 + 1
            by1 = by1 + y1 - y0 + 1

        perim = self.perim_mask
        axial = self.sept_mask
        cyto = self.cyto_mask
        img[bx0:bx1, by0:by1] = color.gray2rgb(
            background[x0:x1 + 1, y0:y1 + 1] * self.cell_mask)
        by0 = by0 + y1 - y0 + 1
        by1 = by1 + y1 - y0 + 1
        img[bx0:bx1, by0:by1] = color.gray2rgb(
            background[x0:x1 + 1, y0:y1 + 1] * perim)
        by0 = by0 + y1 - y0 + 1
        by1 = by1 + y1 - y0 + 1
        img[bx0:bx1, by0:by1] = color.gray2rgb(
            background[x0:x1 + 1, y0:y1 + 1] * cyto)
        if params.find_septum:
            by0 = by0 + y1 - y0 + 1
            by1 = by1 + y1 - y0 + 1
            img[bx0:bx1, by0:by1] = color.gray2rgb(background[x0:x1 + 1,
                                                              y0:y1 + 1] * axial)
        self.image = img_as_int(img)

    def get_cell_image(self, fluor_image):
        x0, y0, x1, y1 = self.box
        img = fluor_image[x1:x0 - 1, y1:y0 - 1] * self.cell_mask

        return img

    def recompute_outline(self, labels):
        ids = self.merged_list
        ids.append(self.label)
        new_outline = []

        for px in self.outline:
            y, x = px
            neigh_pixels = labels[y-1:y+2, x-1:x+2].flatten()

            outline_check = False

            for val in neigh_pixels:
                if val in ids:
                    pass
                else:
                    outline_check = True
            if outline_check:
                new_outline.append(px)

        self.outline = new_outline


class CellManager(object):
    """Main class of the module. Should be used to interact with the rest of
    the modules."""

    def __init__(self, params):
        self.cells = {}
        self.original_cells = {}
        self.merged_cells = []
        self.merged_labels = None

        spmap = plt.cm.get_cmap("hsv", params.cellprocessingparams.cell_colors)
        self.cell_colors = spmap(np.arange(
            params.cellprocessingparams.cell_colors))

        self.base_w_cells = None
        self.fluor_w_cells = None
        self.optional_w_cells = None

    def clean_empty_cells(self):
        """Removes empty cell objects from the cells dict"""
        newcells = {}
        for k in self.cells.keys():
            if self.cells[k].stats["Area"] > 0:
                newcells[k] = self.cells[k]

        self.cells = newcells

    def cell_regions_from_labels(self, labels, pixel_size):
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

        for y in range(1, len(labels[0, :]) - 1):
            old_label = 0
            x1 = -1
            x2 = -1

            for x in range(1, len(labels[:, 0]) - 1):
                l = int(labels[x, y])

                # check if line began or ended, add line
                if l != old_label:
                    if x1 > 0:
                        x2 = x - 1
                        cells[str(old_label)].add_line(y, x1, x2, pixel_size)
                        x1 = -1
                    if l > 0:
                        x1 = x
                    old_label = l

                # check neighbours
                if l > 0:
                    square = labels[x - 1:x + 2, y - 1:y + 2]
                    cells[str(l)].add_frontier_point(x, y, square)

        for key in cells.keys():
            cells[key].stats["Perimeter"] = len(cells[key].outline) * float(pixel_size)
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

    def overlay_cells_w_optional(self, optional_image):
        """Creates an overlay of the cells over the optional image"""
        optional = color.rgb2gray(img_as_float(optional_image))
        optional = exposure.rescale_intensity(optional)
        self.optional_w_cells = cp.overlay_cells(self.cells, optional, self.cell_colors)

    def overlay_cells(self, image_manager):
        """Calls the methods used to create an overlay of the cells
        over the base and fluor images"""
        labels = np.zeros(image_manager.fluor_image.shape)

        for k in self.cells.keys():
            c = self.cells[k]
            labels = cp.paint_cell(c, labels, c.label)

        self.merged_labels = labels
        self.overlay_cells_w_base(image_manager.base_image, image_manager.clip)
        self.overlay_cells_w_fluor(image_manager.fluor_image)

        if image_manager.optional_image is not None:
            self.overlay_cells_w_optional(image_manager.optional_image)

    def compute_box_axes(self, rotations, maskshape, pixel_size):
        for k in self.cells.keys():
            if self.cells[k].stats["Area"] > 0:
                self.cells[k].compute_axes(rotations, maskshape, pixel_size)

    def compute_cells(self, params, image_manager, segments_manager):
        """Creates a cell list that is stored on self.cells as a dict, where
        each cell id is a key of the dict.
        Also creates an overlay of the cells edges over both the base and
        fluor image.
        Requires the loading of the images and the computation of the
        segments"""

        self.cell_regions_from_labels(segments_manager.labels, params.imageloaderparams.pixel_size)
        rotations = cp.rotation_matrices(params.cellprocessingparams.axial_step)

        self.compute_box_axes(rotations, image_manager.mask.shape, params.imageloaderparams.pixel_size)

        self.original_cells = deepcopy(self.cells)

        for k in list(self.cells.keys()):
            try:
                c = self.cells[k]
                if len(c.neighbours) > 0:

                    bestneigh = max(list(iter(c.neighbours.keys())),
                                    key=(lambda key: c.neighbours[key]))
                    bestinterface = c.neighbours[bestneigh]
                    cn = self.cells[str(int(bestneigh))]

                    if cp.check_merge(c, cn, rotations, bestinterface,
                                      image_manager.mask, params.cellprocessingparams):
                        self.merge_cells(c.label, cn.label, params, segments_manager, image_manager)
            except KeyError:
                print("Cell was already merged and deleted")

        for k in self.cells.keys():
            cp.assign_cell_color(self.cells[k], self.cells,
                                 self.cell_colors, params.imageloaderparams.pixel_size)
        self.overlay_cells(image_manager)


    def merge_cells(self, label_c1, label_c2, params, segments_manager, image_manager):
        """merges two cells"""
        label_c1 = int(label_c1)
        label_c2 = int(label_c2)
        self.cells[str(label_c2)].stats["Area"] = self.cells[str(label_c2)].stats[
            "Area"] + self.cells[str(label_c1)].stats["Area"]

        self.cells[str(label_c2)].lines.extend(self.cells[str(label_c1)].lines)

        self.cells[str(label_c2)].merged_list.append(label_c1)

        self.cells[str(label_c2)].outline.extend(self.cells[str(label_c1)].outline)
        self.cells[str(label_c2)].stats["Neighbours"] = self.cells[str(label_c2)].stats["Neighbours"] + self.cells[str(label_c1)].stats["Neighbours"] - 2

        del self.cells[str(label_c1)]

        rotations = cp.rotation_matrices(params.cellprocessingparams.axial_step)
        self.cells[str(label_c2)].compute_axes(rotations, image_manager.mask.shape, params.imageloaderparams.pixel_size)

        self.cells[str(label_c2)].recompute_outline(segments_manager.labels)
        self.cells[str(label_c2)].stats["Perimeter"] = len(self.cells[str(label_c2)].outline) * float(params.imageloaderparams.pixel_size)


        if len(self.cells[str(label_c2)].merged_list) > 0:
            self.cells[str(label_c2)].merged_with = "Yes"

    def split_cells(self, label_c1, params, segments_manager, image_manager):
        """Splits a previously merged cell."""
        merged_cells = self.cells[str(label_c1)].merged_list
        merged_cells.append(label_c1)
        del self.cells[str(label_c1)]

        rotations = cp.rotation_matrices(params.cellprocessingparams.axial_step)
        for id in merged_cells:
            id = int(id)
            self.cells[str(id)] = deepcopy(self.original_cells[str(id)])
            self.cells[str(id)].compute_axes(rotations, image_manager.mask.shape, params.imageloaderparams.pixel_size)
            self.cells[str(id)].recompute_outline(segments_manager.labels)
            self.cells[str(id)].stats["Perimeter"] = len(self.cells[str(id)].outline) * float(params.imageloaderparams.pixel_size)
        if len(self.cells[str(id)].merged_list) == 0:
                self.cells[str(id)].merged_with = "No"

        for k in self.cells.keys():
            cp.assign_cell_color(self.cells[k], self.cells,
                                 self.cell_colors, params.imageloaderparams.pixel_size)

    def mark_cell_as_noise(self, label_c1, image_manager, is_noise):
        """Used to change the selection_state of a cell to 0 (noise)
        or to revert that change if the optional param "is_noise" is marked as
        false."""
        if is_noise:
            self.cells[str(label_c1)].selection_state = 0
            self.cells[str(label_c1)].marked_as_noise = "Yes"
        else:
            self.cells[str(label_c1)].selection_state = 1
            self.cells[str(label_c1)].marked_as_noise = "No"

        self.overlay_cells(image_manager)

    def process_cells(self, params, image_manager):
        """Method used to compute the individual regions of each cell and the
        computation of the stats related to the fluorescence"""
        for k in list(self.cells.keys()):
            try:
                self.cells[k].compute_regions(params, image_manager)
                self.cells[k].compute_fluor_stats(params, image_manager)
            except TypeError:
                del self.cells[k]

        fluorgray = exposure.rescale_intensity(color.rgb2gray(img_as_float(
            image_manager.fluor_image)))
        for k in self.cells.keys():
            self.cells[k].set_image(params, [self.fluor_w_cells], fluorgray)

        self.overlay_cells(image_manager)

    def filter_cells(self, params, image_manager):
        """Gets the list of filters on the parameters [("Stat", min, max)].
        Compares each cell to the filter and only select the ones that pass the filter"""
        for k in self.cells.keys():
            if self.cells[k].selection_state != 0:
                if cp.blocked_by_filter(self.cells[k], params.cell_filters):
                    self.cells[k].selection_state = -1

        self.overlay_cells(image_manager)

    def select_cells_optional(self, signal_ratio, image_manager):
        optional_image = image_manager.optional_image
        thresh = threshold_isodata(optional_image)

        for k in self.cells.keys():
            cell = self.cells[k]
            x0, y0, x1, y1 = cell.box
            cell_mask = cell.cell_mask

            optional_cell = optional_image[x0:x1+1, y0:y1+1]
            optional_signal = (optional_cell * cell_mask) > thresh

            if np.sum(optional_signal)/np.sum(cell_mask) > signal_ratio:
                self.cells[k].selection_state = 1
