import numpy as np
from skimage.draw import line

class FluorLine(object):
    """Class used as a template for each line of the linescan.
    The class is initialized with the coordinates of two points.
    Contains a method to measure the fluorescence along that line in the
    fluor_image and computes the Fluorescence Ratio."""

    def __init__(self, point_1, point_2, point_3):
        self.line_bg_mem = line(point_1[0], point_1[1], point_2[0], point_2[1])
        self.line_cyt_sept = line(point_2[0], point_2[1], point_3[0], point_3[1])
        self.background = None
        self.membrane = None
        self.septum = None
        self.fr = None

    def measure_fluor(self, fluor_image):
        x_line_bg_mem = self.line_bg_mem[0]
        y_line_bg_mem = self.line_bg_mem[1]

        x_line_cyt_sept = self.line_cyt_sept[0]
        y_line_cyt_sept = self.line_cyt_sept[1]

        background_fluorescence = []
        for i in range(3):
            background_fluorescence.append(fluor_image[x_line_bg_mem[i], y_line_bg_mem[i]])

        membrane_fluorescence = []
        for i in range(len(x_line_bg_mem)):
            membrane_fluorescence.append(fluor_image[x_line_bg_mem[i], y_line_bg_mem[i]])

        septum_fluorescence = []
        for i in range(len(x_line_cyt_sept)):
            septum_fluorescence.append(fluor_image[x_line_cyt_sept[i], y_line_cyt_sept[i]])

        self.background = np.mean(background_fluorescence)
        self.membrane = np.max(membrane_fluorescence) - self.background
        self.septum = np.max(septum_fluorescence) - self.background
        self.fr = self.septum / self.membrane

class LineScanManager(object):
    """Class used to perform a manual linescan analysis on a selection of cells.
    Contains the methods to draw those lines, measure the fluorescence and store
    the results."""

    def __init__(self):
        self.lines = {}

    def add_line(self, point_1, point_2, point_3):
        """Creates a line object based on the coordinates of two points,
        (x0, y0) and (x1, y1)"""
        lines_id_list = [0]

        for key in self.lines:
            lines_id_list.append(int(key))

        lines_id_list = sorted(lines_id_list)
        last_id = lines_id_list[len(lines_id_list) - 1]

        self.lines[str(last_id+1)] = FluorLine(point_1, point_2, point_3)

        return last_id+1

    def remove_line(self, line_id):
        """Removes the select line from the dict"""
        new_dict = {}

        for key in self.lines.keys():
            if key != str(line_id):
                new_dict[key] = self.lines[key]

        self.lines = new_dict

    def measure_fluorescence(self, fluor_image):
        """Method used to measure the fluorescence ratios over the defined
        lines.
        Requires the fluorescene image as the first argument."""
        for key in self.lines.keys():
            self.lines[key].measure_fluor(fluor_image)

    def overlay_lines_on_image(self, fluor_img):
        color = (245, 113, 18)

        img = fluor_img

        for ln in self.lines.keys():
            for i in range(len(ln.line_bg_mem[0])):
                img[ln.line_bg_mem[0][i], ln.line_bg_mem[1][i]] = color

            for i in range(len(ln.line_cyt_sept[0])):
                img[ln.line_cyt_sept[0][i], ln.line_cyt_sept[1][i]] = color

        return fluor_img
