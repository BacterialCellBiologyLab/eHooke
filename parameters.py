"""Module containing the parameters for the image analysis.
Creates a class that encapsulates three different classes corresponding
to the parameters of each step of the analysis"""

import ConfigParser as cp
import tkFileDialog
import numpy as np


class ParametersManager(object):
    """Main class. Encapsulates the different stages parameters
    classes"""

    def __init__(self):
        self.imageloaderparams = MaskParameters()
        self.imageprocessingparams = RegionParameters()
        self.cellprocessingparams = CellParameters()

    def load_parameters(self, filename=None):
        """Loads the parameters config file"""
        if filename is None:
            filename = tkFileDialog.askopenfilename()

        parser = cp.ConfigParser()
        parser.read(filename)

        self.imageloaderparams.load_from_parser(parser, "ImageLoader")
        self.imageprocessingparams.load_from_parser(parser,
                                                    "ImageProcessing")
        self.cellprocessingparams.load_from_parser(parser,
                                                   "CellProcessing")

    def save_parameters(self, filename=None):
        """Saves parameters from a configuration file"""

        if filename is None:
            filename = tkFileDialog.asksaveasfilename()

        parser = cp.ConfigParser()

        self.imageloaderparams.save_to_parser(parser, "ImageLoader")
        self.imageprocessingparams.save_to_parser(parser, "ImageProcessing")
        self.cellprocessingparams.save_to_parser(parser, "CellProcessing")

        cfgfile = open(filename, 'w')
        parser.write(cfgfile)
        cfgfile.close()

class MaskParameters(object):
    """Class containing the parameters needed for the image loading and mask
    creation process"""

    def __init__(self):

        # phase image parameters
        self.border = 10  # phase file, including path
        self.invert_base = False
            # if true, phase will be inverted.
            # Useful when using fluorescence or light on dark background

        self.mask_algorithms = ['Local Average', 'Isodata']
        self.mask_algorithm = 'Isodata'

        # used for local average algorithm
        self.mask_blocksize = 100  # block size for moving average
        self.mask_offset = 0.02    # offset for moving average

        # used as mask creation parameters
        self.mask_fill_holes = False
            # fill holes in enclosed regions,
            # useful if cells are not uniform dark blobs
        self.mask_closing = np.ones((5, 5))
            # matrix for removing white and black spots, if empty no removal
        self.mask_dilation = 0  # mask dilation iterations

        self.auto_align = True

        self.x_align = 0
        self.y_align = 0

    def load_from_parser(self, parser, section):
        """Loads frame parameters from a ConfigParser object of the
        configuration file. The section parameters specifies the configuration
        file section"""

        self.border = parser.get(section, "border")
        self.invert_base = parser.get(section, "invert base")
        self.mask_algorithm = parser.get(section, "mask algorithm")
        self.mask_blocksize = parser.get(section, "mask blocksize")
        self.mask_offset = parser.get(section, "mask offset")
        self.mask_fill_holes = parser.get(section, "mask fill holes")
        self.mask_closing = parser.get(section, "mask closing")
        self.mask_dilation = parser.get(section, "mask dilation")
        self.auto_align = parser.get(section, "auto align")
        self.x_align = parser.get(section, "x align")
        self.y_align = parser.get(section, "y align")

    def save_to_parser(self, parser, section):
        """Saves mask parameters to a ConfigParser object of the
        configuration file. It creates the section if it does not
        exist.
        """

        if section not in parser.sections():
            parser.add_section(section)

        parser.set(section, "border", self.border)
        parser.set(section, "invert base", self.invert_base)
        parser.set(section, "mask algorithm", self.mask_algorithm)
        parser.set(section, "mask blocksize", self.mask_blocksize)
        parser.set(section, "mask offset", self.mask_offset)
        parser.set(section, "mask fill holes", self.mask_fill_holes)
        parser.set(section, "mask closing", self.mask_closing)
        parser.set(section, "mask dilation", self.mask_dilation)
        parser.set(section, "auto align", self.auto_align)
        parser.set(section, "x align", self.x_align)
        parser.set(section, "y align", self.y_align)


class RegionParameters(object):
    """Class containing the parameters for the image processing.
    Feature, labels and cell computation"""

    def __init__(self):
        # distance peak parameters
        self.peak_min_distance = 5
        self.peak_min_height = 5
        self.peak_min_distance_from_edge = 10
        self.max_peaks = 1000

        # feature labelling parameters
        self.outline_use_base_mask = False
            # assign fixed height to all in base mask

    def load_from_parser(self, parser, section):
        """Loads frame parameters from a ConfigParser object of the
        configuration file. The section parameters specifies the
        configuration file section"""

        self.peak_min_distance = parser.get(section,
                                            "peak min distance")
        self.peak_min_height = parser.get(section, "peak min height")
        self.peak_min_distance_from_edge = \
            parser.get(section, "peak min distance from edge")
        self.max_peaks = parser.get(section, "max peaks")
        self.outline_use_base_mask = parser.get(section,
                                                "outline use base mask")

    def save_to_parser(self, parser, section):
        """Saves mask parameters to a ConfigParser object of the configuration
        file. It creates the section if it does not exist."""

        if section not in parser.sections():
            parser.add_section(section)

        parser.set(section, "peak min distance", self.peak_min_distance)
        parser.set(section, "peak min height", self.peak_min_height)
        parser.set(section, "peak min distance from edge",
                   self.peak_min_distance_from_edge)
        parser.set(section, "max peaks", self.max_peaks)
        parser.set(section, "outline use base mask",
                   self.outline_use_base_mask)

class CellParameters(object):
    """Class containing the parameters needed for the process of the cells"""

    def __init__(self):
        self.axial_step = 5

        self.find_septum = True
        self.septum_algorithms = ["Box", "Isodata", "Narrowest"]
        self.septum_algorithm = "Isodata"

        # cell filtering criteria
        self.cell_filters = []

        # cell merging parameters
        self.cell_force_merge_below = 150
        self.merge_dividing_cells = False
        self.merge_length_tolerance = 1.1
        self.merge_min_interface = 15

        # cell mask for brightness
        self.inner_mask_thickness = 4

        # margin for local baseline
        self.baseline_margin = 30

        # display
        self.cell_colors = 10

    def load_from_parser(self, parser, section):
        """Loads frame parameters from a ConfigParser object of the
        configuration file. The section parameters specifies the configuration
        file section"""

        self.axial_step = parser.get(section, "axial step")
        self.find_septum = parser.get(section, "find septum")
        self.cell_filters = parser.get(section, "cell filters")
        self.cell_force_merge_below = parser.get(section,
                                                 "cell force merge below")
        self.merge_dividing_cells = parser.get(section, "merge dividing cells")
        self.merge_length_tolerance = parser.get(section,
                                                 "merge length tolerance")
        self.merge_min_interface = parser.get(section, "merge min interface")
        self.inner_mask_thickness = parser.get(section, "inner mask thickness")
        self.baseline_margin = parser.get(section, "baseline margin")
        self.cell_colors = parser.get(section, "cell colors")

    def save_to_parser(self, parser, section):
        """Saves mask parameters to a ConfigParser object of the configuration
        file. It creates the section if it does not exist."""

        if section not in parser.sections():
            parser.add_section(section)

        parser.set(section, "axial step", self.axial_step)
        parser.set(section, "find septum", self.find_septum)
        parser.set(section, "cell filters", self.cell_filters)
        parser.set(section, "cell force merge below",
                   self.cell_force_merge_below)
        parser.set(section, "merge dividing cells", self.merge_dividing_cells)
        parser.set(section, "merge length tolerance",
                   self.merge_length_tolerance)
        parser.set(section, "merge min interface", self.merge_min_interface)
        parser.set(section, "inner mask thickness", self.inner_mask_thickness)
        parser.set(section, "baseline margin", self.baseline_margin)
        parser.set(section, "cell colors", self.cell_colors)
