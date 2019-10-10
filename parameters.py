"""Module containing the parameters for the image analysis.
Creates a class that encapsulates three different classes corresponding
to the parameters of each step of the analysis"""

import configparser as cp
from tkinter import filedialog as tkFileDialog

def check_bool(param):
    if param == "True":
        return True
    elif param == "False":
        return False
    else:
        return param

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
        self.mask_blocksize = 151  # block size for moving average
        self.mask_offset = 0.02    # offset for moving average

        # used as mask creation parameters
        self.mask_fill_holes = False
        # fill holes in enclosed regions,
        # useful if cells are not uniform dark blobs
        self.mask_closing = 1
        # matrix for removing white and black spots, if empty no removal
        self.mask_dilation = 0  # mask dilation iterations

        self.auto_align = True

        self.x_align = 0
        self.y_align = 0

    def load_from_parser(self, parser, section):
        """Loads frame parameters from a ConfigParser object of the
        configuration file. The section parameters specifies the configuration
        file section"""

        self.border = int(parser.get(section, "border"))
        self.invert_base = check_bool(parser.get(section, "invert base"))
        self.mask_algorithm = str(parser.get(section, "mask algorithm"))
        self.mask_blocksize = int(parser.get(section, "mask blocksize"))
        self.mask_offset = float(parser.get(section, "mask offset"))
        self.mask_fill_holes = check_bool(parser.get(section, "mask fill holes"))
        self.mask_closing = int(float(parser.get(section, "mask closing")))
        self.mask_dilation = int(parser.get(section, "mask dilation"))
        self.auto_align = check_bool(parser.get(section, "auto align"))
        self.x_align = int(parser.get(section, "x align"))
        self.y_align = int(parser.get(section, "y align"))

    def save_to_parser(self, parser, section):
        """Saves mask parameters to a ConfigParser object of the
        configuration file. It creates the section if it does not
        exist.
        """

        if section not in parser.sections():
            parser.add_section(section)

        parser.set(section, "border", str(self.border))
        parser.set(section, "invert base", str(self.invert_base))
        parser.set(section, "mask algorithm", str(self.mask_algorithm))
        parser.set(section, "mask blocksize", str(self.mask_blocksize))
        parser.set(section, "mask offset", str(self.mask_offset))
        parser.set(section, "mask fill holes", str(self.mask_fill_holes))
        parser.set(section, "mask closing", str(self.mask_closing))
        parser.set(section, "mask dilation", str(self.mask_dilation))
        parser.set(section, "auto align", str(self.auto_align))
        parser.set(section, "x align", str(self.x_align))
        parser.set(section, "y align", str(self.y_align))


class RegionParameters(object):
    """Class containing the parameters for the image processing.
    Feature, labels and cell computation"""

    def __init__(self):
        # distance peak parameters
        self.peak_min_distance = 5
        self.peak_min_height = 5
        self.peak_min_distance_from_edge = 10
        self.max_peaks = 10000

        # feature labelling parameters
        self.outline_use_base_mask = False
        # assign fixed height to all in base mask

    def load_from_parser(self, parser, section):
        """Loads frame parameters from a ConfigParser object of the
        configuration file. The section parameters specifies the
        configuration file section"""

        self.peak_min_distance = int(parser.get(section,
                                            "peak min distance"))
        self.peak_min_height = int(parser.get(section, "peak min height"))
        self.peak_min_distance_from_edge = \
            int(parser.get(section, "peak min distance from edge"))
        self.max_peaks = int(parser.get(section, "max peaks"))
        self.outline_use_base_mask = check_bool(parser.get(section,
                                                "outline use base mask"))

    def save_to_parser(self, parser, section):
        """Saves mask parameters to a ConfigParser object of the configuration
        file. It creates the section if it does not exist."""

        if section not in parser.sections():
            parser.add_section(section)

        parser.set(section, "peak min distance", str(self.peak_min_distance))
        parser.set(section, "peak min height", str(self.peak_min_height))
        parser.set(section, "peak min distance from edge",
                   str(self.peak_min_distance_from_edge))
        parser.set(section, "max peaks", str(self.max_peaks))
        parser.set(section, "outline use base mask",
                   str(self.outline_use_base_mask))


class CellParameters(object):
    """Class containing the parameters needed for the process of the cells"""

    def __init__(self):
        self.axial_step = 5

        self.find_septum = False
        self.look_for_septum_in_base = False
        self.septum_algorithms = ["Box", "Isodata"]
        self.septum_algorithm = "Isodata"

        # microscope options for cyphid

        self.classify_cells = False
        self.microscope = "Epifluorescence"
        self.microscope_options = ["Epifluorescence", "SIM"]


        # cell filtering criteria
        self.cell_filters = []

        # cell merging parameters
        self.cell_force_merge_below = 150
        self.merge_dividing_cells = False
        self.merge_length_tolerance = 1.1
        self.merge_min_interface = 15

        # cell selection based on optional signal
        self.signal_ratio = 0.5

        # cell mask for brightness
        self.inner_mask_thickness = 4

        # margin for local baseline
        self.baseline_margin = 30

        # display
        self.cell_colors = 10

    def process_filters(self, text):
        filters = []
        if len(text.split(")")) > 1:
            tmp_filters = text.split(")")

            for i in range(len(tmp_filters)-1):
                flt = tmp_filters[i]
                tmp_filter = flt.split("(")[1]
                values = tmp_filter.split(",")
                name = str(values[0].split("'")[1])
                mini = float(values[1])
                maxi = float(values[2])

                filters.append((name, mini, maxi))

        return filters

    def load_from_parser(self, parser, section):
        """Loads frame parameters from a ConfigParser object of the
        configuration file. The section parameters specifies the configuration
        file section"""

        self.axial_step = int(parser.get(section, "axial step"))
        self.find_septum = check_bool(parser.get(section, "find septum"))
        self.classify_cells = check_bool(parser.get(section, "classify cells"))
        self.microscope = str(parser.get(section, "microscope"))
        self.look_for_septum_in_base = check_bool(parser.get(section,
                                                  "look for septum in base"))
        self.cell_filters = self.process_filters(parser.get(section, "cell filters"))
        self.cell_force_merge_below = int(parser.get(section,
                                                 "cell force merge below"))
        self.merge_dividing_cells = check_bool(parser.get(section, "merge dividing cells"))
        self.merge_length_tolerance = float(parser.get(section,
                                                 "merge length tolerance"))
        self.merge_min_interface = int(parser.get(section, "merge min interface"))
        self.inner_mask_thickness = int(parser.get(section, "inner mask thickness"))
        self.baseline_margin = int(parser.get(section, "baseline margin"))
        self.cell_colors = int(parser.get(section, "cell colors"))
        self.signal_ratio = float(parser.get(section, "signal ratio"))

    def save_to_parser(self, parser, section):
        """Saves mask parameters to a ConfigParser object of the configuration
        file. It creates the section if it does not exist."""

        if section not in parser.sections():
            parser.add_section(section)

        parser.set(section, "axial step", str(self.axial_step))
        parser.set(section, "find septum", str(self.find_septum))
        parser.set(section, "classify cells", str(self.classify_cells))
        parser.set(section, "microscope", str(self.microscope))
        parser.set(section, "look for septum in base",
                   str(self.look_for_septum_in_base))
        parser.set(section, "cell filters", str(self.cell_filters))
        parser.set(section, "cell force merge below",
                   str(self.cell_force_merge_below))
        parser.set(section, "merge dividing cells", str(self.merge_dividing_cells))
        parser.set(section, "merge length tolerance",
                   str(self.merge_length_tolerance))
        parser.set(section, "merge min interface", str(self.merge_min_interface))
        parser.set(section, "inner mask thickness", str(self.inner_mask_thickness))
        parser.set(section, "baseline margin", str(self.baseline_margin))
        parser.set(section, "cell colors", str(self.cell_colors))
        parser.set(section, "signal ratio", str(self.signal_ratio))
