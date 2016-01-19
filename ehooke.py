"""Main module of the software.
Controls the flow of the analysis and handles the interaction of the different
modules.
Contains a single class EHooke."""

import tkFileDialog
from parameters import ParametersManager
from images import ImageManager
from segments import SegmentsManager

class EHooke(object):
    """Main class of the software.
    Starts with an instance of the Parameters and Image class.
    Contains the methods needed to perform the analysis"""

    def __init__(self):
        self.parameters = ParametersManager()
        self.image_manager = ImageManager()
        self.segments_manager = None
        self.cell_manager = None
        self.reports_manager = None
        self.base_path = None
        self.fluor_path = None

    def load_base_image(self, filename=None):
        """Calls the load_base_image method from the ImageManager
        Can be called without a filename or by passing one as an arg
        (filename=...)"""
        if filename is None:
            filename = tkFileDialog.askopenfilename()

        self.base_path = filename

        self.image_manager.load_base_image(filename,
                                           self.parameters.imageloaderparams)

        print "Base Image Loaded"

    def compute_mask(self):
        """Calls the compute_mask method from image_manager."""

        self.image_manager.compute_mask(self.parameters.imageloaderparams)

        print "Mask Computation Finished"

    def load_fluor_image(self, filename=None):
        """Calls the load_fluor_image method from the ImageManager
        Can be called without a filename or by passing one as an arg
        (filename=...)"""
        if filename is None:
            filename = tkFileDialog.askopenfilename()

        self.fluor_path = filename

        self.image_manager.load_fluor_image(self.fluor_path,
                                            self.parameters.imageloaderparams)

        print "Fluor Image Loaded"

    def compute_segments(self):
        """Calls the compute_segments method from Segments.
        Requires the prior loading of both the phase and fluor images and
        the computation of the mask"""

        self.segments_manager = SegmentsManager(self.image_manager)
        self.segments_manager.compute_segments(self.parameters.
                                               imageprocessingparams)

        print "Segments Computation Finished"

    def compute_cells(self):
        pass

    def merge_cells(self):
        pass

    def split_cells(self):
        pass

    def filter_cells(self):
        pass

    def generate_reports(self):
        pass

    def plot_data(self):
        pass
