"""Main module of the software.
Controls the flow of the analysis and handles the interaction of the different
modules.
Contains a single class EHooke."""

import tkFileDialog
from parameters import ParametersManager
from images import ImageManager
from segments import SegmentsManager
from cells import CellManager
from reports import ReportManager

class EHooke(object):
    """Main class of the software.
    Starts with an instance of the Parameters and Image class.
    Contains the methods needed to perform the analysis"""

    def __init__(self):
        self.parameters = ParametersManager()
        self.image_manager = ImageManager()
        self.segments_manager = None
        self.cell_manager = None
        self.report_manager = None
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

        self.segments_manager = SegmentsManager()
        self.segments_manager.compute_segments(self.parameters.
                                               imageprocessingparams,
                                               self.image_manager)

        print "Segments Computation Finished"

    def compute_cells(self):
        """Creates an instance of the CellManager class and uses the
        compute_cells_method to create a list of cells based on the labels
        computed by the SegmentsManager instance."""
        self.cell_manager = CellManager(self.parameters)
        self.cell_manager.compute_cells(self.parameters.cellprocessingparams,
                                        self.image_manager,
                                        self.segments_manager)

        print "Cells Computation Finished"

    def merge_cells(self, label_c1, label_c2):
        """Merges two cells using the merge_cells method from the cell_manager
        instance and the compute_merged_cells to create a new list of cells,
        containing a cell corresponding to the merge of the previous two."""
        self.cell_manager.merge_cells(label_c1, label_c2)
        self.cell_manager.compute_merged_cells(
            self.parameters.cellprocessingparams, self.image_manager,
            self.segments_manager)

        print "Merge Finished"

    def split_cells(self, label_c1):
        """Splits a previously merged cell, requires the label of cell to be
        splitted. Calls the split_cells method from the cell_manager instance"""
        self.cell_manager.split_cells(label_c1,
                                      self.parameters.cellprocessingparams,
                                      self.image_manager,
                                      self.segments_manager)

        print "Split Finished"

    def define_as_noise(self, label_c1, noise):
        """Method used to change the state of a cell to noise or to undo it"""
        self.cell_manager.mark_cell_as_noise(label_c1, self.image_manager,
                                             noise)

    def process_cells(self):
        """Process the list of computed cells to identify the different regions
        of each cell and computes the stats related to the fluorescence"""
        self.cell_manager.process_cells(self.parameters.cellprocessingparams,
                                        self.image_manager)

        print "Processing Cells Finished"

    def select_all_cells(self):
        """Method used to mark all the cells as selected"""
        for k in self.cell_manager.cells.keys():
            if self.cell_manager.cells[k].selection_state != 0:
                self.cell_manager.cells[k].selection_state = 1

        self.cell_manager.overlay_cells(self.image_manager)

    def reject_all_cells(self):
        """Method used to mark all the cells as rejected"""
        for k in self.cell_manager.cells.keys():
            if self.cell_manager.cells[k].selection_state != 0:
                self.cell_manager.cells[k].selection_state = -1

        self.cell_manager.overlay_cells(self.image_manager)

    def filter_cells(self):
        """Filters the cell based on the filters defined in the
        params.cellprocessingparams. Calls the filter_cells method from the
        cell_manager instance"""
        self.cell_manager.filter_cells(self.parameters.cellprocessingparams,
                                       self.image_manager)

        print "Finished Filtering Cells"

    def generate_reports(self, filename=None, label=None):
        """Generates the report files by calling the generate_report method
        from Reports"""

        if filename is None:
            filename = tkFileDialog.askdirectory()
        if label is None:
            label = self.fluor_path.split("/")
            label = label[len(label)-1].split(".")[0]

        self.report_manager = ReportManager(self.parameters)
        self.report_manager.generate_report(filename, label, self.cell_manager,
                                            self.parameters)

        print "Reports Generated"
