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
from linescan import LineScanManager


class EHooke(object):
    """Main class of the software.
    Starts with an instance of the Parameters and Image class.
    Contains the methods needed to perform the analysis"""

    def __init__(self, cell_data=True):
        self.parameters = ParametersManager()
        self.image_manager = ImageManager()
        self.segments_manager = None
        self.cell_manager = None
        self.linescan_manager = None
        self.report_manager = None
        self.working_dir = None
        self.base_path = None
        self.fluor_path = None
        self.get_cell_images = cell_data

    def load_base_image(self, filename=None):
        """Calls the load_base_image method from the ImageManager
        Can be called without a filename or by passing one as an arg
        (filename=...)"""
        if filename is None:
            filename = tkFileDialog.askopenfilename(initialdir=self.working_dir)

        self.working_dir = "/".join(filename.split("/")[:len(filename.split("/"))-1])

        self.base_path = filename

        self.image_manager.load_base_image(filename,
                                           self.parameters.imageloaderparams)

        print "Base Image Loaded"

    def compute_mask(self):
        """Calls the compute_mask method from image_manager."""

        self.image_manager.compute_mask(self.parameters.imageloaderparams)

        if self.image_manager.fluor_image is not None:
            self.load_fluor_image(self.fluor_path)

        print "Mask Computation Finished"

    def load_fluor_image(self, filename=None):
        """Calls the load_fluor_image method from the ImageManager
        Can be called without a filename or by passing one as an arg
        (filename=...)"""
        if filename is None:
            filename = tkFileDialog.askopenfilename(initialdir=self.working_dir)

        self.fluor_path = filename

        self.image_manager.load_fluor_image(self.fluor_path,
                                            self.parameters.imageloaderparams)

        print "Fluor Image Loaded"

    def load_option_image(self, filename=None):
        """Calls the load_optional_image method from the ImageManager
        Can be called without a filename or by passing on as an arg"""

        if filename is None:
            filename = tkFileDialog.askopenfilename()

        self.image_manager.load_option_image(filename)

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
        self.cell_manager.merge_cells(label_c1, label_c2,
                                      self.parameters.cellprocessingparams,
                                      self.segments_manager,
                                      self.image_manager)
        self.cell_manager.overlay_cells(self.image_manager)

        print "Merge Finished"

    def split_cells(self, label_c1):
        """Splits a previously merged cell, requires the label of cell to be
        splitted.
        Calls the split_cells method from the cell_manager instance"""
        self.cell_manager.split_cells(int(label_c1),
                                      self.parameters.cellprocessingparams,
                                      self.segments_manager,
                                      self.image_manager)
        self.cell_manager.overlay_cells(self.image_manager)

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
        self.linescan_manager = LineScanManager()

        print "Processing Cells Finished"

    def add_line_linescan(self, point_1, point_2, point_3):
        self.linescan_manager.add_line(point_1, point_2, point_3)

    def remove_line_linescan(self, line_id):
        self.linescan_manager.remove_line(line_id)

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

    def select_from_file(self, filename=None):
        if filename is None:
            filename = tkFileDialog.askopenfilename()
        o_file = open(filename)
        data = o_file.readlines()
        o_file.close()

        data = data[0].split(";")[:len(data[0].split(";")) - 1]

        for key in data:
            self.cell_manager.cells[key].selection_state = 1

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
            filename = tkFileDialog.askdirectory(initialdir=self.working_dir)
        if label is None:
            label = self.fluor_path.split("/")
            label = label[len(label) - 1].split(".")
            if len(label) > 2:
                label = label[len(label)-3] + "." + label[len(label)-2]
            else:
                label = label[len(label)-2]

        if len(self.linescan_manager.lines.keys()) > 0:
            self.linescan_manager.measure_fluorescence(
                self.image_manager.fluor_image)

        self.report_manager = ReportManager(self.parameters)
        self.report_manager.generate_report(filename, label,
                                            self.cell_manager,
                                            self.linescan_manager,
                                            self.parameters)
        if self.get_cell_images:
            self.report_manager.get_cell_images(filename, label,
                                                self.image_manager,
                                                self.cell_manager,
                                                self.parameters)

        print "Reports Generated"
