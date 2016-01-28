"""Module responsible for creating the GUI and handling the ehooke module"""

import tkMessageBox
import Tkinter as tk
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import numpy as np
from ehooke import EHooke
import cellprocessing

class Interface(object):
    """Main class of the module. Used to create the GUI"""

    def __init__(self):
        self.ehooke = EHooke()
        self.default_params = self.ehooke.parameters
        self.images = {}
        self.current_image = None

        self.cid = None
        self.event_connected = False

        self.image_buttons_width = 25
        self.status_bar_width = 40
        self.status_length = 200

        self.main_window = tk.Tk()
        self.main_window.wm_title("eHooke")

        self.top_frame = tk.Frame(self.main_window, width=1200, height=10)
        self.top_frame.pack(fill="x")

        self.middle_frame = tk.Frame(self.main_window)
        self.middle_frame.pack(fill="x")

        self.parameters_panel = tk.Frame(self.middle_frame, width=600)
        self.parameters_panel.pack(side="left", fill="y")

        self.images_frame = tk.Frame(self.middle_frame)
        self.images_frame.pack(side="right", fill="y")

        self.empty_space = tk.Label(self.images_frame, text="")
        self.empty_space.pack(side="top")

        self.fig = plt.figure(figsize=(11, 8), frameon=True)
        self.canvas = FigureCanvasTkAgg(self.fig, self.middle_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side="top")

        self.ax = plt.subplot(111)
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1)
        self.ax.axis("off")
        plt.autoscale(False)

        self.canvas.show()

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.middle_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(fill="both")

        self.status = tk.StringVar()
        self.status.set("Load Base Image")
        self.status_bar = tk.Label(self.parameters_panel, textvariable=self.status, wraplength= self.status_length)
        self.status_bar.pack(side="bottom")

        self.set_imageloader()

    def remove_coord(self, x, y):
        """"Hack" to remove the mpl coordinates"""
        return ""

    def load_parameters(self):
        """Loads a .cfg with the parameters and sets them as the default
        params"""
        self.ehooke.parameters.load_parameters()
        self.default_params = self.ehooke.parameters
        self.load_default_params_imgloader()
        self.load_default_params_segments()
        self.load_default_params_cell_computation()

    def save_parameters(self):
        """Saves the current parameters in a .cfg file"""
        self.ehooke.parameters.save_parameters()

    def load_default_params_imgloader(self):
        """Loads the default params for the image loading"""
        self.border_value.set(self.default_params.imageloaderparams.border)
        self.x_align_value.set(self.default_params.imageloaderparams.x_align)
        self.y_align_value.set(self.default_params.imageloaderparams.y_align)
        self.fluor_as_base_value.set(
            self.default_params.imageloaderparams.invert_base)
        self.mask_blocksize_value.set(
            self.default_params.imageloaderparams.mask_blocksize)
        self.mask_offset_value.set(
            self.default_params.imageloaderparams.mask_offset)
        self.mask_fillholes_value.set(
            self.default_params.imageloaderparams.mask_fill_holes)
        self.mask_closing_value.set(
            self.default_params.imageloaderparams.mask_closing.shape[0])
        self.mask_dilation_value.set(
            self.default_params.imageloaderparams.mask_dilation)

    def load_default_params_segments(self):
        """Loads the default params for the segments computation"""
        self.peak_min_distance_edge_value.set(self.default_params.imageprocessingparams.peak_min_distance_from_edge)
        self.peak_min_distance_value.set(self.default_params.imageprocessingparams.peak_min_distance)
        self.peak_min_height_value.set(self.default_params.imageprocessingparams.peak_min_height)
        self.max_peaks_value.set(self.default_params.imageprocessingparams.max_peaks)
        self.use_base_mask_value.set(self.default_params.imageprocessingparams.outline_use_base_mask)

    def load_default_params_cell_computation(self):
        """Loads the default params for the cell computation"""
        self.find_septum_checkbox_value.set(self.default_params.cellprocessingparams.find_septum)
        self.axial_step_value.set(self.default_params.imageprocessingparams.axial_step)
        self.force_merge_below_value.set(self.default_params.cellprocessingparams.cell_force_merge_below)
        self.merge_dividing_value.set(self.default_params.cellprocessingparams.merge_dividing_cells)
        self.merge_length_tolerance_value.set(self.default_params.cellprocessingparams.merge_length_tolerance)
        self.merge_min_interface_value.set(self.default_params.cellprocessingparams.merge_min_interface)
        self.membrane_thickness_value.set(self.default_params.cellprocessingparams.inner_mask_thickness)

    def show_image(self, image):
        """Method use to display the selected image on the canvas"""

        if self.current_image is None:
            self.ax.cla()
            self.ax.axis("off")

        else:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()

            self.ax.cla()
            self.ax.axis("off")

            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)

        self.current_image = image

        if image == "Base":
            x1, y1, x2, y2 = self.ehooke.image_manager.clip
            self.ax.imshow(self.images[image][x1:x2, y1:y2], cmap=cm.Greys_r)
        else:
            self.ax.imshow(self.images[image], cmap=cm.Greys_r)

        self.canvas.draw()
        self.canvas.show()

    def load_base_image(self):
        """Loads the base image"""
        self.ehooke.parameters.imageloaderparams.border = \
            self.border_value.get()
        self.ehooke.parameters.imageloaderparams.invert_base = \
            self.fluor_as_base_value.get()
        self.ehooke.load_base_image()
        self.images["Base"] = self.ehooke.image_manager.base_image

        self.show_image("Base")
        self.compute_mask_button.config(state="active")
        self.base_button.config(state="active")
        self.status.set("Base Image Loaded. Proceed with mask computation")
        self.main_window.wm_title("eHooke - Base:"+str(self.ehooke.base_path))

    def compute_mask(self):
        """Computes the mask of the cell regions"""
        self.ehooke.parameters.imageloaderparams.mask_algorithm = \
            self.mask_algorithm_value.get()
        self.ehooke.parameters.imageloaderparams.mask_blocksize = \
            self.mask_blocksize_value.get()
        self.ehooke.parameters.imageloaderparams.mask_offset = \
            self.mask_offset_value.get()
        self.ehooke.parameters.imageloaderparams.mask_fill_holes = \
            self.mask_fillholes_value.get()
        tmp_closing = self.mask_closing_value.get()
        self.ehooke.parameters.imageloaderparams.mask_closing = \
            np.ones((tmp_closing, tmp_closing))
        self.ehooke.parameters.imageloaderparams.mask_dilation = \
            self.mask_dilation_value.get()
        self.ehooke.compute_mask()
        self.images["Mask"] = self.ehooke.image_manager.mask
        self.images["Base_mask"] = self.ehooke.image_manager.base_w_mask
        self.show_image("Base_mask")
        self.load_fluorescence_button.config(state="active")
        self.mask_button.config(state="active")
        self.base_with_mask_button.config(state="active")
        self.status.set("Mask computation finished. Load Fluorescence Image")

    def load_fluor(self):
        """Loads the fluor image"""
        self.ehooke.parameters.imageloaderparams.auto_align = \
            self.auto_align_value.get()
        self.ehooke.parameters.imageloaderparams.x_align = \
            self.x_align_value.get()
        self.ehooke.parameters.imageloaderparams.y_align = \
            self.y_align_value.get()
        self.ehooke.load_fluor_image()
        self.images["Fluorescence"] = \
            self.ehooke.image_manager.fluor_image
        self.images["Fluor_mask"] = self.ehooke.image_manager.fluor_w_mask
        self.show_image("Fluor_mask")
        self.next_button.config(state="active")
        self.fluor_button.config(state="active")
        self.fluor_with_mask_button.config(state="active")
        self.status.set("Fluorescence Image Loaded. Proceed to the next step")
        self.main_window.wm_title("eHooke - Base:"+str(self.ehooke.base_path)+ " - Fluorescence: "+str(self.ehooke.fluor_path))

    def set_imageloader(self):
        """Method used to change the interface to the Image Loader Step"""
        plt.clf()
        self.ax = plt.subplot(111)
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1)
        self.ax.axis("off")
        plt.autoscale(False)
        self.ax.format_coord = self.remove_coord

        self.canvas.show()

        for w in self.top_frame.winfo_children():
            w.destroy()

        for w in self.parameters_panel.winfo_children():
            w.destroy()

        for w in self.images_frame.winfo_children():
            w.destroy()

        self.empty_space = tk.Label(self.images_frame, text="")
        self.empty_space.pack(side="top")

        self.load_base_button = tk.Button(self.top_frame,
                                          text="Load Base Image",
                                          command=self.load_base_image)
        self.load_base_button.pack(side="left")

        self.compute_mask_button = tk.Button(self.top_frame,
                                             text="Compute Mask",
                                             command=self.compute_mask)
        self.compute_mask_button.pack(side="left")
        self.compute_mask_button.config(state="disabled")

        self.load_fluorescence_button = tk.Button(self.top_frame,
                                                  text="Load Fluorescence",
                                                  command=self.load_fluor)
        self.load_fluorescence_button.pack(side="left")
        self.load_fluorescence_button.config(state="disabled")

        self.load_params_button = tk.Button(self.parameters_panel,
                                            text="Load Parameters",
                                            command=self.load_parameters)
        self.load_params_button.pack(side="top", fill="x")

        self.base_parameters_label = tk.Label(self.parameters_panel,
                                               text="Load Base Parameters:")
        self.base_parameters_label.pack(side="top", fill="x")

        self.border_frame = tk.Frame(self.parameters_panel)
        self.border_frame.pack(side="top", fill="both")
        self.border_label = tk.Label(self.border_frame, text="Border: ")
        self.border_label.pack(side="left")
        self.border_value = tk.IntVar()
        self.border_entry = tk.Entry(self.border_frame,
                                     textvariable=self.border_value, width=4)
        self.border_entry.pack(side="left")
        self.border_value.set(self.ehooke.parameters.imageloaderparams.border)

        self.fluor_as_base_frame = tk.Frame(self.parameters_panel)
        self.fluor_as_base_frame.pack(side="top", fill="x")
        self.fluor_as_base_label = tk.Label(self.fluor_as_base_frame,
                                             text="Use Fluorescence Image as Base: ")
        self.fluor_as_base_label.pack(side="left")
        self.fluor_as_base_value = tk.BooleanVar()
        self.fluor_as_base_checkbox = tk.Checkbutton(self.fluor_as_base_frame,
                                                      variable=self.fluor_as_base_value, onvalue=True, offvalue=False)
        self.fluor_as_base_checkbox.pack(side="left")
        self.fluor_as_base_value.set(self.ehooke.parameters.imageloaderparams.invert_base)

        self.mask_parameters_label = tk.Label(self.parameters_panel, text="Mask Parameters:")
        self.mask_parameters_label.pack(side="top", fill="x")

        self.mask_algorithm_frame = tk.Frame(self.parameters_panel)
        self.mask_algorithm_frame.pack(side="top", fill="x")
        self.mask_algorithm_label = tk.Label(self.mask_algorithm_frame, text="Mask Algorithm: ")
        self.mask_algorithm_label.pack(side="left")
        self.mask_algorithm_value = tk.StringVar()
        self.mask_algorithm_menu = tk.OptionMenu(self.mask_algorithm_frame, self.mask_algorithm_value,
                                                 'Local Average', 'Isodata')
        self.mask_algorithm_menu.pack(side="left")
        self.mask_algorithm_value.set("Isodata")

        self.mask_blocksize_frame = tk.Frame(self.parameters_panel)
        self.mask_blocksize_frame.pack(side="top", fill="x")
        self.mask_blocksize_label = tk.Label(self.mask_blocksize_frame, text="Local Average Blocksize: ")
        self.mask_blocksize_label.pack(side="left")
        self.mask_blocksize_value = tk.IntVar()
        self.mask_blocksize_entry = tk.Entry(self.mask_blocksize_frame, textvariable=self.mask_blocksize_value, width=4)
        self.mask_blocksize_entry.pack(side="left")
        self.mask_blocksize_value.set(self.ehooke.parameters.imageloaderparams.mask_blocksize)

        self.mask_offset_frame = tk.Frame(self.parameters_panel)
        self.mask_offset_frame.pack(side="top", fill="x")
        self.mask_offset_label = tk.Label(self.mask_offset_frame, text="Local Average Offset: ")
        self.mask_offset_label.pack(side="left")
        self.mask_offset_value = tk.DoubleVar()
        self.mask_offset_entry = tk.Entry(self.mask_offset_frame, textvariable=self.mask_offset_value, width=4)
        self.mask_offset_entry.pack(side="left")
        self.mask_offset_value.set(self.ehooke.parameters.imageloaderparams.mask_offset)

        self.mask_fillholes_frame = tk.Frame(self.parameters_panel)
        self.mask_fillholes_frame.pack(side="top", fill="x")
        self.mask_fillholes_label = tk.Label(self.mask_fillholes_frame, text="Fill Holes: ")
        self.mask_fillholes_label.pack(side="left")
        self.mask_fillholes_value = tk.BooleanVar()
        self.mask_fillholes_checkbox = tk.Checkbutton(self.mask_fillholes_frame,
                                                      variable=self.mask_fillholes_value, onvalue=True, offvalue=False)
        self.mask_fillholes_checkbox.pack(side="left")
        self.mask_fillholes_value.set(self.ehooke.parameters.imageloaderparams.mask_fill_holes)

        self.mask_closing_frame = tk.Frame(self.parameters_panel)
        self.mask_closing_frame.pack(side="top", fill="x")
        self.mask_closing_label = tk.Label(self.mask_closing_frame, text="Mask Closing: ")
        self.mask_closing_label.pack(side="left")
        self.mask_closing_value = tk.DoubleVar()
        self.mask_closing_entry = tk.Entry(self.mask_closing_frame, textvariable=self.mask_closing_value, width=4)
        self.mask_closing_entry.pack(side="left")
        self.mask_closing_value.set(self.ehooke.parameters.imageloaderparams.mask_closing[0, 0])

        self.mask_dilation_frame = tk.Frame(self.parameters_panel)
        self.mask_dilation_frame.pack(side="top", fill="x")
        self.mask_dilation_label = tk.Label(self.mask_dilation_frame, text="Mask Dilation: ")
        self.mask_dilation_label.pack(side="left")
        self.mask_dilation_value = tk.IntVar()
        self.mask_dilation_entry = tk.Entry(self.mask_dilation_frame, textvariable=self.mask_dilation_value, width=4)
        self.mask_dilation_entry.pack(side="left")
        self.mask_dilation_value.set(self.ehooke.parameters.imageloaderparams.mask_dilation)

        self.load_fluor_label = tk.Label(self.parameters_panel,
                                         text="Fluorescence Image Parameters:")
        self.load_fluor_label.pack(side="top", fill="x")

        self.auto_align_frame = tk.Frame(self.parameters_panel)
        self.auto_align_frame.pack(side="top", fill="x")
        self.auto_align_label = tk.Label(self.auto_align_frame, text="Auto-align: ")
        self.auto_align_label.pack(side="left")
        self.auto_align_value = tk.BooleanVar()
        self.auto_align_checkbox = tk.Checkbutton(self.auto_align_frame,
                                                      variable=self.auto_align_value, onvalue=True, offvalue=False)
        self.auto_align_checkbox.pack(side="left")
        self.auto_align_value.set(self.ehooke.parameters.imageloaderparams.auto_align)

        self.x_align_frame = tk.Frame(self.parameters_panel)
        self.x_align_frame.pack(side="top", fill="x")
        self.x_align_label = tk.Label(self.x_align_frame, text="X align: ")
        self.x_align_label.pack(side="left")
        self.x_align_value = tk.IntVar()
        self.x_align_entry = tk.Entry(self.x_align_frame, textvariable=self.x_align_value, width=4)
        self.x_align_entry.pack(side="left")
        self.x_align_value.set(self.ehooke.parameters.imageloaderparams.x_align)

        self.y_align_frame = tk.Frame(self.parameters_panel)
        self.y_align_frame.pack(side="top", fill="x")
        self.y_align_label = tk.Label(self.y_align_frame, text="Y align: ")
        self.y_align_label.pack(side="left")
        self.y_align_value = tk.IntVar()
        self.y_align_entry = tk.Entry(self.y_align_frame, textvariable=self.y_align_value, width=4)
        self.y_align_entry.pack(side="left")
        self.y_align_value.set(self.ehooke.parameters.imageloaderparams.y_align)

        self.imgloader_params_default_button = tk.Button(self.parameters_panel, text="Default Parameters",
                                                         command=self.load_default_params_imgloader)
        self.imgloader_params_default_button.pack(side="top", fill="x")

        self.next_button = tk.Button(self.top_frame, text="Next", command=self.set_segmentscomputation)
        self.next_button.pack(side="right")
        self.next_button.config(state="disabled")

        self.base_button = tk.Button(self.images_frame, text="Base", command=lambda: self.show_image("Base"),
                                      width=self.image_buttons_width)
        self.base_button.pack(side="top", fill="x")
        self.base_button.config(state="disabled")

        self.mask_button = tk.Button(self.images_frame, text="Mask", command=lambda: self.show_image("Mask"),
                                     width=self.image_buttons_width)
        self.mask_button.pack(side="top", fill="x")
        self.mask_button.config(state="disabled")

        self.base_with_mask_button = tk.Button(self.images_frame, text="Base with mask",
                                                command=lambda: self.show_image("Base_mask"),
                                                width=self.image_buttons_width)
        self.base_with_mask_button.pack(side="top", fill="x")
        self.base_with_mask_button.config(state="disabled")

        self.fluor_button = tk.Button(self.images_frame, text="Fluorescence",
                                      command=lambda: self.show_image("Fluorescence"), width=self.image_buttons_width)
        self.fluor_button.pack(side="top", fill="x")
        self.fluor_button.config(state="disabled")

        self.fluor_with_mask_button = tk.Button(self.images_frame, text="Fluorescence with mask",
                                                command=lambda: self.show_image("Fluor_mask"),
                                                width=self.image_buttons_width)
        self.fluor_with_mask_button.pack(side="top", fill="x")
        self.fluor_with_mask_button.config(state="disabled")

        self.status = tk.StringVar()
        self.status.set("Load Base Image")
        self.status_bar = tk.Label(self.parameters_panel, textvariable=self.status, wraplength= self.status_length)
        self.status_bar.pack(side="bottom")

    def compute_features(self):
        """Calls the compute_segments method from ehooke"""
        self.ehooke.parameters.imageprocessingparams.peak_min_distance = self.peak_min_distance_value.get()
        self.ehooke.parameters.imageprocessingparams.peak_min_height = self.peak_min_height_value.get()
        self.ehooke.parameters.imageprocessingparams.peak_min_distance_from_edge = self.peak_min_distance_edge_value.get()
        self.ehooke.parameters.imageprocessingparams.max_peaks = self.max_peaks_value.get()
        self.ehooke.parameters.imageprocessingparams.outline_use_base_mask = self.use_base_mask_value.get()
        self.ehooke.compute_segments()
        self.images["Base_features"] = self.ehooke.segments_manager.base_w_features
        self.images["Fluor_features"] = self.ehooke.segments_manager.fluor_w_features
        self.show_image("Base_features")
        self.next_button.config(state="active")
        self.base_features_button.config(state="active")
        self.fluor_features_button.config(state="active")
        self.status.set("Computation of the features finished. Proceed to the next step")

    def set_segmentscomputation(self):
        """Method used to change the interface to the Segments Computation
        Step"""
        self.ax.cla()
        self.ax.axis("off")
        self.show_image("Fluor_mask")
        self.ax.format_coord = self.remove_coord
        self.canvas.show()

        for w in self.top_frame.winfo_children():
            w.destroy()

        for w in self.parameters_panel.winfo_children():
            w.destroy()

        for w in self.images_frame.winfo_children():
            w.destroy()

        self.empty_space = tk.Label(self.images_frame, text="")
        self.empty_space.pack(side="top")

        self.status = tk.StringVar()
        self.status_bar = tk.Label(self.parameters_panel,
                                   textvariable=self.status,
                                   wraplength= self.status_length)
        self.status_bar.pack(side="bottom")
        self.status.set("Waiting for features computation")

        self.compute_features_button = tk.Button(self.top_frame,
                                                 text="Compute Features",
                                                 command=self.compute_features)
        self.compute_features_button.pack(side="left")

        self.next_button = tk.Button(self.top_frame, text="Next", command=self.set_cellcomputation)
        self.next_button.pack(side="right")
        self.next_button.config(state="disabled")

        self.back_button = tk.Button(self.top_frame, text="Back", command=self.new_analysis)
        self.back_button.pack(side="right")

        self.segments_parameters_label = tk.Label(self.parameters_panel,
                                               text="Segments Computation Parameters:")
        self.segments_parameters_label.pack(side="top", fill="x")

        self.peak_min_distance_frame = tk.Frame(self.parameters_panel)
        self.peak_min_distance_frame.pack(side="top", fill="x")
        self.peak_min_distance_label = tk.Label(self.peak_min_distance_frame,
                                                text="Peak Min Distance: ")
        self.peak_min_distance_label.pack(side="left")
        self.peak_min_distance_value = tk.IntVar()
        self.peak_min_distance_entry = tk.Entry(self.peak_min_distance_frame,
                                                textvariable=self.peak_min_distance_value, width=4)
        self.peak_min_distance_entry.pack(side="left")
        self.peak_min_distance_value.set(int(self.ehooke.parameters.imageprocessingparams.peak_min_distance))

        self.peak_min_height_frame = tk.Frame(self.parameters_panel)
        self.peak_min_height_frame.pack(side="top", fill="x")
        self.peak_min_height_label = tk.Label(self.peak_min_height_frame, text="Peak Min Height: ")
        self.peak_min_height_label.pack(side="left")
        self.peak_min_height_value = tk.IntVar()
        self.peak_min_height_entry = tk.Entry(self.peak_min_height_frame,
                                              textvariable=self.peak_min_height_value, width=4)
        self.peak_min_height_entry.pack(side="left")
        self.peak_min_height_value.set(int(self.ehooke.parameters.imageprocessingparams.peak_min_height))

        self.peak_min_distance_edge_frame = tk.Frame(self.parameters_panel)
        self.peak_min_distance_edge_frame.pack(side="top", fill="x")
        self.peak_min_distance_edge_label = tk.Label(self.peak_min_distance_edge_frame, text="Peak Min Margin: ")
        self.peak_min_distance_edge_label.pack(side="left")
        self.peak_min_distance_edge_value = tk.IntVar()
        self.peak_min_distance_edge_entry = tk.Entry(self.peak_min_distance_edge_frame,
                                                     textvariable=self.peak_min_distance_edge_value, width=4)
        self.peak_min_distance_edge_entry.pack(side="left")
        self.peak_min_distance_edge_value.set(int(
            self.ehooke.parameters.imageprocessingparams.peak_min_distance_from_edge))

        self.max_peaks_frame = tk.Frame(self.parameters_panel)
        self.max_peaks_frame.pack(side="top", fill="x")
        self.max_peaks_label = tk.Label(self.max_peaks_frame, text="Max Peaks: ")
        self.max_peaks_label.pack(side="left")
        self.max_peaks_value = tk.IntVar()
        self.max_peaks_entry = tk.Entry(self.max_peaks_frame,
                                        textvariable=self.max_peaks_value, width=5)
        self.max_peaks_entry.pack(side="left")
        self.max_peaks_value.set(int(self.ehooke.parameters.imageprocessingparams.max_peaks))

        self.use_base_mask_frame = tk.Frame(self.parameters_panel)
        self.use_base_mask_frame.pack(side="top", fill="x")
        self.use_base_mask_label = tk.Label(self.use_base_mask_frame, text="Use Base Mask: ")
        self.use_base_mask_label.pack(side="left")
        self.use_base_mask_value = tk.BooleanVar()
        self.use_base_mask_checkbox = tk.Checkbutton(self.use_base_mask_frame, variable=self.use_base_mask_value,
                                                     onvalue=True, offvalue=False)
        self.use_base_mask_checkbox.pack(side="left")
        self.use_base_mask_value.set(self.ehooke.parameters.imageprocessingparams.outline_use_base_mask)

        self.features_default_button = tk.Button(self.parameters_panel, text="Default Parameters",
                                                  command=self.load_default_params_segments)
        self.features_default_button.pack(side="top", fill="x")

        self.fluor_with_mask_button = tk.Button(self.images_frame, text="Fluorescence with Mask",
                                                command=lambda: self.show_image("Fluor_mask"),
                                                width=self.image_buttons_width)
        self.fluor_with_mask_button.pack(side="top", fill="x")

        self.fluor_features_button = tk.Button(self.images_frame, text="Fluor with Features",
                                               command=lambda: self.show_image("Fluor_features"),
                                               width=self.image_buttons_width)
        self.fluor_features_button.pack(side="top", fill="x")
        self.fluor_features_button.config(state="disabled")

        self.base_features_button = tk.Button(self.images_frame, text="Base with Features",
                                               command=lambda: self.show_image("Base_features"),
                                               width=self.image_buttons_width)
        self.base_features_button.pack(side="top", fill="x")
        self.base_features_button.config(state="disabled")

    def show_cell_info_cellcomputation(self, x, y):
        """Shows the stats of each cell on the side panel"""
        label = int(self.ehooke.cell_manager.merged_labels[int(y), int(x)])

        if 0 < label:
            stats = self.ehooke.cell_manager.cells[str(label)].stats

            self.cellid_value.set(label)
            self.merged_with_value.set(self.ehooke.cell_manager.cells[str(label)].merged_with)
            self.marked_as_noise_value.set(self.ehooke.cell_manager.cells[str(label)].marked_as_noise)
            self.area_value.set(int(stats["Area"]))
            self.perimeter_value.set(int(stats["Perimeter"]))
            self.length_value.set(int(stats["Length"]))
            self.width_value.set(int(stats["Width"]))
            self.eccentricity_value.set(float(str(stats["Eccentricity"])[0:6]))
            self.irregularity_value.set(float(str(stats["Irregularity"])[0:6]))
            self.neighbours_value.set(stats["Neighbours"])

        else:
            self.cellid_value.set(0)
            self.merged_with_value.set("No")
            self.marked_as_noise_value.set("No")
            self.area_value.set(0)
            self.perimeter_value.set(0)
            self.length_value.set(0)
            self.width_value.set(0)
            self.eccentricity_value.set(0)
            self.irregularity_value.set(0)
            self.neighbours_value.set(0)

        return ""

    def compute_cells(self):
        """Method used to compute the cells"""
        self.ehooke.parameters.imageprocessingparams.axial_step = self.axial_step_value.get()
        self.ehooke.parameters.cellprocessingparams.cell_force_merge_below = self.force_merge_below_value.get()
        self.ehooke.parameters.cellprocessingparams.merge_length_tolerance = self.merge_length_tolerance_value.get()
        self.ehooke.parameters.cellprocessingparams.merge_dividing_cells = self.merge_dividing_value.get()
        self.ehooke.parameters.cellprocessingparams.merge_min_interface = self.merge_min_interface_value.get()
        self.status.set("Computing cells...")
        self.ehooke.compute_cells()
        self.ax.format_coord = self.show_cell_info_cellcomputation

        self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells

        self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

        self.show_image("Fluor_cells_outlined")

        self.next_button.config(state="active")
        self.force_merge_button.config(state="active")
        self.split_cell_button.config(state="active")
        self.declare_as_noise_button.config(state="active")
        self.undo_as_noise_button.config(state="active")
        self.base_w_cells_button.config(state="active")
        self.fluor_cells_out_button.config(state="active")
        self.status.set("Cell Computation Finished. Proceed to the next step")

    def merge_on_press(self, event):

        if event.button == 3:
            label = int(self.ehooke.cell_manager.merged_labels[int(event.ydata), int(event.xdata)])

            if label > 0:

                if len(self.merge_list) <1:
                    self.merge_list.append(label)
                    self.status.set("Waiting for second cell")

                elif len(self.merge_list) == 1:
                    if label != self.merge_list[0]:
                        self.status.set("Merging cells...")

                        self.ehooke.merge_cells(self.merge_list[0], label)
                        self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells
                        self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

                        self.show_image(self.current_image)

                        self.canvas.mpl_disconnect(self.cid)
                        self.event_connected = False
                        self.merge_list = []

                        self.status.set("Merge Finished")
                    else:
                        self.canvas.mpl_disconnect(self.cid)
                        self.event_connected = False
                        self.status.set("Same Cell. Repeat Merge.")

            else:
                self.canvas.mpl_disconnect(self.cid)
                self.event_connected = False
                self.status.set("Not a Cell. Repeat Merge")

    def force_merge(self):
        """Method used to force the merge of two cells"""
        if self.event_connected:
            self.canvas.mpl_disconnect(self.cid)
        self.status.set("Right-click on two cells")
        self.merge_list = []
        self.cid = self.canvas.mpl_connect('button_release_event', self.merge_on_press)
        self.event_connected = True

    def splitting_on_press(self, event):
        if event.button == 3:
            self.status.set("Splitting Cells")
            label = int(self.ehooke.cell_manager.merged_labels[int(event.ydata), int(event.xdata)])

            if label > 0:

                for pair in self.ehooke.cell_manager.merged_cells:
                    if label == pair[0] or label == pair[1]:

                        self.ehooke.split_cells(label)
                        self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells
                        self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

                        self.show_image(self.current_image)

                    else:
                        pass

            self.canvas.mpl_disconnect(self.cid)
            self.event_connected = False
            self.status.set("Splitting Finished")

    def split_cell(self):
        """Method used to split a previously merged cell"""
        if self.event_connected:
            self.canvas.mpl_disconnect(self.cid)
        self.status.set("Right-click on cell")
        self.cid = self.canvas.mpl_connect('button_release_event', self.splitting_on_press)
        self.event_connected = True

    def noise_on_press(self, event):
        if event.button == 3:
            self.status.set("Removing Cell")
            label = int(self.ehooke.cell_manager.merged_labels[int(event.ydata), int(event.xdata)])

            if label > 0:
                if self.ehooke.cell_manager.cells[str(label)].selection_state != 0:
                    self.ehooke.define_as_noise(label, True)
                    self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells
                    self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

                    self.show_image(self.current_image)

            self.canvas.mpl_disconnect(self.cid)
            self.event_connected = False
            self.status.set("May Proceed to Cell Processing")

    def declare_as_noise(self):
        """Method used to define a cell object as noise"""
        if self.event_connected:
            self.canvas.mpl_disconnect(self.cid)
        self.status.set("Righ-click noise object")
        self.cid = self.canvas.mpl_connect('button_release_event', self.noise_on_press)
        self.event_connected = True

    def undo_noise_on_press(self, event):
        if event.button == 3:
            self.status.set("Adding Cell")
            label = int(self.ehooke.cell_manager.merged_labels[int(event.ydata), int(event.xdata)])

            if label > 0:
                if self.ehooke.cell_manager.cells[str(label)].selection_state == 0:
                    self.ehooke.define_as_noise(label, False)
                    self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells
                    self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

                    self.show_image(self.current_image)

            self.canvas.mpl_disconnect(self.cid)
            self.event_connected = False
            self.status.set("May Proceed to Cell Processing")

    def undo_as_noise(self):
        """Method used to define a a cell from an object that was previously
        defined as noise"""
        if self.event_connected:
            self.canvas.mpl_disconnect(self.cid)
        self.status.set("Righ-click cell object")
        self.cid = self.canvas.mpl_connect('button_release_event', self.undo_noise_on_press)
        self.event_connected = True

    def set_cellcomputation(self):
        """Method used to change the interface to the Cell Computation
        Step"""
        self.ax.cla()
        self.ax.axis("off")
        self.show_image("Base_features")
        self.ax.format_coord = self.remove_coord
        self.canvas.show()

        for w in self.top_frame.winfo_children():
            w.destroy()

        for w in self.parameters_panel.winfo_children():
            w.destroy()

        for w in self.images_frame.winfo_children():
            w.destroy()

        self.empty_space = tk.Label(self.images_frame, text="")
        self.empty_space.pack(side="top")

        self.status = tk.StringVar()
        self.status_bar = tk.Label(self.parameters_panel, textvariable=self.status, wraplength= self.status_length)
        self.status_bar.pack(side="bottom")
        self.status.set("Waiting for cell computation")

        self.compute_cells_button = tk.Button(self.top_frame, text="Compute Cells", command=self.compute_cells)
        self.compute_cells_button.pack(side="left")

        self.next_button = tk.Button(self.top_frame, text="Next", command=self.set_cellprocessing)
        self.next_button.pack(side="right")
        self.next_button.config(state="disabled")
        self.back_button = tk.Button(self.top_frame, text="Back", command=self.set_segmentscomputation)
        self.back_button.pack(side="right")

        self.cellcomputation_parameters_label = tk.Label(self.parameters_panel,
                                               text="Cell Computation Parameters:")
        self.cellcomputation_parameters_label.pack(side="top", fill="x")

        self.axial_step_frame = tk.Frame(self.parameters_panel)
        self.axial_step_frame.pack(side="top", fill="x")
        self.axial_step_label = tk.Label(self.axial_step_frame, text="Axial Step: ")
        self.axial_step_label.pack(side="left")
        self.axial_step_value = tk.IntVar()
        self.axial_step_entry = tk.Entry(self.axial_step_frame, textvariable=self.axial_step_value, width=4)
        self.axial_step_entry.pack(side="left")
        self.axial_step_value.set(self.ehooke.parameters.cellprocessingparams.axial_step)

        self.force_merge_below_frame = tk.Frame(self.parameters_panel)
        self.force_merge_below_frame.pack(side="top", fill="x")
        self.force_merge_below_label = tk.Label(self.force_merge_below_frame, text="Force Merge If Area Below: ")
        self.force_merge_below_label.pack(side="left")
        self.force_merge_below_value = tk.IntVar()
        self.force_merge_below_entry = tk.Entry(self.force_merge_below_frame,
                                                textvariable=self.force_merge_below_value, width=4)
        self.force_merge_below_entry.pack(side="left")
        self.force_merge_below_value.set(self.ehooke.parameters.cellprocessingparams.cell_force_merge_below)

        self.merge_dividing_frame = tk.Frame(self.parameters_panel)
        self.merge_dividing_frame.pack(side="top", fill="x")
        self.merge_dividing_label = tk.Label(self.merge_dividing_frame, text="Merge Dividing Cells: ")
        self.merge_dividing_label.pack(side="left")
        self.merge_dividing_value = tk.BooleanVar()
        self.merge_dividing_checkbox = tk.Checkbutton(self.merge_dividing_frame, variable=self.merge_dividing_value,
                                                      onvalue=True, offvalue=False)
        self.merge_dividing_checkbox.pack(side="left")
        self.merge_dividing_value.set(self.ehooke.parameters.cellprocessingparams.merge_dividing_cells)

        self.merge_length_tolerance_frame = tk.Frame(self.parameters_panel)
        self.merge_length_tolerance_frame.pack(side="top", fill="x")
        self.merge_length_tolerance_label = tk.Label(self.merge_length_tolerance_frame,
                                                     text="Length Tolerance on Merge: ")
        self.merge_length_tolerance_label.pack(side="left")
        self.merge_length_tolerance_value = tk.DoubleVar()
        self.merge_length_tolerance_entry = tk.Entry(self.merge_length_tolerance_frame,
                                                     textvariable=self.merge_length_tolerance_value, width=4)
        self.merge_length_tolerance_entry.pack(side="left")
        self.merge_length_tolerance_value.set(self.ehooke.parameters.cellprocessingparams.merge_length_tolerance)

        self.merge_min_interface_frame = tk.Frame(self.parameters_panel)
        self.merge_min_interface_frame.pack(side="top", fill="x")
        self.merge_min_interface_label = tk.Label(self.merge_min_interface_frame, text="Min Interface for Merge: ")
        self.merge_min_interface_label.pack(side="left")
        self.merge_min_interface_value = tk.IntVar()
        self.merge_min_interface_entry = tk.Entry(self.merge_min_interface_frame,
                                                  textvariable=self.merge_min_interface_value, width=4)
        self.merge_min_interface_entry.pack(side="left")
        self.merge_min_interface_value.set(self.ehooke.parameters.cellprocessingparams.merge_min_interface)

        self.cellprocessing_default_button = tk.Button(self.parameters_panel, text="Default Parameters",
                                                       command=self.load_default_params_cell_computation)
        self.cellprocessing_default_button.pack(side="top", fill="x")

        self.force_merge_button = tk.Button(self.parameters_panel, text="Force Merge",
                                            command=self.force_merge)
        self.force_merge_button.pack(side="top", fill="x")
        self.force_merge_button.config(state="disabled")

        self.split_cell_button = tk.Button(self.parameters_panel, text="Split Cell", command=self.split_cell)
        self.split_cell_button.pack(side="top", fill="x")
        self.split_cell_button.config(state="disabled")

        self.declare_as_noise_button = tk.Button(self.parameters_panel, text="Define as Noise",
                                                 command=self.declare_as_noise)
        self.declare_as_noise_button.config(state="disabled")
        self.declare_as_noise_button.pack(side="top", fill="x")

        self.undo_as_noise_button = tk.Button(self.parameters_panel, text="Undo as Noise",
                                                 command=self.undo_as_noise)
        self.undo_as_noise_button.config(state="disabled")
        self.undo_as_noise_button.pack(side="top", fill="x")

        self.cell_info_frame = tk.Frame(self.images_frame)
        self.cell_info_frame.pack(side="bottom", fill="x")

        self.cellid_frame = tk.Frame(self.cell_info_frame)
        self.cellid_frame.pack(side="top", fill="x")
        self.cellid_label = tk.Label(self.cellid_frame, text="Cell ID: ")
        self.cellid_label.pack(side="left")
        self.cellid_value = tk.IntVar()
        self.cellid_value_label = tk.Label(self.cellid_frame, textvariable=self.cellid_value)
        self.cellid_value_label.pack(side="left")

        self.merged_with_frame = tk.Frame(self.cell_info_frame)
        self.merged_with_frame.pack(side="top", fill="x")
        self.merged_with_label = tk.Label(self.merged_with_frame, text="Merged Cell: ")
        self.merged_with_label.pack(side="left")
        self.merged_with_value = tk.StringVar()
        self.merged_with_value_label = tk.Label(self.merged_with_frame, textvariable=self.merged_with_value)
        self.merged_with_value_label.pack(side="left")

        self.marked_as_noise_frame = tk.Frame(self.cell_info_frame)
        self.marked_as_noise_frame.pack(side="top", fill="x")
        self.marked_as_noise_label = tk.Label(self.marked_as_noise_frame, text="Marked as Noise: ")
        self.marked_as_noise_label.pack(side="left")
        self.marked_as_noise_value = tk.StringVar()
        self.marked_as_noise_value_label = tk.Label(self.marked_as_noise_frame, textvariable=self.marked_as_noise_value)
        self.marked_as_noise_value_label.pack(side="left")

        self.area_frame = tk.Frame(self.cell_info_frame)
        self.area_frame.pack(side="top", fill="x")
        self.area_label = tk.Label(self.area_frame, text="Area: ")
        self.area_label.pack(side="left")
        self.area_value = tk.IntVar()
        self.area_value_label = tk.Label(self.area_frame, textvariable=self.area_value)
        self.area_value_label.pack(side="left")

        self.perimeter_frame = tk.Frame(self.cell_info_frame)
        self.perimeter_frame.pack(side="top", fill="x")
        self.perimeter_label = tk.Label(self.perimeter_frame, text="Perimeter: ")
        self.perimeter_label.pack(side="left")
        self.perimeter_value = tk.IntVar()
        self.perimeter_value_label = tk.Label(self.perimeter_frame, textvariable=self.perimeter_value)
        self.perimeter_value_label.pack(side="left")

        self.length_frame = tk.Frame(self.cell_info_frame)
        self.length_frame.pack(side="top", fill="x")
        self.length_label = tk.Label(self.length_frame, text="Length: ")
        self.length_label.pack(side="left")
        self.length_value = tk.IntVar()
        self.length_value_label = tk.Label(self.length_frame, textvariable=self.length_value)
        self.length_value_label.pack(side="left")

        self.width_frame = tk.Frame(self.cell_info_frame)
        self.width_frame.pack(side="top", fill="x")
        self.width_label = tk.Label(self.width_frame, text="Width: ")
        self.width_label.pack(side="left")
        self.width_value = tk.IntVar()
        self.width_value_label = tk.Label(self.width_frame, textvariable=self.width_value)
        self.width_value_label.pack(side="left")

        self.eccentricity_frame = tk.Frame(self.cell_info_frame)
        self.eccentricity_frame.pack(side="top", fill="x")
        self.eccentricity_label = tk.Label(self.eccentricity_frame, text="Eccentricity: ")
        self.eccentricity_label.pack(side="left")
        self.eccentricity_value = tk.IntVar()
        self.eccentricity_value_label = tk.Label(self.eccentricity_frame, textvariable=self.eccentricity_value)
        self.eccentricity_value_label.pack(side="left")

        self.irregularity_frame = tk.Frame(self.cell_info_frame)
        self.irregularity_frame.pack(side="top", fill="x")
        self.irregularity_label = tk.Label(self.irregularity_frame, text="Irregularity: ")
        self.irregularity_label.pack(side="left")
        self.irregularity_value = tk.IntVar()
        self.irregularity_value_label = tk.Label(self.irregularity_frame, textvariable=self.irregularity_value)
        self.irregularity_value_label.pack(side="left")

        self.neighbours_frame = tk.Frame(self.cell_info_frame)
        self.neighbours_frame.pack(side="top", fill="x")
        self.neighbours_label = tk.Label(self.neighbours_frame, text="Neighbours: ")
        self.neighbours_label.pack(side="left")
        self.neighbours_value = tk.IntVar()
        self.neighbours_value_label = tk.Label(self.neighbours_frame, textvariable=self.neighbours_value)
        self.neighbours_value_label.pack(side="left")

        self.empty_label = tk.Label(self.cell_info_frame)
        self.empty_label.pack(side="top")

        self.base_button = tk.Button(self.images_frame, text="Base", command=lambda: self.show_image("Base"),
                                      width=self.image_buttons_width)
        self.base_button.pack(side="top", fill="x")
        self.base_button.config(state="active")

        self.base_features_button = tk.Button(self.images_frame, text="Base with Features",
                                               command=lambda: self.show_image("Base_features"),
                                               width=self.image_buttons_width)
        self.base_features_button.pack(side="top", fill="x")

        self.base_w_cells_button = tk.Button(self.images_frame, text="Base With Cells Outlined", command=lambda: self.show_image("Base_cells_outlined"),
                                      width=self.image_buttons_width)
        self.base_w_cells_button.pack(side="top", fill="x")
        self.base_w_cells_button.config(state="disabled")

        self.fluor_button = tk.Button(self.images_frame, text="Fluorescence", command=lambda: self.show_image("Fluorescence"),
                                      width=self.image_buttons_width)
        self.fluor_button.pack(side="top", fill="x")

        self.fluor_cells_out_button = tk.Button(self.images_frame, text="Fluorescence with Cells Outlined",
                                                command=lambda: self.show_image("Fluor_cells_outlined"),
                                                width=self.image_buttons_width)
        self.fluor_cells_out_button.pack(side="top", fill="x")
        self.fluor_cells_out_button.config(state="disabled")

    def set_cellcomputation_from_cellprocessing(self):
        """Method to go back to cell computation"""
        self.canvas.mpl_disconnect(self.cid)
        self.event_connected = False

        self.set_cellcomputation()

    def on_press(self, event):
        if event.button == 3:

            label = str(int(self.ehooke.cell_manager.merged_labels[int(event.ydata), int(event.xdata)]))

            if int(label) > 0:

                self.ehooke.cell_manager.cells[label].selection_state *= -1
                self.ehooke.cell_manager.overlay_cells(self.ehooke.image_manager)

                self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells

                self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

                self.show_image(self.current_image)

    def show_cell_info_cellprocessing(self, x, y):
        """Shows the stats of each cell (including fluor stats) on the side
        panel"""
        """Shows the stats of each cell on the side panel"""
        label = int(self.ehooke.cell_manager.merged_labels[int(y), int(x)])

        if 0 < label:
            stats = self.ehooke.cell_manager.cells[str(label)].stats

            self.cellid_value.set(label)
            self.merged_with_value.set(self.ehooke.cell_manager.cells[str(label)].merged_with)
            self.marked_as_noise_value.set(self.ehooke.cell_manager.cells[str(label)].marked_as_noise)
            self.area_value.set(int(stats["Area"]))
            self.perimeter_value.set(int(stats["Perimeter"]))
            self.length_value.set(int(stats["Length"]))
            self.width_value.set(int(stats["Width"]))
            self.eccentricity_value.set(float(str(stats["Eccentricity"])[0:6]))
            self.irregularity_value.set(float(str(stats["Irregularity"])[0:6]))
            self.neighbours_value.set(stats["Neighbours"])
            self.baseline_value.set(float(str(stats["Baseline"])[0:6]))
            self.cellmedian_value.set(float(str(stats["Cell Median"])[0:6]))
            self.permedian_value.set(float(str(stats["Membrane Median"])[0:6]))
            self.septmedian_value.set(float(str(stats["Septum Median"])[0:6]))
            self.cytomedian_value.set(float(str(stats["Cytoplasm Median"])[0:6]))
            self.fr_value.set(float(str(stats["Fluor Ratio"])[0:6]))
            self.fr75_value.set(float(str(stats["Fluor Ratio 75%"])[0:6]))
            self.fr25_value.set(float(str(stats["Fluor Ratio 25%"])[0:6]))
            self.fr10_value.set(float(str(stats["Fluor Ratio 10%"])[0:6]))

        else:
            self.cellid_value.set(0)
            self.merged_with_value.set("No")
            self.marked_as_noise_value.set("No")
            self.area_value.set(0)
            self.perimeter_value.set(0)
            self.length_value.set(0)
            self.width_value.set(0)
            self.eccentricity_value.set(0)
            self.irregularity_value.set(0)
            self.neighbours_value.set(0)
            self.baseline_value.set(0)
            self.cellmedian_value.set(0)
            self.permedian_value.set(0)
            self.septmedian_value.set(0)
            self.cytomedian_value.set(0)
            self.fr_value.set(0)
            self.fr75_value.set(0)
            self.fr25_value.set(0)
            self.fr10_value.set(0)

        return ""

    def process_cells(self):
        """Method used to process the individual regions of each cell
        aswell as their fluor stats"""
        self.ehooke.parameters.cellprocessingparams.find_septum = self.find_septum_checkbox_value.get()
        self.ehooke.parameters.cellprocessingparams.septum_algorithm = self.septum_algorithm_value.get()
        self.ehooke.parameters.cellprocessingparams.inner_mask_thickness = self.membrane_thickness_value.get()
        self.status.set("Processing cells...")
        self.ehooke.process_cells()

        self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells

        self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

        self.show_image(self.current_image)

        self.filter_cells_button.config(state="active")
        self.generate_report_button.config(state="active")
	self.select_all_button.config(state="active")
	self.unselect_all_button.config(state="active")
        self.status.set("Cell Processing Finished")

    def select_all_cells(self):
        """Method used to mark all cells as selected"""
        self.ehooke.select_all_cells()

        self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells
        self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

        self.show_image(self.current_image)

    def reject_all_cells(self):
        """Method used to mark all cells as rejected"""
        self.ehooke.reject_all_cells()

        self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells
        self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

        self.show_image(self.current_image)

    def change_cell_stat(self, label_c1):
        """Method used to change the state of a cell"""
        self.ehooke.cell_manager.cells[str(label_c1)].selection_state *= -1

    def filter_cells(self):
        """Method used to filter the cells based on a set of params"""
        filters = []

        if self.areafilter_checkbox_value.get():
            filters.extend([("Area", self.areafilter_min_value.get(), self.areafilter_max_value.get())])

        if self.perimeterfilter_checkbox_value.get():
            filters.extend([("Perimeter", self.perimeterfilter_min_value.get(), self.perimeterfilter_max_value.get())])

        if self.eccentricityfilter_checkbox_value.get():
            filters.extend([("Eccentricity", self.eccentricityfilter_min_value.get(),
                           self.eccentricityfilter_max_value.get())])

        if self.irregularityfilter_checkbox_value.get():
            filters.extend([("Irregularity", self.irregularityfilter_min_value.get(),
                             self.irregularityfilter_max_value.get())])

        if self.neighboursfilter_checkbox_value.get():
            filters.extend([("Neighbours", self.neighboursfilter_min_value.get(),
                             self.neighboursfilter_max_value.get())])

        self.ehooke.parameters.cellprocessingparams.cell_filters = filters
        self.ehooke.filter_cells()

        self.images["Fluor_cells_outlined"] = self.ehooke.cell_manager.fluor_w_cells
        self.images["Base_cells_outlined"] = self.ehooke.cell_manager.base_w_cells

        self.show_image(self.current_image)

    def generate_report(self):
        """Method used to save a report with the cell stats"""
        self.ehooke.generate_reports()

    def set_cellprocessing(self):
        """Method used to change the interface to the Cell Processing
        Step"""
        self.ax.cla()
        self.ax.axis("off")
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1)
        self.show_image("Fluor_cells_outlined")
        self.ax.format_coord = self.show_cell_info_cellprocessing
        self.canvas.show()

        self.cid = self.canvas.mpl_connect('button_release_event', self.on_press)

        for w in self.top_frame.winfo_children():
            w.destroy()

        for w in self.parameters_panel.winfo_children():
            w.destroy()

        for w in self.images_frame.winfo_children():
            w.destroy()

        self.empty_space = tk.Label(self.images_frame, text="")
        self.empty_space.pack(side="top")

        self.status = tk.StringVar()
        self.status.set("Load Phase Image")
        self.status_bar = tk.Label(self.parameters_panel, textvariable=self.status, wraplength= self.status_length)
        self.status_bar.pack(side="bottom")
        self.status.set("Right click to change selection state of a cell")

        self.process_cells_button = tk.Button(self.top_frame, text="Process Cells", command=self.process_cells)
        self.process_cells_button.pack(side="left")

        self.filter_cells_button = tk.Button(self.top_frame, text="Apply Filters", command=self.filter_cells)
        self.filter_cells_button.pack(side="left")
        self.filter_cells_button.config(state="disabled")

        self.generate_report_button = tk.Button(self.top_frame, text="Save Report", command=self.generate_report)
        self.generate_report_button.pack(side="left")
        self.generate_report_button.config(state="disabled")

        self.new_analysis_button = tk.Button(self.top_frame, text="New Analysis", command=self.new_analysis)
        self.new_analysis_button.pack(side="right")

        self.back_button = tk.Button(self.top_frame, text="Back", command=self.set_cellcomputation_from_cellprocessing)
        self.back_button.pack(side="right")

        self.cellprocessing_label = tk.Label(self.parameters_panel, text="Cell Processing Parameters: ")
        self.cellprocessing_label.pack(side="top")

        self.find_septum_frame = tk.Frame(self.parameters_panel)
        self.find_septum_frame.pack(side="top", fill="x")
        self.find_septum_label = tk.Label(self.find_septum_frame, text="Find Septum: ")
        self.find_septum_label.pack(side="left")
        self.find_septum_checkbox_value = tk.BooleanVar()
        self.find_septum_checkbox = tk.Checkbutton(self.find_septum_frame, variable=self.find_septum_checkbox_value,
                                                  onvalue=True, offvalue=False)
        self.find_septum_checkbox_value.set(False)
        self.find_septum_checkbox.pack(side="left")

        self.septum_algorithm_frame = tk.Frame(self.parameters_panel)
        self.septum_algorithm_frame.pack(side="top", fill="x")
        self.septum_algorithm_label = tk.Label(self.septum_algorithm_frame, text="Septum Algorithm: ")
        self.septum_algorithm_label.pack(side="left")
        self.septum_algorithm_value = tk.StringVar()
        self.septum_algorithm_menu = tk.OptionMenu(self.septum_algorithm_frame, self.septum_algorithm_value,
                                                 'Box', 'Isodata')
        self.septum_algorithm_menu.pack(side="left")
        self.septum_algorithm_value.set("Isodata")

        self.membrane_thickness_frame = tk.Frame(self.parameters_panel)
        self.membrane_thickness_frame.pack(side="top", fill="x")
        self.membrane_thickness_label = tk.Label(self.membrane_thickness_frame, text="Membrane Thickness: ")
        self.membrane_thickness_label.pack(side="left")
        self.membrane_thickness_value = tk.IntVar()
        self.membrane_thickness_entry = tk.Entry(self.membrane_thickness_frame, textvariable=self.membrane_thickness_value, width=4)
        self.membrane_thickness_entry.pack(side="left")
        self.membrane_thickness_value.set(self.ehooke.parameters.cellprocessingparams.inner_mask_thickness)

        self.filters_label = tk.Label(self.parameters_panel, text="Cell Filters: ")
        self.filters_label.pack(side="top")

        self.areafilter_frame = tk.Frame(self.parameters_panel)
        self.areafilter_frame.pack(side="top", fill="x")
        self.areafilter_checkbox_value = tk.BooleanVar()
        self.areafilter_checkbox = tk.Checkbutton(self.areafilter_frame, variable=self.areafilter_checkbox_value,
                                                  onvalue=True, offvalue=False)
        self.areafilter_checkbox_value.set(False)
        self.areafilter_checkbox.pack(side="left")
        self.areafilter_label = tk.Label(self.areafilter_frame, text="Area: ")
        self.areafilter_label.pack(side="left")
        self.areafilter_min_value = tk.IntVar()
        self.areafilter_max_value = tk.IntVar()
        self.areafilter_min_entry = tk.Entry(self.areafilter_frame, textvariable=self.areafilter_min_value, width=5)
        self.areafilter_min_entry.pack(side="left")
        self.areafilter_max_entry = tk.Entry(self.areafilter_frame, textvariable=self.areafilter_max_value, width=5)
        self.areafilter_max_entry.pack(side="left")
        self.areafilter_min_value.set(0)
        self.areafilter_max_value.set(1000)

        self.perimeterfilter_frame = tk.Frame(self.parameters_panel)
        self.perimeterfilter_frame.pack(side="top", fill="x")
        self.perimeterfilter_checkbox_value = tk.BooleanVar()
        self.perimeterfilter_checkbox = tk.Checkbutton(self.perimeterfilter_frame,
                                                       variable=self.perimeterfilter_checkbox_value,
                                                       onvalue=True, offvalue=False)
        self.perimeterfilter_checkbox_value.set(False)
        self.perimeterfilter_checkbox.pack(side="left")
        self.perimeterfilter_label = tk.Label(self.perimeterfilter_frame, text="Perimeter: ")
        self.perimeterfilter_label.pack(side="left")
        self.perimeterfilter_min_value = tk.IntVar()
        self.perimeterfilter_max_value = tk.IntVar()
        self.perimeterfilter_min_entry = tk.Entry(self.perimeterfilter_frame,
                                                  textvariable=self.perimeterfilter_min_value, width=5)
        self.perimeterfilter_min_entry.pack(side="left")
        self.perimeterfilter_max_entry = tk.Entry(self.perimeterfilter_frame,
                                                  textvariable=self.perimeterfilter_max_value, width=5)
        self.perimeterfilter_max_entry.pack(side="left")
        self.perimeterfilter_min_value.set(0)
        self.perimeterfilter_max_value.set(500)

        self.eccentricityfilter_frame = tk.Frame(self.parameters_panel)
        self.eccentricityfilter_frame.pack(side="top", fill="x")
        self.eccentricityfilter_checkbox_value = tk.BooleanVar()
        self.eccentricityfilter_checkbox = tk.Checkbutton(self.eccentricityfilter_frame,
                                                          variable=self.eccentricityfilter_checkbox_value,
                                                          onvalue=True, offvalue=False)
        self.eccentricityfilter_checkbox_value.set(False)
        self.eccentricityfilter_checkbox.pack(side="left")
        self.eccentricityfilter_label = tk.Label(self.eccentricityfilter_frame, text="Eccentricity: ")
        self.eccentricityfilter_label.pack(side="left")
        self.eccentricityfilter_min_value = tk.DoubleVar()
        self.eccentricityfilter_max_value = tk.DoubleVar()
        self.eccentricityfilter_min_entry = tk.Entry(self.eccentricityfilter_frame,
                                                     textvariable=self.eccentricityfilter_min_value, width=5)
        self.eccentricityfilter_min_entry.pack(side="left")
        self.eccentricityfilter_max_entry = tk.Entry(self.eccentricityfilter_frame,
                                                     textvariable=self.eccentricityfilter_max_value, width=5)
        self.eccentricityfilter_max_entry.pack(side="left")
        self.eccentricityfilter_min_value.set(-10)
        self.eccentricityfilter_max_value.set(10)

        self.irregularityfilter_frame = tk.Frame(self.parameters_panel)
        self.irregularityfilter_frame.pack(side="top", fill="x")
        self.irregularityfilter_checkbox_value = tk.BooleanVar()
        self.irregularityfilter_checkbox = tk.Checkbutton(self.irregularityfilter_frame,
                                                          variable=self.irregularityfilter_checkbox_value,
                                                  onvalue=True, offvalue=False)
        self.irregularityfilter_checkbox_value.set(False)
        self.irregularityfilter_checkbox.pack(side="left")
        self.irregularityfilter_label = tk.Label(self.irregularityfilter_frame, text="Irregularity: ")
        self.irregularityfilter_label.pack(side="left")
        self.irregularityfilter_min_value = tk.DoubleVar()
        self.irregularityfilter_max_value = tk.DoubleVar()
        self.irregularityfilter_min_entry = tk.Entry(self.irregularityfilter_frame,
                                                     textvariable=self.irregularityfilter_min_value, width=5)
        self.irregularityfilter_min_entry.pack(side="left")
        self.irregularityfilter_max_entry = tk.Entry(self.irregularityfilter_frame,
                                                     textvariable=self.irregularityfilter_max_value, width=5)
        self.irregularityfilter_max_entry.pack(side="left")
        self.irregularityfilter_min_value.set(0)
        self.irregularityfilter_max_value.set(20)

        self.neighboursfilter_frame = tk.Frame(self.parameters_panel)
        self.neighboursfilter_frame.pack(side="top", fill="x")
        self.neighboursfilter_checkbox_value = tk.BooleanVar()
        self.neighboursfilter_checkbox = tk.Checkbutton(self.neighboursfilter_frame,
                                                        variable=self.neighboursfilter_checkbox_value,
                                                  onvalue=True, offvalue=False)
        self.neighboursfilter_checkbox_value.set(False)
        self.neighboursfilter_checkbox.pack(side="left")
        self.neighboursfilter_label = tk.Label(self.neighboursfilter_frame, text="Neighbours: ")
        self.neighboursfilter_label.pack(side="left")
        self.neighboursfilter_min_value = tk.IntVar()
        self.neighboursfilter_max_value = tk.IntVar()
        self.neighboursfilter_min_entry = tk.Entry(self.neighboursfilter_frame,
                                                   textvariable=self.neighboursfilter_min_value, width=5)
        self.neighboursfilter_min_entry.pack(side="left")
        self.neighboursfilter_max_entry = tk.Entry(self.neighboursfilter_frame,
                                                   textvariable=self.neighboursfilter_max_value, width=5)
        self.neighboursfilter_max_entry.pack(side="left")
        self.neighboursfilter_min_value.set(0)
        self.neighboursfilter_max_value.set(10)

        self.select_all_button = tk.Button(self.parameters_panel, text="Select All Cells",
                                           command=self.select_all_cells)
        self.select_all_button.pack(side="top", fill="x")
	self.select_all_button.config(state="disabled")

        self.unselect_all_button = tk.Button(self.parameters_panel, text="Reject All Cells",
                                             command=self.reject_all_cells)
        self.unselect_all_button.pack(side="top", fill="x")
	self.unselect_all_button.config(state="disabled")

        self.cellprocessing_default_button = tk.Button(self.parameters_panel, text="Default Parameters",
                                                       command=self.load_default_params_cell_computation)
        self.cellprocessing_default_button.pack(side="top", fill="x")

        self.save_parameters_button = tk.Button(self.parameters_panel, text="Save Parameters", command=self.save_parameters)
        self.save_parameters_button.pack(side="top", fill="x")

        self.cell_info_frame = tk.Frame(self.images_frame)
        self.cell_info_frame.pack(side="bottom", fill="x")

        self.cellid_frame = tk.Frame(self.cell_info_frame)
        self.cellid_frame.pack(side="top", fill="x")
        self.cellid_label = tk.Label(self.cellid_frame, text="Cell ID: ")
        self.cellid_label.pack(side="left")
        self.cellid_value = tk.IntVar()
        self.cellid_value_label = tk.Label(self.cellid_frame, textvariable=self.cellid_value)
        self.cellid_value_label.pack(side="left")

        self.merged_with_frame = tk.Frame(self.cell_info_frame)
        self.merged_with_frame.pack(side="top", fill="x")
        self.merged_with_label = tk.Label(self.merged_with_frame, text="Merged Cell: ")
        self.merged_with_label.pack(side="left")
        self.merged_with_value = tk.StringVar()
        self.merged_with_value_label = tk.Label(self.merged_with_frame, textvariable=self.merged_with_value)
        self.merged_with_value_label.pack(side="left")

        self.marked_as_noise_frame = tk.Frame(self.cell_info_frame)
        self.marked_as_noise_frame.pack(side="top", fill="x")
        self.marked_as_noise_label = tk.Label(self.marked_as_noise_frame, text="Marked as Noise: ")
        self.marked_as_noise_label.pack(side="left")
        self.marked_as_noise_value = tk.StringVar()
        self.marked_as_noise_value_label = tk.Label(self.marked_as_noise_frame, textvariable=self.marked_as_noise_value)
        self.marked_as_noise_value_label.pack(side="left")

        self.area_frame = tk.Frame(self.cell_info_frame)
        self.area_frame.pack(side="top", fill="x")
        self.area_label = tk.Label(self.area_frame, text="Area: ")
        self.area_label.pack(side="left")
        self.area_value = tk.IntVar()
        self.area_value_label = tk.Label(self.area_frame, textvariable=self.area_value)
        self.area_value_label.pack(side="left")

        self.perimeter_frame = tk.Frame(self.cell_info_frame)
        self.perimeter_frame.pack(side="top", fill="x")
        self.perimeter_label = tk.Label(self.perimeter_frame, text="Perimeter: ")
        self.perimeter_label.pack(side="left")
        self.perimeter_value = tk.IntVar()
        self.perimeter_value_label = tk.Label(self.perimeter_frame, textvariable=self.perimeter_value)
        self.perimeter_value_label.pack(side="left")

        self.length_frame = tk.Frame(self.cell_info_frame)
        self.length_frame.pack(side="top", fill="x")
        self.length_label = tk.Label(self.length_frame, text="Length: ")
        self.length_label.pack(side="left")
        self.length_value = tk.IntVar()
        self.length_value_label = tk.Label(self.length_frame, textvariable=self.length_value)
        self.length_value_label.pack(side="left")

        self.width_frame = tk.Frame(self.cell_info_frame)
        self.width_frame.pack(side="top", fill="x")
        self.width_label = tk.Label(self.width_frame, text="Width: ")
        self.width_label.pack(side="left")
        self.width_value = tk.IntVar()
        self.width_value_label = tk.Label(self.width_frame, textvariable=self.width_value)
        self.width_value_label.pack(side="left")

        self.eccentricity_frame = tk.Frame(self.cell_info_frame)
        self.eccentricity_frame.pack(side="top", fill="x")
        self.eccentricity_label = tk.Label(self.eccentricity_frame, text="Eccentricity: ")
        self.eccentricity_label.pack(side="left")
        self.eccentricity_value = tk.IntVar()
        self.eccentricity_value_label = tk.Label(self.eccentricity_frame, textvariable=self.eccentricity_value)
        self.eccentricity_value_label.pack(side="left")

        self.irregularity_frame = tk.Frame(self.cell_info_frame)
        self.irregularity_frame.pack(side="top", fill="x")
        self.irregularity_label = tk.Label(self.irregularity_frame, text="Irregularity: ")
        self.irregularity_label.pack(side="left")
        self.irregularity_value = tk.IntVar()
        self.irregularity_value_label = tk.Label(self.irregularity_frame, textvariable=self.irregularity_value)
        self.irregularity_value_label.pack(side="left")

        self.neighbours_frame = tk.Frame(self.cell_info_frame)
        self.neighbours_frame.pack(side="top", fill="x")
        self.neighbours_label = tk.Label(self.neighbours_frame, text="Neighbours: ")
        self.neighbours_label.pack(side="left")
        self.neighbours_value = tk.IntVar()
        self.neighbours_value_label = tk.Label(self.neighbours_frame, textvariable=self.neighbours_value)
        self.neighbours_value_label.pack(side="left")

        self.baseline_frame = tk.Frame(self.cell_info_frame)
        self.baseline_frame.pack(side="top", fill="x")
        self.baseline_label = tk.Label(self.baseline_frame, text="Baseline: ")
        self.baseline_label.pack(side="left")
        self.baseline_value = tk.IntVar()
        self.baseline_value_label = tk.Label(self.baseline_frame, textvariable=self.baseline_value)
        self.baseline_value_label.pack(side="left")

        self.cellmedian_frame = tk.Frame(self.cell_info_frame)
        self.cellmedian_frame.pack(side="top", fill="x")
        self.cellmedian_label = tk.Label(self.cellmedian_frame, text="Cell Median: ")
        self.cellmedian_label.pack(side="left")
        self.cellmedian_value = tk.IntVar()
        self.cellmedian_value_label = tk.Label(self.cellmedian_frame, textvariable=self.cellmedian_value)
        self.cellmedian_value_label.pack(side="left")

        self.permedian_frame = tk.Frame(self.cell_info_frame)
        self.permedian_frame.pack(side="top", fill="x")
        self.permedian_label = tk.Label(self.permedian_frame, text="Perimeter Median: ")
        self.permedian_label.pack(side="left")
        self.permedian_value = tk.IntVar()
        self.permedian_value_label = tk.Label(self.permedian_frame, textvariable=self.permedian_value)
        self.permedian_value_label.pack(side="left")

        self.septmedian_frame = tk.Frame(self.cell_info_frame)
        self.septmedian_frame.pack(side="top", fill="x")
        self.septmedian_label = tk.Label(self.septmedian_frame, text="Septum Median: ")
        self.septmedian_label.pack(side="left")
        self.septmedian_value = tk.IntVar()
        self.septmedian_value_label = tk.Label(self.septmedian_frame, textvariable=self.septmedian_value)
        self.septmedian_value_label.pack(side="left")

        self.cytomedian_frame = tk.Frame(self.cell_info_frame)
        self.cytomedian_frame.pack(side="top", fill="x")
        self.cytomedian_label = tk.Label(self.cytomedian_frame, text="Cyto Median: ")
        self.cytomedian_label.pack(side="left")
        self.cytomedian_value = tk.IntVar()
        self.cytomedian_value_label = tk.Label(self.cytomedian_frame, textvariable=self.cytomedian_value)
        self.cytomedian_value_label.pack(side="left")

        self.fr_frame = tk.Frame(self.cell_info_frame)
        self.fr_frame.pack(side="top", fill="x")
        self.fr_label = tk.Label(self.fr_frame, text="FR: ")
        self.fr_label.pack(side="left")
        self.fr_value = tk.IntVar()
        self.fr_value_label = tk.Label(self.fr_frame, textvariable=self.fr_value)
        self.fr_value_label.pack(side="left")

        self.fr75_frame = tk.Frame(self.cell_info_frame)
        self.fr75_frame.pack(side="top", fill="x")
        self.fr75_label = tk.Label(self.fr75_frame, text="FR 75%: ")
        self.fr75_label.pack(side="left")
        self.fr75_value = tk.IntVar()
        self.fr75_value_label = tk.Label(self.fr75_frame, textvariable=self.fr75_value)
        self.fr75_value_label.pack(side="left")

        self.fr25_frame = tk.Frame(self.cell_info_frame)
        self.fr25_frame.pack(side="top", fill="x")
        self.fr25_label = tk.Label(self.fr25_frame, text="FR 25%: ")
        self.fr25_label.pack(side="left")
        self.fr25_value = tk.IntVar()
        self.fr25_value_label = tk.Label(self.fr25_frame, textvariable=self.fr25_value)
        self.fr25_value_label.pack(side="left")

        self.fr10_frame = tk.Frame(self.cell_info_frame)
        self.fr10_frame.pack(side="top", fill="x")
        self.fr10_label = tk.Label(self.fr10_frame, text="FR 10%: ")
        self.fr10_label.pack(side="left")
        self.fr10_value = tk.IntVar()
        self.fr10_value_label = tk.Label(self.fr10_frame, textvariable=self.fr10_value)
        self.fr10_value_label.pack(side="left")

        self.empty_label = tk.Label(self.cell_info_frame)
        self.empty_label.pack(side="top")

        self.base_button = tk.Button(self.images_frame, text="Base", command=lambda: self.show_image("Base"),
                                      width=self.image_buttons_width)
        self.base_button.pack(side="top", fill="x")
        self.base_button.config(state="active")

        self.base_features_button = tk.Button(self.images_frame, text="Base with Features",
                                               command=lambda: self.show_image("Base_features"),
                                               width=self.image_buttons_width)
        self.base_features_button.pack(side="top", fill="x")

        self.base_w_cells_button = tk.Button(self.images_frame, text="Base With Cells Outlined", command=lambda: self.show_image("Base_cells_outlined"),
                                      width=self.image_buttons_width)
        self.base_w_cells_button.pack(side="top", fill="x")
        self.base_w_cells_button.config(state="active")

        self.fluor_button = tk.Button(self.images_frame, text="Fluorescence", command=lambda: self.show_image("Fluorescence"),
                                      width=self.image_buttons_width)
        self.fluor_button.pack(side="top", fill="x")

        self.fluor_cells_out_button = tk.Button(self.images_frame, text="Fluorescence with Cells Outlined",
                                                command=lambda: self.show_image("Fluor_cells_outlined"),
                                                width=self.image_buttons_width)
        self.fluor_cells_out_button.pack(side="top", fill="x")
        self.fluor_cells_out_button.config(state="active")

    def new_analysis(self):
        """Restarts ehooke to conduct a new analysis"""
        self.ehooke = EHooke()
        self.default_params = self.ehooke.parameters
        self.images = {}
        self.current_image = None
        self.set_imageloader()

    def on_closing(self):
        """Creates a prompt when trying to close the main windows"""
        if tkMessageBox.askokcancel("Quit", "Do you want to quit?"):
            self.main_window.destroy()

if __name__ == "__main__":
    interface = Interface()
    interface.main_window.protocol("WM_DELETE_WINDOW", interface.on_closing)
    interface.main_window.mainloop()
