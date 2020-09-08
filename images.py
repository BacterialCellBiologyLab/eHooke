"""This module is responsible for the loading and processing of
both the base and fluor image.
The base image corresponds to the image that is going to be used for
the computation of the mask and the fluor image corresponds to the
image that is used to measure the fluorescence.
Contains a single class: the ImageManager class which is responsible for
the handling of the images and the computation of the mask.
This class should also be responsible for the connection of this module
with the main module of the software."""

from tkinter.filedialog import asksaveasfilename
import numpy as np
from copy import deepcopy
from skimage.segmentation import mark_boundaries
from skimage.io import imsave, imread
from skimage.util import img_as_float, img_as_uint
from skimage.filters import threshold_isodata, threshold_local
from skimage import exposure, color, morphology
from scipy import ndimage


class ImageManager(object):
    """Main class of the module. This class is responsible for the loading of
    the base image and the fluor image aswell as the computation of the masks.
    The class contains the load_base_image, set_clip, compute_masks,
    load_fluor_image and overlay_mask methods,
    that should take care of all the needed functions of this module."""

    def __init__(self):
        self.base_image = None
        self.clip = None
        self.base_mask = None
        self.mask = None
        self.fluor_image = None
        self.original_fluor_image = None
        self.optional_image = None
        self.base_w_mask = None
        self.fluor_w_mask = None
        self.optional_w_mask = None
        self.align_values = (0, 0)

    def clear_all(self):
        """Sets the class back to the __init__ state"""

        self.base_image = None
        self.clip = None
        self.base_mask = None
        self.mask = None
        self.fluor_image = None
        self.original_fluor_image = None
        self.optional_image = None
        self.base_w_mask = None
        self.fluor_w_mask = None
        self.optional_w_mask = None
        self.align_values = (0, 0)

    def set_clip(self, margin):
        """ Defines the clipping size based on the base image dimensions
        and the border margin. Stores the clipping size on self.clip in the
        (x0, y0, x1, y1) format
        In case the margin is less than 10px, this value is going to be set to
        10, to make sure there is room for the alignments of the fluor image"""

        if margin < 10:
            border = 10
        else:
            border = margin

        x_length, y_length = self.base_image.shape
        self.clip = (border, border, x_length - border,
                     y_length - border)

    def load_base_image(self, filename, params):
        """This method is responsible for the loading of the base image and
        its storage on the self.base_image of the class instance.
        Requires a filename (or path) of the image to be open (should be
        compatible with the imread function of the scikit-image.io module)
        and an instance of the ImageLoadingParams of the parameters module."""

        image = img_as_float(imread(filename))

        # note: changed order, rescale_intensity was called
        # before the rgb2gray, test to see which one is better
        image = color.rgb2gray(image)
        image = exposure.rescale_intensity(image)

        if params.invert_base:
            image = 1 - image

        self.base_image = image
        self.set_clip(params.border)

    def compute_base_mask(self, params):
        """Creates the base mask for the base image.
        Needs the base image, an instance of imageloaderparams
        and the clip area, which should be already defined
        by the load_base_image method
        To create the base mask, two algorithms are available, on based on the
        threshold_isodata and the other one on the threshold_local functions
        of the scikit-image.threshold module.
        """
        x0, y0, x1, y1 = self.clip
        base_mask = np.copy(self.base_image[x0:x1, y0:y1])

        if params.mask_algorithm == "Isodata":
            isodata_threshold = threshold_isodata(base_mask)
            base_mask = img_as_float(base_mask <= isodata_threshold)

        elif params.mask_algorithm == "Local Average":
            # need to invert because threshold_adaptive sets dark parts to 0
            block_size = params.mask_blocksize

            if block_size%2 == 0:
                block_size += 1

            threshold = threshold_local(base_mask,
                                              block_size,
                                              method="gaussian",
                                              offset=params.mask_offset)

            base_mask = 1.0 - (base_mask > threshold)
        
        elif params.mask_algorithm == "Absolute":
            value = float(raw_input("Insert Threshold Value: "))
            print(value)

            base_mask = img_as_float(base_mask <= value)

        else:
            print("Not a valid mask algorithm")

        self.base_mask = 1 - base_mask

    def compute_mask(self, params):
        """Creates the mask for the base image.
        Needs the base image, an instance of imageloaderparams
        and the clip area, which should be already defined
        by the load_base_image method.
        Creates the mask by improving the base mask created by the
        compute_base_mask method. Applies the mask closing, dilation and
        fill holes parameters.
        """
        self.compute_base_mask(params)

        mask = np.copy(self.base_mask)
        closing_matrix = np.ones((int(params.mask_closing), int(params.mask_closing)))

        if params.mask_closing > 0:
            # removes small dark spots and then small white spots
            mask = img_as_float(morphology.closing(
                mask, closing_matrix))
            mask = 1 - \
                img_as_float(morphology.closing(
                    1 - mask, closing_matrix))

        for f in range(params.mask_dilation):
            mask = morphology.erosion(mask, np.ones((3, 3)))

        if params.mask_fill_holes:
            # mask is inverted
            mask = 1 - img_as_float(ndimage.binary_fill_holes(1.0 - mask))

        self.mask = mask

        self.overlay_mask_base_image()

    def load_fluor_image(self, filename, params):
        """This method is responsible for the loading of the fluor image and
        its storage on the self.fluor_image of the class instance.
        Requires a filename (or path) of the image to be open (should be
        compatible with the imread function of the scikit-image.io module)
        and an instance of the ImageLoadingParams of the parameters module."""

        inverted_mask = 1 - self.mask

        fluor_image = imread(filename)

        if len(fluor_image.shape) > 2:
            fluor_image = color.rgb2gray(fluor_image)

        self.original_fluor_image = deepcopy(fluor_image)

        fluor_image = img_as_float(fluor_image)

        best = (0, 0)
        x0, y0, x1, y1 = self.clip

        if params.auto_align:
            minscore = 0
            width = params.border
            for dx in range(-width, width):
                for dy in range(-width, width):
                    tot = -np.sum(np.multiply(inverted_mask,
                                              fluor_image[x0 + dx:x1 + dx,
                                                          y0 + dy:y1 + dy]))
                                                          
                    if tot < minscore:
                        minscore = tot
                        best = (dx, dy)

        else:
            best = (params.x_align, params.y_align)

        self.align_values = best

        dx, dy = best
        self.original_fluor_image = self.original_fluor_image[x0 + dx:x1 + dx, y0 + dy:y1 + dy]
        self.fluor_image = fluor_image[x0 + dx:x1 + dx, y0 + dy:y1 + dy]

        self.overlay_mask_fluor_image()

    def load_option_image(self, filename, params):
        """Loads an option image that can be used to look for the septum
        and to help classify the cell cycle phases. No fluorescence is measured
        on this image"""

        inverted_mask = 1 - self.mask

        optional_image = imread(filename)

        if len(optional_image.shape) > 2:
            optional_image = color.rgb2gray(optional_image)

        optional_image = img_as_float(optional_image)

        best = (0, 0)
        x0, y0, x1, y1 = self.clip

        if params.auto_align:
            minscore = 0
            width = params.border
            for dx in range(-width, width):
                for dy in range(-width, width):
                    tot = -np.sum(np.multiply(inverted_mask,
                                              optional_image[x0 + dx:x1 + dx,
                                                             y0 + dy:y1 + dy]))
                                                          
                    if tot < minscore:
                        minscore = tot
                        best = (dx, dy)

        else:
            best = (params.x_align, params.y_align)

        dx, dy = best

        self.optional_image = optional_image[x0 + dx:x1 + dx, y0 + dy:y1 + dy]

    def overlay_mask_base_image(self):
        """ Creates a new image with an overlay of the mask
        over the base image"""

        x0, y0, x1, y1 = self.clip

        self.base_w_mask = mark_boundaries(self.base_image[x0:x1, y0:y1],
                                           img_as_uint(self.mask),
                                           color=(0, 1, 1),
                                           outline_color=None)

    def overlay_mask_fluor_image(self):
        """ Creates a new image with an overlay of the mask
        over the fluor image"""

        fluor_image = color.rgb2gray(self.fluor_image)
        fluor_image = exposure.rescale_intensity(fluor_image)
        fluor_image = img_as_float(fluor_image)

        self.fluor_w_mask = mark_boundaries(fluor_image, img_as_uint(self.mask),
                                            color=(0, 1, 1),
                                            outline_color=None)

    def overlay_mask_optional_image(self):
        """Creates a new image with an overlay of the mask over the fluor
        image"""

        optional_image = color.rgb2gray(self.optional_image)
        optional_image = exposure.rescale_intensity(optional_image)
        optional_image = img_as_float(optional_image)

        self.optional_w_mask = mark_boundaries(optional_image, img_as_uint(
            self.mask), color=(0, 1, 1), outline_color=None)

    def save_image(self, image_to_save, filename=None):
        """Saves the choosen image as a .png file.
        Can be called with a filname(path) or without one, in which case a
        the method calls a tkFileDialog asksaveasfilename window.
        The options for the image_to_save are: Base, Fluor, Base With Mask
        and Fluor With Mask"""
        if filename is None:
            filename = asksaveasfilename()

        if image_to_save == "Base":
            imsave(filename + ".png", self.base_image)
        elif image_to_save == "Fluor":
            imsave(filename + ".png", self.fluor_image)
        elif image_to_save == "Base With Mask":
            imsave(filename + ".png", self.base_w_mask)
        elif image_to_save == "Fluor With Mask":
            imsave(filename + ".png", self.fluor_w_mask)
        else:
            print("Not a valid image selection.")
            print("Choose between:")
            print("Base, Fluor, Base With Mask, Fluor With Mask")

    def save_mask(self, filename=None):

        if filename is None:
            filename = asksaveasfilename()

        imsave(filename + ".png", self.mask)

