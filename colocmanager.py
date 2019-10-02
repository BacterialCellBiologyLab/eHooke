import os
import numpy as np
from tkinter import filedialog as fd
from scipy.stats import pearsonr

class ColocManager(object):

    def __int__(self):
        self.report = {}

    def save_report(self, sept=False):

        sorted_keys = sorted(self.report.keys())

        header = ["Whole Cell", "Membrane", "Cytoplasm"]
        if sept:
            header.extend(["Septum", "MembSept"])

        results = "Cell ID;"
        results += ";".join(header)
        results += ";\n"

        for key in sorted_keys:
            results += key + ";"
            for measurement in header:
                results += str(self.report[key][measurement]) + ";"

            results += "\n"

        save_directory = fd.askdirectory()
        open(save_directory + os.sep + "pcc_report.csv", "w").writelines(results)


    def pearsons_score(self, channel_1, channel_2, mask):

        filtered_1 = (channel_1 * mask).flatten()
        filtered_1 = filtered_1[filtered_1 > 0.0]
        filtered_2 = (channel_2 * mask).flatten()
        filtered_2 = filtered_2[filtered_2 > 0.0]

        return pearsonr(filtered_1, filtered_2)

    def compute_pcc(self, cell_manager, image_manager, parameters):
        self.report = {}

        fluor_image = image_manager.original_fluor_image
        optional_image = image_manager.optional_image

        for key in cell_manager.cells.keys():

            if cell_manager.cells[key].selection_state == 1:

                self.report[key] = {}

                cell = cell_manager.cells[key]

                x0, y0, x1, y1 = cell.box

                fluor_box = fluor_image[x0:x1+1, y0:y1+1]
                optional_box = optional_image[x0:x1+1, y0:y1+1]

                self.report[key]["Channel 1"] = fluor_box
                self.report[key]["Channel 2"] = optional_box

                self.report[key]["Whole Cell"] = self.pearsons_score(fluor_box, optional_box, cell.cell_mask)[0]
                self.report[key]["Membrane"] = self.pearsons_score(fluor_box, optional_box, cell.perim_mask)[0]
                self.report[key]["Cytoplasm"] = self.pearsons_score(fluor_box, optional_box, cell.cyto_mask)[0]

                if parameters.cellprocessingparams.find_septum:
                    self.report[key]["Septum"] = self.pearsons_score(fluor_box, optional_box, cell.sept_mask)[0]
                    self.report[key]["MembSept"] = self.pearsons_score(fluor_box, optional_box, cell.membsept_mask)[0]

        self.save_report(sept=parameters.cellprocessingparams.find_septum)

