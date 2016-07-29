import Tkinter as tk
import os
import tkFileDialog as FD
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from skimage.io import imread, imsave
from skimage.color import rgb2gray
from skimage.util import img_as_float


class Interface(object):
    def __init__(self):
        self.identifier = Identifier()
        self.current_index = 0
        self.cells_keys = None

        self.main_window = tk.Tk()
        self.main_window.wm_title("cyphID")

        self.top_frame = tk.Frame(self.main_window)
        self.top_frame.pack(fill="x")

        self.middle_frame = tk.Frame(self.main_window)
        self.middle_frame.pack(fill="x")

        self.bottom_frame = tk.Frame(self.main_window)
        self.bottom_frame.pack(fill="x")

        self.askdirectory_button = tk.Button(self.top_frame, text="Choose Images Directory",
                                             command=self.choose_report_directory)
        self.askdirectory_button.pack()

        self.increase_figsize_button = tk.Button(self.top_frame, text=" + ", command=self.increase_figsize)
        self.increase_figsize_button.pack(side="right")

        self.decrease_figsize_button = tk.Button(self.top_frame, text="  -  ", command=self.decrease_figsize)
        self.decrease_figsize_button.pack(side="right")

        self.label_text = tk.StringVar()
        self.label_text.set("No Cells Selected")
        self.number_of_cells_label = tk.Label(
            self.top_frame, textvariable=self.label_text)
        self.number_of_cells_label.pack()

        # creates the figure canvas
        self.fig_width = 6
        self.fig_height = 3
        self.fig = plt.figure(figsize=(self.fig_width, self.fig_height), frameon=True)
        self.canvas = FigureCanvasTkAgg(self.fig, self.middle_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(fill="both")

        self.ax = plt.subplot(111)
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1)
        self.ax.axis("off")
        plt.autoscale(False)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.middle_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(fill="both")

        self.ax.format_coord = self.show_nothing

        self.canvas.show()

        self.phase1_button = tk.Button(
            self.bottom_frame, text="Phase 1", command=lambda: self.choose_phase(1))
        self.phase1_button.pack(side="left")

        self.phase2_button = tk.Button(
            self.bottom_frame, text="Phase 2", command=lambda: self.choose_phase(2))
        self.phase2_button.pack(side="left")

        self.phase3_button = tk.Button(
            self.bottom_frame, text="Phase 3", command=lambda: self.choose_phase(3))
        self.phase3_button.pack(side="left")

        self.discard_button = tk.Button(
            self.bottom_frame, text="Discard Cell", command=lambda: self.choose_phase("Discard"))
        self.discard_button.pack(side="left")

        self.back_button = tk.Button(self.bottom_frame, text="Back", command=self.previous_cell)
        self.back_button.pack(side="right")

        self.luminance = tk.StringVar()
        self.luminance.set("0")
        self.luminance_label = tk.Label(
            self.bottom_frame, textvariable=self.luminance)
        self.luminance_label.pack(side="right")

    def key(self, event):
        if event.char == "1":
            self.choose_phase(1)
        elif event.char == "2":
            self.choose_phase(2)
        elif event.char == "3":
            self.choose_phase(3)
        elif event.char == "d":
            self.choose_phase("Discard")
        elif event.char == "b":
            self.previous_cell()

    def show_nothing(self, x, y):

        return ""

    def show_luminance(self, x, y):
        current_image = self.identifier.cells[self.cells_keys[self.current_index]].image
        return "Luminance: " + str(rgb2gray(current_image)[int(y), int(x)])

    def choose_report_directory(self):
        self.current_index = 0
        self.identifier.load_cells()
        self.cells_keys = self.identifier.cells.keys()
        self.total_cells = str(len(self.cells_keys))
        self.label_text.set(str(self.current_index+1) + "of" + self.total_cells + " total")
        self.main_window.bind("<Key>", self.key)
        self.show_image(self.current_index)

    def show_image(self, index):
        self.ax.cla()
        self.ax.imshow(self.identifier.cells[self.cells_keys[self.current_index]].image)
        self.ax.format_coord = self.show_luminance
        self.label_text.set(str(self.current_index+1) + " of " + self.total_cells + " total")
        self.canvas.show()

    def choose_phase(self, phase):
        self.identifier.choose_phase(self.cells_keys[self.current_index], str(phase))

        if self.current_index >= len(self.cells_keys)-1:
            self.identifier.save_report()
            self.plot_stats()
        else:
            self.current_index += 1
            self.show_image(self.current_index)

    def previous_cell(self):
        if self.current_index >= len(self.cells_keys)-1:
            self.identifier.discarded_count = 0
            self.identifier.phase1_count = 0
            self.identifier.phase2_count = 0
            self.identifier.phase3_count = 0
            self.identifier.phase1_ids = ""
            self.identifier.phase2_ids = ""
            self.identifier.phase3_ids = ""

        self.current_index -= 1
        if self.current_index < 0:
            self.current_index = 0
        self.show_image(self.current_index)

    def increase_figsize(self):

        self.fig_width += 1
        self.fig_height += 1

        self.fig.set_size_inches(self.fig_width, self.fig_height, forward=True)

        self.show_image(self.current_index)

    def decrease_figsize(self):
        if self.fig_height > 2:
            self.fig_width -= 1
            self.fig_height -= 1
        self.fig.set_size_inches(self.fig_width, self.fig_height, forward=True)

        self.show_image(self.current_index)

    def plot_stats(self):
        self.label_text.set("All Images Done")
        plt.clf()
        self.ax = plt.subplot(111)
        self.ax.format_coord = self.show_nothing
        plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8)
        x = np.arange(4)
        width = 0.35
        self.ax.set_ylabel("Cell Count")
        self.ax.set_xlabel("Cell Cycle Phase")
        self.ax.set_xticks(x +width)
        self.ax.set_xticklabels(("Discarded", "Phase 1", "Phase 2", "Phase 3"))
        self.ax.bar(x, [self.identifier.discarded_count, self.identifier.phase1_count, self.identifier.phase2_count, self.identifier.phase3_count], width)
        self.ax.set_title("Phase 1:" + str(self.identifier.phase1_count) + " - " + "Phase 2: " + str(self.identifier.phase2_count) + " - " + "Phase 3: " + str(self.identifier.phase3_count))
        self.canvas.show()

    def new_analysis(self):
        self.identifier = Identifier()
        self.current_index = 0
        self.cells_keys = None
        self.ax.cla()
        self.canvas.show()


class Cell(object):
    def __init__(self):
        self.image = None
        self.id = None
        self.phase = None


class Identifier(object):
    def __init__(self):
        self.report_path = None
        self.path = None
        self.cells = {}
        self.phase1_ids = ""
        self.phase2_ids = ""
        self.phase3_ids = ""
        self.phase1_count = 0
        self.phase2_count = 0
        self.phase3_count = 0
        self.discarded_count = 0

    def load_cells(self):
        self.phase1_ids = ""
        self.phase2_ids = ""
        self.phase3_ids = ""
        self.phase1_count = 0
        self.phase2_count = 0
        self.phase3_count = 0
        self.discarded_count = 0
        self.report_path = FD.askdirectory()
        self.path = self.report_path + "/_cyphid_images/"
        images_list = sorted(os.listdir(self.report_path + "/_images"))

        for image in images_list:
            id = str(image.split(".")[0])
            self.cells[id] = Cell()
            self.cells[id].id = id
            self.cells[id].image = imread(self.report_path + "/_images/"+image)
            self.cells[id].image = self.cells[id].image[:, len(self.cells[id].image[0]) / 5:len(self.cells[id].image[0]) * 2 / 5:]
            self.cells[id].image = img_as_float(self.cells[id].image)


    def choose_phase(self, cell_id, phase):
        self.cells[cell_id].phase = phase

    def save_report(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)

        for key in self.cells.keys():
            if self.cells[key].phase == "1":
                self.phase1_count += 1
                self.phase1_ids += str(self.cells[key].id) + ";"
                imsave(self.path + "1_" + self.cells[key].id + ".png", self.cells[key].image)
            elif self.cells[key].phase == "2":
                self.phase2_count += 1
                self.phase2_ids += str(self.cells[key].id) + ";"
                imsave(self.path + "2_" + self.cells[key].id + ".png", self.cells[key].image)
            elif self.cells[key].phase == "3":
                self.phase3_count += 1
                self.phase3_ids += str(self.cells[key].id) + ";"
                imsave(self.path + "/" + "3_" + self.cells[key].id + ".png", self.cells[key].image)
            else:
                self.discarded_count += 1

        report = ["Discarded Cells:;" + str(self.discarded_count) + "\n", "Phase 1 Cells:;" + str(self.phase1_count) + "\n", "Phase 2 Cells:;" + str(self.phase2_count) + "\n", "Phase 3 Cells:;" + str(self.phase3_count) + "\n"]
        open(self.report_path + "/cyphID_report.csv", 'w').writelines(report)

        open(self.report_path + "/phase1_cells.txt",
             "w").writelines(self.phase1_ids)
        open(self.report_path + "/phase2_cells.txt",
             "w").writelines(self.phase2_ids)
        open(self.report_path + "/phase3_cells.txt",
             "w").writelines(self.phase3_ids)

if __name__ == "__main__":
    identifier = Interface()
    identifier.main_window.mainloop()
