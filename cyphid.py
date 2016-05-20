import Tkinter as tk
import os
import tkFileDialog as FD
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from skimage.io import imread, imsave
from skimage.util import img_as_int, img_as_float
from skimage import exposure
from skimage.color import rgb2gray


class Cyphid:

    def __init__(self):
        self.images_list = []
        self.phase1_ids = ""
        self.phase2_ids = ""
        self.phase3_ids = ""
        self.current_image = None
        self.current_image_index = 0
        self.phase1_count = 0
        self.phase2_count = 0
        self.phase3_count = 0
        self.discarded_count = 0

        # GUI
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

        self.label_text = tk.StringVar()
        self.label_text.set("No Cells Selected")
        self.number_of_cells_label = tk.Label(
            self.top_frame, textvariable=self.label_text)
        self.number_of_cells_label.pack()

        # creates the figure canvas
        self.fig = plt.figure(figsize=(8, 5), frameon=True)
        self.canvas = FigureCanvasTkAgg(self.fig, self.middle_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side="top")

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
        
        self.luminance = tk.StringVar()
        self.luminance.set("0")
        self.luminance_label = tk.Label(self.bottom_frame, textvariable=self.luminance)
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

    def new_analysis(self):
        self.images_list = []
        self.current_image = None
        self.current_image_index = 0
        self.phase1_count = 0
        self.phase2_count = 0
        self.phase3_count = 0
        self.discarded_count = 0

    def show_nothing(self, x, y):
        
        return ""

    def show_luminance(self, x, y):
        
        return "Luminance: " + str(rgb2gray(self.current_image)[int(y), int(x)])

    def choose_report_directory(self):
        self.new_analysis()
        self.report_path = FD.askdirectory()
        self.images_list = sorted(os.listdir(self.report_path + "/_images"))
        self.label_text.set(
            "1 Cell of " + str(len(self.images_list)) + " total.")

        self.path = self.report_path + "/cyphid_images"
        if not os.path.exists(self.path):
            os.makedirs(self.path)

        plt.clf()
        self.ax = plt.subplot(111)
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1)
        self.ax.axis("off")
        self.ax.format_coord = self.show_luminance

        self.show_image(0)
        self.main_window.bind("<Key>", self.key)

    def show_image(self, image):
        self.ax.cla()
        self.current_image = imread(
            self.report_path + "/_images/" + self.images_list[image])
        self.current_image = self.current_image[
            :, len(self.current_image[0]) / 5:len(self.current_image[0]) * 2 / 5:]
        
        self.ax.imshow(self.current_image)
        self.ax.format_coord = self.show_luminance
        self.canvas.show()

    def choose_phase(self, phase):
        if phase == 1:
            self.phase1_count += 1
            imsave(self.path + "/1_" +
                   self.images_list[self.current_image_index], self.current_image)
            self.phase1_ids += self.images_list[self.current_image_index].split(".")[0] + ";"
            
        elif phase == 2:
            self.phase2_count += 1
            imsave(self.path + "/2_" +
                   self.images_list[self.current_image_index], self.current_image)
            self.phase2_ids += self.images_list[self.current_image_index].split(".")[0] + ";"
             
        elif phase == 3:
            self.phase3_count += 1
            imsave(self.path + "/3_" +
                   self.images_list[self.current_image_index], self.current_image)
            self.phase3_ids += self.images_list[self.current_image_index].split(".")[0] + ";"

        else:
            self.discarded_count += 1

        if self.current_image_index == len(self.images_list) - 1:
            self.label_text.set("All Images Done")
            plt.clf()
            self.ax = plt.subplot(111)
            self.ax.format_coord = self.show_nothing
            plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8)
            x = np.arange(4)
            width = 0.35
            self.ax.set_ylabel("Cell Count")
            self.ax.set_xlabel("Cell Cycle Phase")
            self.ax.set_xticks(x + width)
            self.ax.set_xticklabels(
                ("Discarded", "Phase 1", "Phase 2", "Phase 3"))
            self.ax.bar(x, [self.discarded_count, self.phase1_count,
                            self.phase2_count, self.phase3_count], width)
            self.ax.set_title("Phase 1: " + str(self.phase1_count) + " - " + "Phase 2: " + str(self.phase2_count) + " - " + "Phase 3: " + str(self.phase3_count))
            #self.canvas.draw()
            self.canvas.show()

            report = ["Discarded Cells:;" + str(self.discarded_count) + "\n", "Phase 1 Cells:;" + str(self.phase1_count) + "\n",
                      "Phase 2 Cells:;" + str(self.phase2_count) + "\n", "Phase 3 Cells:;" + str(self.phase3_count) + "\n"]
            open(self.report_path + "/cyphID_report.csv", 'w').writelines(report)
            
            open(self.report_path + "/phase1_cells.txt", "w").writelines(self.phase1_ids)
            open(self.report_path + "/phase2_cells.txt", "w").writelines(self.phase2_ids)
            open(self.report_path + "/phase3_cells.txt", "w").writelines(self.phase3_ids)

        else:
            self.current_image_index += 1
            self.label_text.set(str(self.current_image_index + 1) +
                                " of " + str(len(self.images_list)) + " total.")

            self.show_image(self.current_image_index)

if __name__ == "__main__":
    cyphid = Cyphid()
    cyphid.main_window.mainloop()
