"""Module used to create the report of the cell identification"""
from skimage.io import imsave
from skimage.util import img_as_float, img_as_uint
import cellprocessing as cp
import numpy as np
import os

CELL_SELECTED = 1
CELL_REJECTED = -1


class ReportManager:

    def __init__(self, parameters):
        self.keys = cp.stats_format(parameters.cellprocessingparams)

        self.cell_data_filename = None

    def csv_report(self, filename, image_name, cell_manager):

        cells = cell_manager.cells

        if len(cells) > 0:
            header = 'Cell ID '

            for k in self.keys:
                label, digits = k
                header = header + ';' + label
            selects = [header + '\n']
            rejects = [header + '\n']
            noise = [header + '\n']

            sorted_keys = []
            for k in sorted(cells.keys()):
                sorted_keys.append(int(k))

            sorted_keys = sorted(sorted_keys)

            for k in sorted_keys:
                cell = cells[str(k)]
                if cell.selection_state == CELL_SELECTED:
                    lin = str(int(cell.label))
                    for stat in self.keys:
                        lin = lin + ';' + str(cell.stats[stat[0]])
                    selects.append(lin + '\n')

                elif cell.selection_state == CELL_REJECTED:
                    lin = str(int(cell.label))

                    for stat in self.keys:
                        lin = lin + ";" + str(cell.stats[stat[0]])
                    rejects.append(lin + "\n")

                elif cell.selection_state == 0:
                    lin = str(int(cell.label))

                    for stat in self.keys:
                        lin = lin + ";" + str(cell.stats[stat[0]])
                    noise.append(lin + "\n")

            if len(selects) > 1:
                open(filename + "/csv_selected_" + image_name + ".csv", 'w').writelines(selects)

            if len(rejects) > 1:
                open(filename + "/csv_rejected_" + image_name + ".csv", "w").writelines(rejects)

            if len(noise) > 1:
                open(filename + "/csv_noise_" + image_name + ".csv", "w").writelines(noise)

    def html_report(self, filename, image_name, cell_manager):
        """generates an html report with the all the cell stats from the
        selected cells"""

        cells = cell_manager.cells

        HTML_HEADER = """<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN"
                        "http://www.w3.org/TR/html4/strict.dtd">
                    <html lang="en">
                      <head>
                        <meta http-equiv="content-type" content="text/html; charset=utf-8">
                        <title>title</title>
                        <link rel="stylesheet" type="text/css" href="style.css">
                        <script type="text/javascript" src="script.js"></script>
                      </head>
                      <body>\n"""

        report = [HTML_HEADER]

        if len(cells) > 0:
            header = '<table>\n<th>Cell ID</th><th>Images'
            for k in self.keys:
                label, digits = k
                header = header + '</th><th>' + label
            header += '</th>\n'
            selects = ['\n<h1>Selected cells:</h1>\n' + header + '\n']
            rejects = ['\n<h1>Rejected cells:</h1>\n' + header + '\n']
            noise = ['\n<h1>Noise:</h1>\n' + header + '\n']

            count = 0
            count2 = 0
            count3 = 0

            print("Total Cells: " + str(len(cells)))

            sorted_keys = []
            for k in sorted(cells.keys()):
                sorted_keys.append(int(k))

            sorted_keys = sorted(sorted_keys)

            for k in sorted_keys:
                cell = cells[str(k)]
                if cell.selection_state == CELL_SELECTED:
                    cellid = str(int(cell.label))
                    img = img_as_float(cell.image)
                    imsave(filename + "/_images" +
                           os.sep + cellid + '.png', img)
                    lin = '<tr><td>' + cellid + '</td><td><img src="./' + '_images/' + \
                          cellid + '.png" alt="pic" width="200"/></td>'

                    count += 1

                    for stat in self.keys:
                        lbl, digits = stat
                        lin = lin + '</td><td>' + \
                            ("{0:." + str(digits) +
                             "f}").format(cell.stats[lbl])

                    lin += '</td></tr>\n'
                    selects.append(lin)

                elif cell.selection_state == CELL_REJECTED:
                    cellid = str(int(cell.label))
                    img = img_as_float(cell.image)
                    imsave(filename + "/_rejected_images" +
                           os.sep + cellid + '.png', img)
                    lin = '<tr><td>' + cellid + '</td><td><img src="./' + '_rejected_images/' + \
                          cellid + '.png" alt="pic" width="200"/></td>'

                    count2 += 1

                    for stat in self.keys:
                        lbl, digits = stat
                        lin = lin + '</td><td>' + \
                            ("{0:." + str(digits) +
                             "f}").format(cell.stats[lbl])

                    lin += '</td></tr>\n'
                    rejects.append(lin)

                elif cell.selection_state == 0:
                    cellid = str(int(cell.label))
                    img = img_as_float(cell.image)
                    imsave(filename + "/_noise_images" +
                           os.sep + cellid + '.png', img)
                    lin = '<tr><td>' + cellid + '</td><td><img src="./' + '_noise_images/' + \
                          cellid + '.png" alt="pic" width="200"/></td>'

                    count3 += 1

                    for stat in self.keys:
                        lbl, digits = stat
                        lin = lin + '</td><td>' + \
                            ("{0:." + str(digits) +
                             "f}").format(cell.stats[lbl])

                    lin += '</td></tr>\n'
                    noise.append(lin)

            print("Selected Cells: " + str(count))
            print("Rejected Cells: " + str(count2))
            print("Noise objects: " + str(count3))

            if len(selects) > 1:
                report.extend(selects)
                report.append('</table>\n')

            if len(rejects) > 1:
                report.extend(rejects)
                report.append('</table>\n')

            if len(noise) > 1:
                report.extend(noise)
                report.append('</table>\n')

            report.append('</body>\n</html>')

        open(filename + '/html_report_' + image_name + '.html', 'w').writelines(report)

    def linescan_report(self, filename, image_name, linescan_manager):
        if len(linescan_manager.lines.keys()) > 0:
            HTML_HEADER = """<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN"
                            "http://www.w3.org/TR/html4/strict.dtd">
                        <html lang="en">
                          <head>
                            <meta http-equiv="content-type" content="text/html; charset=utf-8">
                            <title>title</title>
                            <link rel="stylesheet" type="text/css" href="style.css">
                            <script type="text/javascript" src="script.js"></script>
                          </head>
                          <body>\n"""

            report = [HTML_HEADER]
            table = '<table cellpadding="10" cellspacing="10">\n<th>Line ID</th><th>Images</th><th>Background</th><th>Membrane</th><th>Septum</th><th>Fluorescence Ratio</th>'

            for key in sorted(linescan_manager.lines.keys()):
                try:
                    lin = linescan_manager.lines[key]
                    img = linescan_manager.lines[key].image
                    w, h, dummy = img.shape
                    if w > h:
                        img = np.rot90(img)
                    imsave(filename + "/_linescan_images" +
                           os.sep + key + '.png', img)
                    row = '<tr style="text-align:center"><td>' + key + '</td><td><img src="./' + '_linescan_images/' + key + '.png" alt="pic" width="200"/></td>' + \
                        "<td>" + str(lin.background) + "</td>" + "<td>" + str(lin.membrane) + "</td>" + \
                        "<td>" + str(lin.septum) + "</td>" + \
                        "<td>" + str(lin.fr) + "</td></tr>"
                    table += row
                except IndexError:
                    print("One Line Not Saved: too close to edge of image")

            report += table
            report += "</table></body></html>"

            open(filename + '/linescan_report_' + image_name + '.html', 'w').writelines(report)

    def check_filename(self, filename):
        if os.path.exists(filename):
            tmp = ""
            split_path = filename.split("_")
            tmp = "_".join(split_path[:len(split_path)-1])
            tmp += "_" + str(int(split_path[-1])+1)
            return self.check_filename(tmp)

        else:
            return filename

    def generate_report(self, path, label, cell_manager, linescan_manager, params):
        if label is None:
            filename = path + "/Report_1"
            filename = self.check_filename(filename)
            self.cell_data_filename = filename

            if not os.path.exists(filename + "/_images"):
                os.makedirs(filename + "/_images")
            if not os.path.exists(filename + "/_rejected_images"):
                os.makedirs(filename + "/_rejected_images")
            if not os.path.exists(filename + "_noise_images"):
                os.makedirs(filename + "/_noise_images")
            if not os.path.exists(filename + "/_linescan_images"):
                os.makedirs(filename + "/_linescan_images")
        else:
            filename = path + "/Report_" + label + "_1"
            filename = self.check_filename(filename)
            self.cell_data_filename = filename

            if not os.path.exists(filename + "/_images"):
                os.makedirs(filename + "/_images")
            if not os.path.exists(filename + "/_rejected_images"):
                os.makedirs(filename + "/_rejected_images")
            if not os.path.exists(filename + "/_noise_images"):
                os.makedirs(filename + "/_noise_images")
            if not os.path.exists(filename + "/_linescan_images"):
                os.makedirs(filename + "/_linescan_images")

        selected_cells = ""
        for cell in cell_manager.cells.keys():
            if cell_manager.cells[cell].selection_state == CELL_SELECTED:
                selected_cells += cell + ";"

        self.csv_report(filename, label, cell_manager)
        self.html_report(filename, label, cell_manager)
        self.linescan_report(filename, label, linescan_manager)
        imsave(filename + "/selected_cells.png", cell_manager.fluor_w_cells)
        params.save_parameters(filename + "/params")
        open(filename + "/selected_cells.txt", "w").writelines(selected_cells)

    def get_cell_images(self, path, label, image_manager, cell_manager, params):
        if label is None:
            filename = self.cell_data_filename
            if not os.path.exists(filename + "/_cell_data/fluor"):
                os.makedirs(filename + "/_cell_data/fluor")

            if image_manager.optional_image is not None:
                if not os.path.exists(filename + "/_cell_data/optional"):
                    os.makedirs(filename + "/_cell_data/optional")
        else:
            filename = self.cell_data_filename
            if not os.path.exists(filename + "/_cell_data/fluor"):
                os.makedirs(filename + "/_cell_data/fluor")

            if image_manager.optional_image is not None:
                if not os.path.exists(filename + "/_cell_data/optional"):
                    os.makedirs(filename + "/_cell_data/optional")

        x_align, y_align = params.imageloaderparams.x_align, params.imageloaderparams.y_align

        fluor_img = image_manager.fluor_image
        fluor_w_cells = cell_manager.fluor_w_cells
        optional_image = image_manager.optional_image
        optional_w_cells = cell_manager.optional_w_cells

        for key in cell_manager.cells.keys():
            x0, y0, x1, y1 = cell_manager.cells[key].box
            fluor_cell = np.concatenate((fluor_img[x0:x1+1, y0:y1+1], fluor_img[x0:x1+1, y0:y1+1] * cell_manager.cells[key].cell_mask), axis=1)
            imsave(filename + "/_cell_data/fluor/" + key + ".png",
                   img_as_uint(fluor_cell))

            if optional_image is not None:
                optional_cell = np.concatenate((optional_image[x0:x1+1, y0:y1+1], optional_image[x0:x1+1, y0:y1+1] * cell_manager.cells[key].cell_mask), axis=1)
                imsave(filename + "/_cell_data/optional/" + key + ".png",
                       img_as_uint(optional_cell))
