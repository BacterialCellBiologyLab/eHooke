"""Module used to create the report of the cell identification"""
from skimage.io import imsave
from skimage.util import img_as_int, img_as_float
import cellprocessing as cp
import numpy as np
import os

CELL_SELECTED = 1
CELL_REJECTED = -1


class ReportManager:
    def __init__(self, parameters):
        self.keys = cp.stats_format(parameters.cellprocessingparams)

    def csv_report(self, filename, cell_manager):

        cells = cell_manager.cells

        if len(cells) > 0:
            header = 'Cell ID '

            for k in self.keys:
                label, digits = k
                header = header+';'+label
            selects = [header+'\n']
            rejects = [header+'\n']
            noise = [header+'\n']

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
                    selects.append(lin+'\n')

                elif cell.selection_state == CELL_REJECTED:
                    lin = str(int(cell.label))

                    for stat in self.keys:
                        lin = lin + ";" + str(cell.stats[stat[0]])
                    rejects.append(lin+"\n")

                elif cell.selection_state == 0:
                    lin = str(int(cell.label))

                    for stat in cell.stats:
                        lin = lin + ";" + str(cell.stats[stat[0]])
                    noise.append(lin+"\n")

        if len(selects) > 1:
            open(filename + "csv_selected.csv", 'w').writelines(selects)

        if len(rejects) > 1:
            open(filename + "csv_rejected.csv", "w").writelines(rejects)

        if len(noise) > 1:
            open(filename + "csv_noise.csv", "w").writelines(noise)

    def html_report(self, filename, cell_manager):
        """generates an html report with the all the cell stats from the selected cells"""

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
                header = header+'</th><th>'+label
            header += '</th>\n'
            selects = ['\n<h1>Selected cells:</h1>\n'+header+'\n']
            rejects = ['\n<h1>Rejected cells:</h1>\n'+header+'\n']
            noise = ['\n<h1>Noise:</h1>\n'+header+'\n']

            count = 0
            count2 = 0
            count3 = 0

            print "Total Cells: " + str(len(cells))

            sorted_keys = []
            for k in sorted(cells.keys()):
                sorted_keys.append(int(k))

            sorted_keys = sorted(sorted_keys)

            for k in sorted_keys:
                cell = cells[str(k)]
                if cell.selection_state == CELL_SELECTED:
                    cellid = str(int(cell.label))
                    img = img_as_float(cell.image)
                    imsave(filename+"/_images"+os.sep+cellid+'.png', img)
                    lin = '<tr><td>'+cellid+'</td><td><img src="./'+'_images/' + \
                          cellid+'.png" alt="pic" width="200"/></td>'

                    count += 1

                    for stat in self.keys:
                        lbl, digits = stat
                        lin = lin + '</td><td>' + ("{0:."+str(digits)+"f}").format(str(cell.stats[lbl]))

                    lin += '</td></tr>\n'
                    selects.append(lin)

                elif cell.selection_state == CELL_REJECTED:
                    cellid = str(int(cell.label))
                    img = img_as_float(cell.image)
                    imsave(filename+"/_rejected_images"+os.sep+cellid+'.png', img)
                    lin = '<tr><td>'+cellid+'</td><td><img src="./'+'_rejected_images/' + \
                          cellid+'.png" alt="pic" width="200"/></td>'

                    count2 += 1

                    for stat in self.keys:
                        lbl, digits = stat
                        lin = lin + '</td><td>' + ("{0:."+str(digits)+"f}").format(str(cell.stats[lbl]))

                    lin += '</td></tr>\n'
                    rejects.append(lin)

                elif cell.selection_state == 0:
                    cellid = str(int(cell.label))
                    img = img_as_float(cell.image)
                    imsave(filename+"/_noise_images"+os.sep+cellid+'.png', img)
                    lin = '<tr><td>'+cellid+'</td><td><img src="./'+'_noise_images/' + \
                          cellid+'.png" alt="pic" width="200"/></td>'

                    count3 += 1

                    for stat in self.keys:
                        lbl, digits = stat
                        lin = lin + '</td><td>' + ("{0:."+str(digits)+"f}").format(str(cell.stats[lbl]))

                    lin += '</td></tr>\n'
                    noise.append(lin)

            print "Selected Cells: " + str(count)
            print "Rejected Cells: " + str(count2)
            print "Noise objects: " + str(count3)

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

        open(filename+'html_report.html', 'w').writelines(report)

    def generate_report(self, path, label, cell_manager, params):
        if label is None:
            filename = path+"/Report/"
            if not os.path.exists(filename+"_images"):
                os.makedirs(filename+"/_images")
            if not os.path.exists(filename+"_rejected_images"):
                os.makedirs(filename+"/_rejected_images")
            if not os.path.exists(filename+"_noise_images"):
                os.makedirs(filename+"/_noise_images")
        else:
            filename = path+"/Report_"+label+"/"
            if not os.path.exists(filename+"_images"):
                os.makedirs(filename+"/_images")
            if not os.path.exists(filename+"_rejected_images"):
                os.makedirs(filename+"/_rejected_images")
            if not os.path.exists(filename+"_noise_images"):
                os.makedirs(filename+"/_noise_images")

        self.csv_report(filename, cell_manager)
        self.html_report(filename, cell_manager)
        cell_manager.fluor_with_cells_outlined.save_image(filename+"selected_cells")
        params.save_parameters(filename+"params")
