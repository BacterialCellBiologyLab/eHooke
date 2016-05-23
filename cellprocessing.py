"Module used to encapsulate some functions used in the cells module"

import numpy as np
from skimage import color
from skimage.util import img_as_int
from skimage.segmentation import mark_boundaries


def rotation_matrices(step):
    """ returns a list of rotation matrixes over 180 deg
    matrixes are transposed to use with 2 column point arrays (x,y),
    multiplying after the array
    TODO: optimize with np vectors
    """

    result = []
    ang = 0

    while ang < 180:
        sa = np.sin(ang / 180.0 * np.pi)
        ca = np.cos(ang / 180.0 * np.pi)
        # note .T, for column points
        result.append(np.matrix([[ca, -sa], [sa, ca]]).T)
        ang = ang + step

    return result


def bounded_value(minval, maxval, currval):
    """ returns the value or the extremes if outside
    """

    if currval < minval:
        return minval

    elif currval > maxval:
        return maxval

    else:
        return currval


def bounded_point(x0, x1, y0, y1, p):
    tx, ty = p
    tx = bounded_value(x0, x1, tx)
    ty = bounded_value(y0, y1, ty)
    return tx, ty


def bound_rectangle(points):
    """ returns a tuple (x0,y0,x1,y1,width) of the bounding rectangle
    points must be a N,2 array of x,y coords
    """

    x0, y0 = np.amin(points, axis=0)
    x1, y1 = np.amax(points, axis=0)
    a = np.min([(x1 - x0), (y1 - y0)])
    return x0, y0, x1, y1, a


def stats_format(params):
    """Returns the list of cell stats to be displayed on the report,
    depending on the computation of the septum"""
    result = []
    result.append(('Area', 0))
    result.append(('Perimeter', 0))
    result.append(('Length', 1))
    result.append(('Width', 1))
    result.append(('Eccentricity', 3))
    result.append(('Irregularity', 3))
    result.append(('Neighbours', 0))
    result.append(('Baseline', 4))
    result.append(('Cell Median', 4))
    result.append(('Membrane Median', 4))
    result.append(('Cytoplasm Median', 4))

    if params.find_septum:
        result.append(('Septum Median', 4))
        result.append(("Fluor Ratio", 4))
        result.append(("Fluor Ratio 75%", 4))
        result.append(("Fluor Ratio 25%", 4))
        result.append(("Fluor Ratio 10%", 4))

    return result


def overlay_cells(cells, image, colors):
    "Overlay the edges of each individual cell in the provided image"

    tmp = color.gray2rgb(image)

    for k in cells.keys():
        c = cells[k]
        if c.selection_state == 1:
            col = colors[c.color_i][:3]

            for px in c.outline:
                x, y = px
                tmp[x, y] = col

            if c.sept_mask is not None:
                try:
                    x0, y0, x1, y1 = c.box
                    tmp[x0:x1, y0:y1] = mark_boundaries(tmp[x0:x1, y0:y1],
                                                        img_as_int(
                                                            c.sept_mask),
                                                        color=col)
                except IndexError:
                    c.selection_state = -1

    return tmp


def assign_cell_color(cell, cells, cell_colors):
    """ assigns an index to cell.color that is different from the neighbours """

    neighcols = []

    for neigh in cell.neighbours.iterkeys():
        col = cells[str(int(neigh))].color_i
        if col not in neighcols:
            neighcols.append(col)

    cell.color_i = cell.stats["Area"] % len(
        cell_colors)  # each cell has a preferred color

    while len(neighcols) < len(cell_colors) and (cell.color_i in neighcols):
        cell.color_i += 1

        if cell.color_i >= len(cell_colors):
            cell.color_i = 0


def update_neighbours(cells, oldlabel, newlabel):
    """ updates the neighbour list when merging cells """
    oc = cells[str(oldlabel)]
    nc = cells[str(newlabel)]

    for nei in oc.neighbours.iterkeys():

        tc = cells[str(int(nei))]
        inter = tc.neighbours[oldlabel]
        del tc.neighbours[oldlabel]

        if int(nei) != newlabel:

            if newlabel in tc.neighbours:
                tc.neighbours[newlabel] = tc.neighbours[newlabel] + inter
                nc.neighbours[tc.label] = nc.neighbours[tc.label] + inter

            else:
                tc.neighbours[newlabel] = inter
                nc.neighbours[tc.label] = inter


def check_merge(cell1, cell2, rotations, interface, mask, params):

    if cell1.stats["Area"] <= 0 or cell2.stats["Area"] <= 0:  # check if both cells exist
        return False

    # check if any cell is small enough for automatic merge
    if cell1.stats["Area"] < params.cell_force_merge_below or \
            cell2.stats["Area"] < params.cell_force_merge_below:
        return True

    # check if dividing cells
    if params.merge_dividing_cells and \
            interface >= params.merge_min_interface:
        tmp = cells.Cell(0)
        tmp.outline.extend(cell1.outline)
        tmp.outline.extend(cell2.outline)
        tmp.compute_axes(rotations, mask.shape)
        tmpshort = axis_length(tmp.short)
        maxshort = max(axis_length(cell1.short), axis_length(cell2.short))

        if tmpshort <= maxshort * params.merge_length_tolerance:
            return True

        else:
            return False

    else:
        return False


def paint_cell(cell, image, newval):
    """ paints the lines of the cell into the image """

    for li in cell.lines:
        y, x0, x1 = li
        image[x0:x1 + 1, y] = newval

    return image


def blocked_by_filter(cell, list_of_filters):
    """ returns true if cell is blocked by any filter
        [("stat", min, max), ("stat2", min, max)]
    """
    for filt in list_of_filters:
        val = cell.stats[filt[0]]
        if (val < filt[1]) or (val > filt[2]):
            return True

    return False
