"Module used to encapsulate some functions used in the cells module"

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
