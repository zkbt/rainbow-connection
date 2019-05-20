'''
These tools are pulled *very* strongly from the great
ingredients in the `colour` package. They've been
modified and simplified here, in order to work better
with the design of therainbowconnection.
'''

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from six.moves import reduce

from colour.algebra import LinearInterpolator
from colour.colorimetry import (
    ILLUMINANTS, ILLUMINANTS_SDS, LIGHTNESS_METHODS, LUMINANCE_METHODS,
    MultiSpectralDistribution, SpectralShape, sd_blackbody, sd_ones, sd_to_XYZ,
    wavelength_to_XYZ)
from colour.plotting import (
    ColourSwatch, COLOUR_STYLE_CONSTANTS, XYZ_to_plotting_colourspace, artist,
    filter_passthrough, filter_cmfs, filter_illuminants, override_style,
    render, plot_single_colour_swatch, plot_multi_functions)
from colour.utilities import (domain_range_scale, first_item,
                              normalise_maximum, tstack)


CMF = first_item(filter_cmfs('CIE 1931 2 Degree Standard Observer').values())

def plot_rainbow(axes=None,
                 cmfs='CIE 1931 2 Degree Standard Observer',
                 **kwargs):
    """
    Plot a simple horizontal rainbow,
    based off colour's plot_single_sd.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        The ax into which the rainbow should be drawn.
    cmfs : string
        The color matching function(s?) to use.

    Returns
    -------
    ax : matplotlib.axes._subplots.AxesSubplot
        The ax into which the rainbow was drawn.
    """

    if axes is None:
        axes = plt.gca()

    # pull out the CMFs
    cmfs = first_item(filter_cmfs(cmfs).values())

    # create a grid of wavelengths (at which CMFs are useful)
    wavelengths = cmfs.wavelengths


    # create colors at theose wavelengths
    colours = XYZ_to_plotting_colourspace(
        wavelength_to_XYZ(wavelengths, cmfs),
        ILLUMINANTS['CIE 1931 2 Degree Standard Observer']['E'],
        apply_encoding_cctf=False)

    # normalize the colors to their maximum?
    colours = COLOUR_STYLE_CONSTANTS.colour.colourspace.encoding_cctf(
        normalise_maximum(colours))

    # create y values that will be plotted (these could be spectrum)
    values = np.ones_like(wavelengths)
    x_min, x_max = min(wavelengths), max(wavelengths)
    y_min, y_max = 0, max(values) + max(values) * 0.05

    # create a polygon to define the top of the bars?
    polygon = Polygon(
        np.vstack([
            (x_min, 0),
            tstack([wavelengths, values]),
            (x_max, 0),
        ]),
        facecolor='none',
        edgecolor='none')
    axes.add_patch(polygon)

    # draw bars, with the colors at each vertical stripe
    padding = 0.1
    axes.bar(
        x=wavelengths - padding,
        height=max(values),
        width=1 + padding,
        color=colours,
        align='edge')

    # plot the actual spectrum?
    # axes.plot(wavelengths, values, color=COLOUR_STYLE_CONSTANTS.colour.dark)

    return axes
