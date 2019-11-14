'''
These tools are pulled *very* strongly from the great
ingredients in the `colour` package. They've been
modified and simplified here, in order to work better
with the design of rainbowconnection.
'''

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from six.moves import reduce

from colour.algebra import LinearInterpolator
from colour.colorimetry import (
    ILLUMINANTS, ILLUMINANTS_SDS, LIGHTNESS_METHODS, LUMINANCE_METHODS,
    SpectralShape, sd_blackbody, sd_ones, sd_to_XYZ,
    wavelength_to_XYZ)
from colour.plotting import (
    ColourSwatch, COLOUR_STYLE_CONSTANTS, XYZ_to_plotting_colourspace, artist,
    filter_passthrough, filter_cmfs, filter_illuminants, override_style,
    render, plot_single_colour_swatch, plot_multi_functions, plot_single_sd)
from colour.utilities import (domain_range_scale, first_item,
                              normalise_maximum, tstack)

from colour import SpectralDistribution

# define a standard set of color-matching functions
CMFs = first_item(filter_cmfs('CIE 1931 2 Degree Standard Observer').values())

def plot_rainbow(axes=None,
                 wavelength=None, flux=None,
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

    # pull out the CMFss
    cmfs = first_item(filter_cmfs(cmfs).values())

    if wavelength is None:
        # create a grid of wavelengths (at which CMFss are useful)
        wavelength = cmfs.wavelengths

    ok = ((wavelength >= np.min(cmfs.wavelengths)) &
          (wavelength <= np.max(cmfs.wavelengths)))

    # create colors at theose wavelengths
    colours = XYZ_to_plotting_colourspace(
        wavelength_to_XYZ(wavelength[ok], cmfs),
        ILLUMINANTS['CIE 1931 2 Degree Standard Observer']['E'],
        apply_encoding_cctf=False)

    # normalize the colors to their maximum?
    colours = COLOUR_STYLE_CONSTANTS.colour.colourspace.encoding_cctf(
        normalise_maximum(colours))

    # create y values that will be plotted (these could be spectrum)
    if flux is None:
        flux = np.ones_like(wavelength)
    x_min, x_max = min(wavelength), max(wavelength)
    y_min, y_max = 0, max(flux) + max(flux) * 0.05

    # create a polygon to define the top of the bars?
    polygon = Polygon(
        np.vstack([
            (x_min, 0),
            tstack([wavelength[ok], flux[ok]]),
            (x_max, 0),
        ]),
        facecolor='none',
        edgecolor='none')
    axes.add_patch(polygon)

    # draw bars, with the colors at each vertical stripe
    padding = 0.2
    axes.bar(
        x=wavelength[ok] - padding/2,
        height=max(flux[ok]),
        width=1 + padding,
        color=colours,
        align='edge')

    # plot the actual spectrum?
    # axes.plot(wavelength, values, color=COLOUR_STYLE_CONSTANTS.colour.dark)

    return axes

def rainbow_spectrum(axes=None,
                 wavelength=None, flux=None,
                 cmfs='CIE 1931 2 Degree Standard Observer',
                 rainbowtop=np.inf,
                 **kwargs):
    """
    (This is still *real* blarg-y*.)
    Plot a spectrum, with a rainbow underneath.

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

    if wavelength is None:
        # create a grid of wavelengths (at which CMFss are useful)
        wavelength = CMFs.wavelengths

    # create y values that will be plotted (these could be spectrum)
    if flux is None:
        flux = np.ones_like(wavelength)

    # pull out only the values that *can* be converted to colors
    ok = ((wavelength >= np.min(CMFs.wavelengths)) &
          (wavelength <= np.max(CMFs.wavelengths)))

    w, f = wavelength[ok], flux[ok]

    # clip the top of the box?
    f = np.minimum(f, rainbowtop)

    XYZ = wavelength_to_XYZ(w)

    # create colors at theose wavelengths
    colours = XYZ_to_plotting_colourspace(XYZ)

    # normalize the colors to their maximum?
    #colours = COLOUR_STYLE_CONSTANTS.colour.colourspace.encoding_cctf(
    #    normalise_maximum(colours))
    colours = np.maximum(0, colours/np.max(colours))

    x_min, x_max = min(w), max(w)
    y_min, y_max = 0, max(f) + max(f) * 0.05

    # create a polygon to define the top of the bars?
    polygon = Polygon(
        np.vstack([
            (x_min, 0),
            tstack([w, f]),
            (x_max, 0),
        ]),
        facecolor='none',
        edgecolor='none')
    axes.add_patch(polygon)

    # draw bars, with the colors at each vertical stripe
    padding = 0.0
    dw = np.mean(np.diff(w))
    axes.bar(
        x=w,
        height=f,
        width=dw,
        color=colours,
        clip_path=polygon,
        align='edge',
        clip_on=True)

    return axes
